"""
Movi analysis
"""

#Requirements:
#In data directory: genome assemblies in fasta format
#In config directory: assemblylist.txt and gff_table.txt

from pathlib import Path
with open("config/assemblylist.txt", "r", encoding="utf-8") as assemblylist:
    samplenames = [Path(filename.strip()).stem for filename in assemblylist]

configfile: "config/config.yaml" 

rule all:
    input:
        expand("results/annotated/{sample}.gff", sample=samplenames),
        "results/assess/out_checkm.tab",
        "results/out_fastani.tsv.matrix",
        "sidekick.done",
        "results/tree/T3.raxml.support",
        "results/tree/Tree_auspice.json"

rule assess:
    input:
        datafolder = "data"
    output:
        checkmout = "results/assess/out_checkm.tab"
    log:
        "logs/checkm.txt"
    conda:
        "envs/checkm.yaml"
    shell:
        """
        checkm lineage_wf \
        {input.datafolder} \
        results/assess \
        --tab_table \
        -f {output.checkmout} &> {log}
        """

rule ani:
    input:
        assemblylist = "config/assemblylist.txt"
    output:
        ani_values = "results/out_fastani.tsv.matrix"
    log:
        "logs/ani.txt"
    conda:
        "envs/fastani.yaml"
    shell:
        """
        fastANI \
        --ql {input.assemblylist} \
        --rl {input.assemblylist} \
        -o results/out_fastani.tsv \
        --matrix  &> {log}
        """   

rule annotate:
    input:
        assemblies = "data/{sample}.fna" 
    output:
        annotated = "results/annotated/{sample}.gff" 
    log:
        "logs/{sample}_annotate.txt"
    conda:
        "envs/prokka.yaml"
    params:
        genetic_code = config["genetic_code"]
    shell:
        """
        prokka \
        --outdir results/annotated \
        --kingdom Bacteria \
        --gcode {params.genetic_code} \
        --force \
        --prefix {wildcards.sample} \
        {input.assemblies}  &> {log} 
        """

rule pangenome:
    input: 
        annotated = expand("results/annotated/{sample}.gff", sample=samplenames)
    output: 
        alignment = "results/pangenome/core_gene_alignment.aln"
    log:
        "logs/pangenome.txt"
    conda:
        "envs/roary.yaml"
    params:
        genetic_code = config["genetic_code"],
        percent_identity = config["percent_identity"]
    threads: 4
    shell:
        """
        rm -r results/pangenome; \
        roary \
        -f results/pangenome  \
        -t {params.genetic_code} \
        -i {params.percent_identity} \
        -e -z -v \
        -p {threads} \
        {input.annotated}  &> {log} 
        """

rule sidekick:
    input:
        datafolder = "results/pangenome/pan_genome_sequences"
    output:
        touch("sidekick.done")
    log:
        "logs/sidekick.txt"
    conda:
        "envs/sidekick.yaml"
    shell:
        """
        python3 config/sidekick.py \
        config/gff_table.txt \
        --alns {input.datafolder} \
        --output results/pangenome/sidekickgenes  &> {log}
        """

rule tree_ml:
    input:
        alignment = "results/pangenome/core_gene_alignment.aln"
    output:
        besttree = "results/tree/T1.raxml.bestTree" 
    log:
        "logs/tree_ml.txt"
    conda:
        "envs/raxml.yaml"
    params:
       model = "GTR+G",
       prefix = "results/tree/T1",
       seed = "2"
    threads: 4
    shell:
        """
        raxml-ng \
        --msa {input.alignment} \
        --model {params.model} \
        --prefix {params.prefix} \
        --threads {threads} \
        --seed {params.seed} \
        --tree pars{{25}},rand{{25}} &> {log}
        """

rule tree_bootstrap:
    input:
        alignment = "results/pangenome/core_gene_alignment.aln"
    output:
        bootstraps = "results/tree/T2.raxml.bootstraps" 
    log:
        "logs/tree_bootstrap.txt"
    conda:
        "envs/raxml.yaml"
    params:
        model = "GTR+G",
        prefix = "results/tree/T2",
        seed = "2"
    threads: 4
    shell:
        """
        raxml-ng --bootstrap \
        --msa results/pangenome/core_gene_alignment.aln \
        --model {params.model} \
        --prefix {params.prefix} \
        --seed {params.seed} \
        --threads {threads}  &> {log}
        """

rule tree_support:
    input:
        besttree = "results/tree/T1.raxml.bestTree",
        bootstraps = "results/tree/T2.raxml.bootstraps" 
    output:
        bootstraptree = "results/tree/T3.raxml.support" 
    log:
        "logs/tree_support.txt"
    conda:
        "envs/raxml.yaml"
    params:
        prefix = "results/tree/T3"
    threads: 4
    shell:
        """
        raxml-ng \
        --support \
        --tree {input.besttree} \
        --bs-trees {input.bootstraps} \
        --prefix {params.prefix} \
        --threads {threads} &> {log}
        """

rule nextstrain_refine:
    input:
        besttree = "results/tree/T1.raxml.bestTree",
        alignment = "results/pangenome/core_gene_alignment.aln"
    output:
        augurtree = "results/tree/T1.raxml.bestTree_augur.nwk",
        nodes = "results/tree/Tree_NodeData.json"
    log:
        "logs/nextstrain_refine.txt"
    conda:
        "envs/nextstrain.yaml"
    threads: 4
    shell:
        """
        augur refine \
        --tree {input.besttree} \
        --alignment {input.alignment} \
        --output-tree {output.augurtree} \
        --output-node-data {output.nodes} \
        --root mid_point &> {log}; 
        """

rule nextstrain_export:
    input:
        augurtree = "results/tree/T1.raxml.bestTree_augur.nwk",
        nodes = "results/tree/Tree_NodeData.json"
    output:
        auspicetree = "results/tree/Tree_auspice.json"
    log:
        "logs/nextstrain_export.txt"
    conda:
        "envs/nextstrain.yaml"
    threads: 4
    shell:
        """
        augur export v2 \
        --tree {input.augurtree}  \
        --node-data {input.nodes}  \
        --output {output.auspicetree}  &> {log}; 
        """

