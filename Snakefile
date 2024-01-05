"""
Movi analysis
"""

#Requirements:
#In config folder: assemblylist.txt and gff_table.txt
#In data folder: genome assemblies in fasta format

from pathlib import Path
assemblylist = open("config/assemblylist.txt")
lines=assemblylist.readlines()
filenames=[line.strip() for line in lines]
samplenames = [Path(file).stem for file in filenames]

configfile: "config/config.yaml" 

rule all:
    input:
        expand("results/annotated/{sample}.gff", sample=samplenames),
        "results/assess/out_checkm.tab",
        "results/out_fastani.tsv.matrix",
        "results/pangenome/core_gene_alignment.aln",
        "results/T1.raxml.bestTree",
        "results/T2.raxml.bootstraps",
        "results/T3.raxml.support",
        "results/T1.raxml.bestTree_augur.nwk"


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

checkpoint pangenome:
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
        percent_identity = config["percent_identity"],
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

def sidekick_input(wildcards):
    checkpoint_output = checkpoints.pangenome.get(**wildcards).output[0]
    return expand("results/pangenome/pan_genome_sequences/{i}.txt",
                  i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fa.aln")).i)

def sidekick_output(wildcards):
    checkpoint_output = checkpoints.pangenome.get(**wildcards).output[0]
    return expand("results/pangenome/sidekickgenes/{i}.txt",
                  i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fa.aln")).i)

rule sidekick:
    input:
        sidekick_input
#    output:
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

rule tree1:
    input:
        alignment = "results/pangenome/core_gene_alignment.aln"
    output:
        besttree = "results/T1.raxml.bestTree" 
    log:
        "logs/tree1.txt"
    conda:
        "envs/raxml.yaml"
    params:
       model = "GTR+G",
       prefix = "T1",
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
        --tree pars{{25}},rand{{25}} &> {log}; \
        mv T1* results 
        """

rule tree2:
    input:
        alignment = "results/pangenome/core_gene_alignment.aln"
    output:
        bootstraps = "results/T2.raxml.bootstraps" 
    log:
        "logs/tree2.txt"
    conda:
        "envs/raxml.yaml"
    threads: 4
    shell:
        """
        raxml-ng --bootstrap \
        --msa results/pangenome/core_gene_alignment.aln \
        --model GTR+G \
        --prefix T2 \
        --seed 2 \
        --threads {threads}  &> {log}; \
        mv T2* results  
        """

rule tree3:
    input:
        besttree = "results/T1.raxml.bestTree",
        bootstraps = "results/T2.raxml.bootstraps" 
    output:
        bootstraptree = "results/T3.raxml.support" 
    log:
        "logs/tree3.txt"
    conda:
        "envs/raxml.yaml"
    threads: 4
    shell:
        """
        raxml-ng \
        --support \
        --tree {input.besttree} \
        --bs-trees {input.bootstraps} \
        --prefix T3 \
        --threads {threads} &> {log}; \
        mv T3* results  
        """

rule nextstrain1:
    input:
        besttree = "results/T1.raxml.bestTree",
        alignment = "results/pangenome/core_gene_alignment.aln"
    output:
        augurtree = "results/T1.raxml.bestTree_augur.nwk",
        nodes = "results/Tree_NodeData.json"
    log:
        "logs/nextstrain1.txt"
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



