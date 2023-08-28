## Pipeline for comparative genomics of *Mycoplasma ovipneumoniae*

### Clean raw Illumina sequence data using HTStream

Example command:
```
hts_Stats -L ./01-cleaned/Sample01_stats.log -1 ./00-RawData/Sample01_R1.fastq.gz -2 ./00-RawData/Sample01_R2.fastq.gz | hts_SeqScreener -k 12 -AL ./01-cleaned/Sample01_stats.log | hts_SuperDeduper -e 100000 -AL ./01-cleaned/Sample01_stats.log | hts_AdapterTrimmer -m 250 -AL ./01-cleaned/Sample01_stats.log | hts_NTrimmer -n -m 250 -AL ./01-cleaned/Sample01_stats.log | hts_Stats -AL ./01-cleaned/Sample01_stats.log -fgp ./01-cleaned/Sample01
```

### Detect any potential contamination using Centrifuge 
Example command:
```
centrifuge -p 40 -t -x ./databases/centrifuge/nt -1 ./01-cleaned/Sample01_R1.fastq.gz -2 ./01-cleaned/Sample01_R2.fastq.gz -S ./02-Centrifuge/Sample01_results.tsv --report-file ./02-Centrifuge/Sample01_report.tsv

centrifuge-kreport -x ./databases/centrifuge/nt ./02-Centrifuge/Sample01_results.tsv > ./02-Centrifuge/Sample01_results.Kreport
```

### Illumina assembly using SPAdes (for samples with no Oxford Nanopore data)
Example command:
```
spades.py --careful -t 60 -k 33,77,127 -1 ./01-cleaned/Sample01_R1.fastq.gz -2 ./01-cleaned/Sample01_R2.fastq.gz -o ./02-assembled/Sample01
```
### Trycycler assembly for Oxford Nanopore sequence data

#### 1. Remove short and low-quality reads
Example command:
```
/bin/Filtlong/bin/filtlong --min_length 500 --keep_percent 95 ./00-RawData_merged/Sample01.fastq > ./02-assembly_Trycycler/01-filter/Sample01_filtered.fastq
```
#### 2. Subsample the reads
Example command:
```
trycycler subsample --reads ./02-assembly_Trycycler/01-filter/Sample01_filtered.fastq --out_dir ./02-assembly_Trycycler/02-subsamples/Sample01 --genome_size 1m --count 18
```
#### 3. Flye assemblies
Example command:
```
flye --nano-raw ./02-assembly_Trycycler/02-subsamples/Sample01/subsample_01.fastq --threads 16 --plasmids --out-dir ./02-assembly_Trycycler/assembly_01_Sample01 && cp ./02-assembly_Trycycler/assembly_01_Sample01/assembly.fasta ./02-assembly_Trycycler/03-assemblies/Sample01/assembly_01.fasta && cp ./02-assembly_Trycycler/assembly_01_Sample01/assembly_graph.gfa ./02-assembly_Trycycler/03-assemblies/Sample01/assembly_01.gfa && rm -r ./02-assembly_Trycycler/assembly_01_Sample01
```
#### 4. Miniasm assemblies
Example command:
```
/mnt/lfs2/kandrews/bin/miniasm_and_minipolish.sh ./02-assembly_Trycycler/02-subsamples/Sample01/subsample_07.fastq 16 > ./02-assembly_Trycycler/assembly_07_Sample01.gfa && /mnt/lfs2/kandrews/bin/any2fasta ./02-assembly_Trycycler/assembly_07_Sample01.gfa > ./02-assembly_Trycycler/03-assemblies/Sample01/assembly_07.fasta && cp ./02-assembly_Trycycler/assembly_07_Sample01.gfa ./02-assembly_Trycycler/03-assemblies/Sample01/assembly_07.gfa && rm ./02-assembly_Trycycler/assembly_07_Sample01.gfa
```
#### 5. Raven assemblies
Example command:
```
raven --threads 16 ./02-assembly_Trycycler/02-subsamples/Sample01/subsample_13.fastq > ./02-assembly_Trycycler/03-assemblies/Sample01/assembly_13.fasta && rm raven.cereal
```
#### 6. Clustering
Example command:
```
trycycler cluster --assemblies ./02-assembly_Trycycler/03-assemblies/Sample01/*.fasta --reads ./02-assembly_Trycycler/01-filter/Sample01_filtered.fastq --out_dir ./02-assembly_Trycycler/04-clustering/Sample01
```
#### 7. Reconciliation
Example command:
```
trycycler reconcile --reads ./02-assembly_Trycycler/01-filter/Sample01_filtered.fastq --cluster_dir ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001
```
#### 8. Multiple sequence alignment
Example command:
```
trycycler msa --cluster_dir ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001
```
#### 9. Partition reads
Example command:
```
trycycler partition --reads ./02-assembly_Trycycler/01-filter/Sample01_filtered.fastq --cluster_dirs ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001
```
#### 10. Generate a consensus
Example command:
```
trycycler consensus --cluster_dir ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001
```
#### 11. Polish with Oxford Nanopore reads using Medaka
Example command:
```
medaka_consensus -i ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001/4_reads.fastq -d ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001/7_final_consensus.fasta -o ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001/medaka -m r941_min_high_g360 && mv ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001/medaka/consensus.fasta ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001/8_medaka.fasta && rm -r ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001/medaka ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001/*.fai ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001/*.mmi
```
#### 12. Polish with Illumina reads
Example command:
```
bowtie2 -1 ./00-RawData_Illumina/Sample01_R1.fastq.gz -2 ./00-RawData_Illumina/Sample01_R2.fastq.gz -x ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001/8_medaka  --threads 16 -I 305 -X 912 --local --very-sensitive-local | samtools sort > ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001/illumina_alignments.bam && samtools index ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001/illumina_alignments.bam && pilon --genome ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001/8_medaka.fasta --frags ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001/illumina_alignments.bam --output ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001/8_medaka_pilon1 --changes && rm ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001/illumina_alignments.bam && rm ./02-assembly_Trycycler/04-clustering/Sample01/cluster_001/illumina_alignments.bam.bai
```

### Evaluate assembly completeness and contamination using CheckM

```
checkm lineage_wf -t 8 ./02-assembled ./02-assembly_checkm
```
### Annotate genomes using Prokka
Example command
```
prokka --outdir ./03-annotated --kingdom Bacteria --gcode 4 --prefix Sample01 ./00-assemblies_all/Sample01.fna
```
### Pangenome analysis using Roary
```
roary -f ./04-Roary/Prank_all -e -z -v -t 4 ./03-annotated/*gff -p 24 -i 92
```
### Annotate the Roary output file pan_genome_reference.fa using EggNOG-mapper 
```
emapper.py -i ./04-Roary/Prank_all/pan_genome_reference.fa --itype CDS --data_dir ./ref_eggnog2 --trans_table 4 --output_dir ./04-Roary/Prank_all/eggnog_annotations -o pan_genome_eggnog
```

### PanGWAS using Scoary
```
scoary -g RoaryGenes_forScoary.csv -t ScoaryTraits.csv
```
### Phylogenetic analysis using RAxML
```
raxml-ng --msa ./04-Roary/Prank_all/core_gene_alignment.aln --model GTR+G --prefix T1 --threads 20 --seed 2 --tree pars{25},rand{25}

raxml-ng --bootstrap --msa ./04-Roary/Prank_all/core_gene_alignment.aln --model GTR+G --prefix T2 --seed 2 --threads 20

raxml-ng --support --tree T1.raxml.bestTree --bs-trees T2.raxml.bootstraps --prefix T3 --threads 2 
```

### Calculate Average Nucleotide Identity with FastANI
```
fastANI --ql AssemblyList.txt --rl AssemblyList.txt -o ./05-FastANI/AllAssemblies/FastANIout_AllAssemblies.tsv --matrix
```

### Find core genes under positive or diversifying selection with PAML
```
codeml Movi_codeml.ctl
```

### Add the correct sample names to the Roary gene alignments for fastGEAR analysis
```
python3 sidekick.py ./06-fastGEAR/Prank_all/gff_table.txt --alns ./04-Roary/Prank_all/pan_genome_sequences --output ./06-fastGEAR/Prank_all/sidekick_genes 
```
### Find CRISPR regions with CrisprCASFinder
Example command:
```
perl CRISPRCasFinder.pl -in ./02-assembled/Sample01.fna  -cas -keep -log -gcode 4
```
### Create interactive phylogeny with Nextstrain
```
augur refine --tree ./04-Roary/Prank_all/T3.raxml.bestTree --alignment ./04-Roary/Prank_all/core_gene_alignment.aln --output-tree ./07-Nextstrain/T3.raxml.bestTree_forAuspice.nwk --output-node-data  ./07-Nextstrain/Tree_NodeData.json --root mid_point

augur export v2 --tree ./07-Nextstrain/T3.raxml.bestTree_forAuspice.nwk  --node-data ./07-Nextstrain/Tree_NodeData.json  --output ./07-Nextstrain/Movi_phylogeny.json --metadata ./07-Nextstrain/MoviMetadataForAuspice.txt --auspice-config ./07-Nextstrain/Movi_auspice_config.json --colors metadata_colors.tsv 
```