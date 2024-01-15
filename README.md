## Pipeline for comparative genomics of *Mycoplasma ovipneumoniae*

### Analysis summary
* Genome assembly:
    * Clean raw Illumina sequence data using [HTStream](https://github.com/s4hts/HTStream)
    * Detect any potential contamination using [Centrifuge](https://github.com/DaehwanKimLab/centrifuge) 
    * [SPAdes](https://github.com/ablab/spades) assembly for samples with only Illumina data
    * [Trycycler](https://github.com/rrwick/Trycycler) assembly for samples with Oxford Nanopore sequence data
* Snakemake workflow:
    * Evaluate assembly completeness and contamination using [CheckM](https://github.com/Ecogenomics/CheckM)
    * Calculate Average Nucleotide Identity using [FastANI](https://github.com/ParBLiSS/FastANI)
    * Annotate genomes using [Prokka](https://github.com/tseemann/prokka)
    * Pangenome analysis using [Roary](https://github.com/sanger-pathogens/Roary)
    * Add correct sample names to the Roary gene alignments for [fastGEAR](https://academic.oup.com/mbe/article/34/5/1167/2983515?login=true) analysis
    * Phylogenetic analysis of core genome using [RAxML](https://github.com/amkozlov/raxml-ng)
    * Create interactive phylogeny using [Nextstrain](https://github.com/nextstrain/augur)
* Other analyses:
    * Annotate Roary output file pan_genome_reference.fa using [EggNOG-mapper with Diamond reference database](https://github.com/eggnogdb/eggnog-mapper)
    * Find core genes under positive or diversifying selection using [PAML](https://github.com/abacus-gene/paml)
    * PanGWAS using [Scoary](https://github.com/AdmiralenOla/Scoary)
    * Find CRISPR regions using [CrisprCASFinder](https://github.com/dcouvin/CRISPRCasFinder)

### Quick Start for Snakemake workflow

#### Requirements

* Genome assemblies should be moved into "data" directory: Each assembly should be a separate file in fasta format
* assemblylist.txt in "config" directory: Text file containing one column with relative paths to each assembly
* gff_table.txt in "config" directory: Text file containing 2 columns, one with the relative paths to each gff file that will be created by Prokka, and one with sample names

#### The example snakemake workflow provided here can be run after downloading and decompressing *M. ovipneumoniae* genome assemblies to the "data" directory using the following commands:
```
curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/714/895/GCA_001714895.1_ASM171489v1/GCA_001714895.1_ASM171489v1_genomic.fna.gz | gunzip -c > ./data/GCA_001714895.1_ASM171489v1_genomic.fna

curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/218/525/GCA_000218525.2_ASM21852v2/GCA_000218525.2_ASM21852v2_genomic.fna.gz  | gunzip -c > ./data/GCA_000218525.2_ASM21852v2_genomic.fna

curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/590/415/GCA_000590415.1_Version_1.0/GCA_000590415.1_Version_1.0_genomic.fna.gz | gunzip -c > ./data/GCA_000590415.1_Version_1.0_genomic.fna

curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/715/025/GCA_001715025.1_ASM171502v1/GCA_001715025.1_ASM171502v1_genomic.fna.gz | gunzip -c > ./data/GCA_001715025.1_ASM171502v1_genomic.fna

```
#### Running the Snakemake workflow
```
snakemake --cores 4 --use-conda --conda-frontend conda
```

### Genome assembly commands

#### Clean raw Illumina sequence data using [HTStream](https://github.com/s4hts/HTStream)

Example command:
```
hts_Stats -L ./01-cleaned/Sample01_stats.log -1 ./00-RawData/Sample01_R1.fastq.gz -2 ./00-RawData/Sample01_R2.fastq.gz | hts_SeqScreener -k 12 -AL ./01-cleaned/Sample01_stats.log | hts_SuperDeduper -e 100000 -AL ./01-cleaned/Sample01_stats.log | hts_AdapterTrimmer -m 250 -AL ./01-cleaned/Sample01_stats.log | hts_NTrimmer -n -m 250 -AL ./01-cleaned/Sample01_stats.log | hts_Stats -AL ./01-cleaned/Sample01_stats.log -fgp ./01-cleaned/Sample01
```

#### Detect any potential contamination using [Centrifuge](https://github.com/DaehwanKimLab/centrifuge) 
Example command:
```
centrifuge -p 40 -t -x ./databases/centrifuge/nt -1 ./01-cleaned/Sample01_R1.fastq.gz -2 ./01-cleaned/Sample01_R2.fastq.gz -S ./02-Centrifuge/Sample01_results.tsv --report-file ./02-Centrifuge/Sample01_report.tsv

centrifuge-kreport -x ./databases/centrifuge/nt ./02-Centrifuge/Sample01_results.tsv > ./02-Centrifuge/Sample01_results.Kreport
```

#### [SPAdes](https://github.com/ablab/spades) assembly for samples with only Illumina data
Example command:
```
spades.py --careful -t 60 -k 33,77,127 -1 ./01-cleaned/Sample01_R1.fastq.gz -2 ./01-cleaned/Sample01_R2.fastq.gz -o ./02-assembled/Sample01
```
#### [Trycycler](https://github.com/rrwick/Trycycler) assembly for samples with Oxford Nanopore sequence data

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

### Other analysis commands (after Snakemake workflow)

#### Annotate Roary output file pan_genome_reference.fa using [EggNOG-mapper with Diamond reference database](https://github.com/eggnogdb/eggnog-mapper)
```
emapper.py -i ./04-Roary/Prank_all/pan_genome_reference.fa --itype CDS --data_dir ./ref_eggnog2 --trans_table 4 --output_dir ./04-Roary/Prank_all/eggnog_annotations -o pan_genome_eggnog
```

#### Find core genes under positive or diversifying selection using [PAML](https://github.com/abacus-gene/paml)
```
codeml Movi_codeml.ctl
```

#### PanGWAS using [Scoary](https://github.com/AdmiralenOla/Scoary)
```
scoary -g RoaryGenes_forScoary.csv -t ScoaryTraits.csv
```

#### Find CRISPR regions using [CrisprCASFinder](https://github.com/dcouvin/CRISPRCasFinder)
Example command:
```
perl CRISPRCasFinder.pl -in ./02-assembled/Sample01.fna  -cas -keep -log -gcode 4
```
