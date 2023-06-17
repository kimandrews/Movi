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

### Illumina assembly using SPAdes 
Example command:
```
spades.py --careful -t 60 -k 33,77,127 -1 ./01-cleaned/Sample01_R1.fastq.gz -2 ./01-cleaned/Sample01_R2.fastq.gz -o ./02-assembled/Sample01
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
raxml-ng --msa core_gene_alignment.aln --model GTR+G --prefix T1 --threads 20 --seed 2 --tree pars{25},rand{25}

raxml-ng --bootstrap --msa core_gene_alignment.aln --model GTR+G --prefix T2 --seed 2 --threads 20

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

