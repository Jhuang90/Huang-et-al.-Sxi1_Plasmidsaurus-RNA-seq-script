# Huang-et-al.-Sxi1_Plasmidsaurus-RNA-seq-script
This repository contains scripts used to analyze Plasmidsaurus 3′ end RNA sequencing data (https://plasmidsaurus.com/rna
). These scripts were applied to investigate transcriptional differences in Ménage à trois mating experiments, as described in our preprint: 
https://www.biorxiv.org/content/10.1101/2025.02.11.637763v3 and formal publication https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1012084

The analysis workflow implemented here was informed by and adapted from the Plasmidsaurus recommended RNA sequencing analysis guidelines (https://plasmidsaurus.com/technical-documentation/rna
).

### The main analysis steps include

1. Read filtering using FastP v0.24.0
2. Alignment to the appropriate reference genome using STAR v2.7.11
3. Coordinate sorting of BAM files using samtools v1.22.1
4. UMI based de duplication
5. Mapping quality control
6. Gene expression quantification using featureCounts from the Subread package v2.1.1
7. Principal component analysis
8. Differential gene expression analysis

