#!/bin/bash
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10

wget -c http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004131/harmonised/28067908-GCST004131-EFO_0003767.h.tsv.gz -O IBD_hg38_harmonized.tsv.gz
