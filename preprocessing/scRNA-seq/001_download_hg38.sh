#!/bin/bash -l
#SBATCH --job-name=download     # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=marliette.rodriguezmatos@einsteinmed.edu
#SBATCH -p pe2 
#SBATCH -t 7-00:00:00
#SBATCH --cpus-per-task=4                   
#SBATCH --mem=32G 
#SBATCH --output=download-%j.out

## downloading latest gencode human reference genome 38
echo "Downloading fasta"
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz

echo "Downloading gft annotation"
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz

echo "Downloading assembly report"
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_assembly_report.txt

echo "Downloading annotation for pseudo autosomal region PAR in ChrY"
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_assembly_structure/Primary_Assembly/pseudoautosomal_region/par_align.gff


echo "Decompressing files"
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v44.primary_assembly.annotation.gtf.gz


