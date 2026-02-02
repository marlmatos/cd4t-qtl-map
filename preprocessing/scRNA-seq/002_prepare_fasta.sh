#!/bin/bash -l
#SBATCH --job-name=prepare_fasta    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=marliette.rodriguezmatos@einsteinmed.edu
#SBATCH -p pe2 
#SBATCH -t 7-00:00:00
#SBATCH --cpus-per-task=4                   
#SBATCH --mem=32G 
#SBATCH --output=002_prepare_fasta-%j.out


###########################################################
##                                                       ##
##   This script is for preparing the FASTA file         ##
##    of hg38 for transcriptomic alignment               ##
##    removing alternative contigs and patches           ##
##   not integrated into the primary build               ##
###   Adapted from advice found                          ##
##    "https://www.biostars.org/p/342482/"               ##
## https://manual.omicsbox.biobam.com/user-manual/omicsbox-modules/ module-transcriptomics/rna-seq-alignment/rna-seq-star/            ##
##                                                       ##
###########################################################

root="/gchm/resources/genome/hg38_gencode_raw"
wdir="/gchm/resources/genome/hg38_gencode_PRI_align"

module load samtools/1.15

echo "indexing fasta"
samtools faidx $root/GRCh38.primary_assembly.genome.fa

#It is strongly recommended to include major chromosomes as well as un-placed and un-localized scaffolds since a substantial number of reads may map to these scaffolds (e.g. ribosomal RNA). 

echo "subsetting annotation for primary and mitochondrial contigs"
sort -k1,1V $root/GCA_000001405.29_GRCh38.p14_assembly_report.txt | awk -v FS="\t" '$2 == "assembled-molecule" || $8 == "non-nuclear" {print $10}' > $wdir/test1

echo "subsetting annotation for unlocalized fragments"
sort -k1,1V $root/GCA_000001405.29_GRCh38.p14_assembly_report.txt | awk -v FS="\t" '$2 == "unlocalized-scaffold" || $2 == "unplaced-scaffold" {print $5}' > $wdir/test2

echo "merging annotation"

cat $wdir/test1 $wdir/test2 > $wdir/subset_ids.txt

rm test1 test2
#command to convert the file format from DOS/Windows to Unix/Linux
#dos2unix $root/GCA_000001405.29_GRCh38.p14_assembly_report.txt

echo "Subset Fasta"
samtools faidx $root/GRCh38.primary_assembly.genome.fa -r $wdir/subset_ids.txt -o $wdir/GRCh38.primary_assembly_subset.genome.fa

echo "Masking the PAR on chromosome Y"
sed -E 's/.*Target=([^;]+).*/\1/g' $root/par_align.gff| awk -v OFS="\t" '$0 !~ "^#" {print $1, $2-1, $3}'  > $wdir/parY.bed

bedtools maskfasta -fi $wdir/GRCh38.primary_assembly_subset.genome.fa -bed $wdir/parY.bed -fo $wdir/GRCh38.primary_assembly_subset_masked.genome.fa

echo "indexing alignment ready fasta"
samtools faidx $wdir/GRCh38.primary_assembly_subset_masked.genome.fa

grep -e ">" GRCh38.primary_assembly_subset_masked.genome.fa