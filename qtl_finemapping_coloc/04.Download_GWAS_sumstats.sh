#!/bin/bash
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8

mkdir -p GWAS_sumstats

base_url="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"

# Telomere length
#wget -P GWAS_sumstats "${base_url}GCST90435001-GCST90436000/GCST90435144/harmonised/GCST90435144.h.tsv.gz"

# COVID GWAS
wget -P GWAS_sumstats "${base_url}GCST90134001-GCST90135000/GCST90134601/harmonised/GCST90134601.h.tsv.gz"

# Load immune GWAS lists
while read line; do

	STUDY=$(echo $line | awk '{print$1}')
	# Remove the GCST part and strip leading zeros
        num=$(echo "${STUDY/GCST/}" | sed 's/^0*//')	

	# Round down to the nearest 1000 and plus 1
	lower=$((num / 1000 * 1000 + 1))
	
	# Round up to the nearest 1000
	upper=$(( (num + 999) / 1000 * 1000 ))


        full_url="${base_url}GCST$(printf "%08d" $lower)-GCST$(printf "%08d" $upper)/${STUDY}/harmonised/${STUDY}.h.tsv.gz"

	echo "Downloading ${full_url}"
	
	wget -P GWAS_sumstats ${full_url}
	
	echo "Done!"

done < <(tail -n +2 GWAS_sumstats/autoimmune_GWAS_list.tsv)
