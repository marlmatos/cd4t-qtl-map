#!/bin/bash
#SBATCH --job-name=CalculateLD                 # Job name
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

CHR=$1
LOWER=$2
UPPER=$3
RANGE=${LOWER}.${UPPER}

DIR="1kg_LD/chr${CHR}/${RANGE}"

mkdir -p $DIR

echo "Calculating LD for region Chr ${CHR}:${RANGE}"
/path/to/software/plink/plink-2.0a5.10/plink2_64bit --vcf "/path/to/lab_group/data/1kg/phase3_GRCh38/ALL.chr${CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz" \
--keep 1kg_eu_samples.txt \
--make-bed \
--set-all-var-ids @:#[b38]\$r,\$a \
--new-id-max-allele-len 20 missing \
--extract "${DIR}/variants_to_extract.txt" \
--out "${DIR}/${RANGE}"

/path/to/software/plink/plink-1.90b6.24/plink --bfile "${DIR}/${RANGE}" \
--keep-allele-order \
--r square \
--out "${DIR}/${RANGE}"

echo "We made it!"
