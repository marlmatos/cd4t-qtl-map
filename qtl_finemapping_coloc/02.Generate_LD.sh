#!/bin/bash
#SBATCH --job-name=CalculateLD                 # Job name
#SBATCH --mem=80G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=7-00:00:00     # 7 days

module purge
module load plink/1.9

REGION_FILE=$1
LD_FILE=$2
NAME=$3

echo "Extracting variants and calculating LD with 02.Generate_LD.sh file"

tail -n +2 "$REGION_FILE" | while read line
do

	CHR=$(echo $line | awk '{print$1}')
	LOWER=$(echo $line | awk '{print$2}')
	UPPER=$(echo $line | awk '{print$3}')
	RANGE=${LOWER}-${UPPER}

	DIR="${NAME}_LD_matrices/chr${CHR}/${RANGE}"
		
	mkdir -p ${DIR}
	
	plink --bfile $LD_FILE \
        --chr ${CHR} \
        --from-bp ${LOWER} \
        --to-bp ${UPPER} \
        --list-duplicate-vars suppress-first \
        --make-bed \
	--out "${DIR}/${RANGE}.with.dup"


	plink --bfile "${DIR}/${RANGE}.with.dup" \
	--exclude "${DIR}/${RANGE}.with.dup.dupvar" \
	--r2 square \
	--make-bed \
	--out "${DIR}/${RANGE}"

	echo "We made it!"

done
