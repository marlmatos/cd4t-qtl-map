#!/bin/bash
#SBATCH --job-name=atac_nextflow         # create a short name for your job
#SBATCH -p pe2 
#SBATCH -t 7-00:00:00
#SBATCH --cpus-per-task=4                   
#SBATCH --mem=5G                     # Job memory reques
#SBATCH --mail-type=END,FAIL             # send email when job begins
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH --output=002_nextflow-%j.out

module load miniconda3/3.22.0
module load singularity/3.8.6
module load java

source activate nextflow.personal
nextflow run main.groovy -c nextflow.config  -resume sleepy_murdock




