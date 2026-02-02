#!/bin/bash
#SBATCH --job-name=scrnaseq_pip         # create a short name for your job
#SBATCH -p pe2 
#SBATCH -t 7-00:00:00
#SBATCH --cpus-per-task=4                   
#SBATCH --mem=16G                     # Job memory reques
#SBATCH --mail-type=END,FAIL             # send email when job begins
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH --output=002_nextflow-%j.out

source ~/.bashrc
module load singularity/3.8.6

conda activate nextflow.personal
nextflow run scRNAseq_preprocessing2.groovy -c nextflow.config -resume 
