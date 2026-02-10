#-------------------------------------------------------
# Perform TensorQTL for CD4_AGING_COHORT
# Author: Angli XUE
# Adapted by: Marliette Matos
# Description: Mapping cis-eQTls for all cells with
#               Genotypes filtered for MAF>5 and Adding nReads to the covariate list
# adapted from https://github.com/powellgenomicslab/sc-veQTL/tree/main/TensorQTL_cis
#-------------------------------------------------------

import pandas as pd
import numpy as np
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print(f"PyTorch {torch.__version__}")
print(f"Pandas {pd.__version__}")
import os
import sys
from sys import argv

inputArgs = sys.argv[1:]
#Convert the argument to integer
chr_num = int(inputArgs[0])
celltype = "allcells"

os.makedirs(f"./analysis/cis_eQTLs/all_CD4T_cells_MAF5", exist_ok = True)
os.chdir(f"./analysis/cis_eQTLs/all_CD4T_cells_MAF5")

# define paths to data
plink_prefix_path = f"/gpfs/commons/home/mmatos/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/analysis/002.v4_calculating_pcs_MAF5/per_chr_plink_files_scRNAseq/chr{chr_num}_ashkenazi.363.AF1.QC.BA.king2.hwe.scrnaseq.annot"
expression_bed = f"/gpfs/commons/home/mmatos/cd4_QTL_analysis/02_Gene_expression/analysis/004_ProcessExpression/{celltype}_pseudo_cells_mean_mx.bed"
covariates_file = f"/gpfs/commons/home/mmatos/cd4_QTL_analysis/02_Gene_expression/analysis/005_PreparecovsFile_ncells/{celltype}_covs_PC_nFeatures_ncells_encoded_batch.txt"


prefix = f"{celltype}"

# load phenotypes and covariates
print("Load phenotypes and covariates")
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep="\t", index_col=0).T
#covariates_df.columns = 'X'+covariates_df.columns


# PLINK reader for genotypes
print("load Genotype")
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index("snp")[["chrom", "pos"]]
# Chromosome name changed format
variant_df.chrom = "chr" + variant_df.chrom


##------------HOUSEKEEPING FORMATTING CHECKS------------#
## The lengths must match to compare 
print("Match all the input files")
print(phenotype_df.columns)
print(covariates_df.columns)

# Convert column names of covariates_df to strings
covariates_df.columns = covariates_df.columns.astype(str)

## The lengths must match to compare 
print("Match all the input files")
print(phenotype_df.columns)
print(covariates_df.columns)

covariates_df = covariates_df.loc[: ,phenotype_df.columns]
covariates_df = covariates_df.T
## female = 1; male = 2
covariates_df['Sex'] = pd.factorize(covariates_df['Sex'])[0] + 1

## A2 = 1; A3 = 2; etc
#covariates_df['scRNAseq_batch'] = pd.factorize(covariates_df['scRNAseq_batch'])[0] + 1 #only use if the batches are in one colunm
covariates_df = covariates_df.apply(pd.to_numeric)

#colnames = list(map(lambda st: str.replace(st, "-", "_"), genotype_df.columns))
#genotype_df.columns = colnames
#genotype_df = genotype_df.loc[: ,phenotype_df.columns]

###-------------END OF HOUSEKEEPPING-------------------#

# genes on chr'{chr_num}'
print("genes on chr{}".format(chr_num))

## cis-QTL: nominal p-values for all variant-phenotype pairs
# map all cis-associations (results for each chromosome are written to file)

# all genes
# usage: cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, prefix)
cis.map_nominal(genotype_df, variant_df,
                phenotype_df.loc[phenotype_pos_df["chr"] == "chr" + str(chr_num)],
                phenotype_pos_df.loc[phenotype_pos_df["chr"] == "chr" + str(chr_num)],
                prefix, covariates_df = covariates_df,
                write_top = True, write_stats = True,
		maf_threshold = 0.05)

# load results
pairs_df = pd.read_parquet(f"{prefix}.cis_qtl_pairs.chr{chr_num}.parquet")
pairs_df.head()

print("Saving the all association pairs")
pairs_df.to_csv(f"{prefix}.cis_qtl_pairs.chr{chr_num}.csv", sep = "\t")

##------------------------------------------------------------------------------------------

## cis-QTL: empirical p-values for phenotypes
# all genes
# usage: cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df)
print("Start permutation test")
# genes on a specific chromosome
# This only output the top signal per gene
cis_df = cis.map_cis(genotype_df, variant_df, 
                     phenotype_df.loc[phenotype_pos_df["chr"] == "chr" + str(chr_num)],
                     phenotype_pos_df.loc[phenotype_pos_df["chr"] == "chr" + str(chr_num)],
                     covariates_df = covariates_df, seed = 123456,
                     maf_threshold = 0.05, nperm = 10000)

# calculate chromosome-wide FDR
# need to understand why the author sets lamdba = 0.85
tensorqtl.calculate_qvalues(cis_df, qvalue_lambda = 0.85)

cis_df.head()

# save the results
print("Saving the significant association pairs")
cis_df.to_csv(f"{prefix}.sig_cis_qtl_pairs.chr{chr_num}.csv", sep = "\t")

##-----------------------------------------------------------------------------------

# conditionally independent QTLs
if any(cis_df.qval < 0.05):
	print("Identify conditionally independent QTLs")
	indep_df = cis.map_independent(genotype_df, variant_df, cis_df,
                               phenotype_df, phenotype_pos_df, 
                               covariates_df = covariates_df,
                               maf_threshold = 0.05, nperm = 10000,
                               fdr=0.05, fdr_col="qval")

	# index = indep_df.index
	# indep_df = pd.merge(indep_df,allele_freq_df,on="variant_id")
	# indep_df.index = index

	print("Saving the conditionally independent QTLs")
	indep_df.to_csv(f"{prefix}.independent_cis_qtl_pairs.chr{chr_num}.csv", sep = "\t")

else :
	print("No significant eQTLs with qval < 0.05. Do not identify conditionally independent QTLs.")

print("Analaysis finished!")


####
