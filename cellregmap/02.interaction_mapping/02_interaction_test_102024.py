#-------------------------------------------------------
#### Perform CellRegmap for CD4_AGING_COHORT 
#Author: Authors of CellRegMap
#Adapted by: Marliette Matos
#Description: Script for expanding genotype and kinship matrix across pseudocells
#-------------------------------------------------------

#------------Import Modules --------------------------
import os
import re
import sys
import time
import datetime
import pandas as pd
import numpy as np
import xarray as xr
from numpy import ones
from numpy.linalg import cholesky
from pandas_plink import read_plink1_bin
import anndata
from limix.qc import quantile_gaussianize
from sklearn.preprocessing import StandardScaler
from cellregmap import run_interaction

start_time=time.time()
current_date = datetime.datetime.now().date()
run_date = current_date.strftime('%Y_%m')

#------------------------------------------------------
#--------------Set input arguments --------------------
index = int(sys.argv[1]) 
res = 34
fvf_filename = f"~/cd4_CellRegMap/002_interaction_analysis/results/results_042426/cis_eQTLS_allcells_independent_lead_snps_per_gene_res{res}.csv"

cd4_eqtl = pd.read_csv(fvf_filename, index_col = False)
egene=cd4_eqtl.loc[cd4_eqtl.index == index, "phenotype_id"].values[0]
var=cd4_eqtl.loc[cd4_eqtl.index == index, 'variant_id'].values[0]

#-------------Input and Working direcotries----------
#set outdirs
outdir = f"~/cd4_CellRegMap/002_interaction_analysis/results/results_042426/many_lead_snps_per_gene_res{res}_5factors/p_values"

try:
    os.makedirs(outdir, exist_ok=True)
    print(f"Directory created successfully: {outdir}")
except Exception as e:
    print(f"Error creating directory: {e}")


#Convert variant name to a format w/o special characters that can be used as a filename

def convert_str_pattern(input_str):
    # Use regular expressions to match and replace the pattern
    pattern = r"(\d+):(\d+)\[([a-zA-Z0-9]+)\](.*)"
    replacement = r"\1_\2_\3_\4"
    output_str = re.sub(pattern, replacement, input_str)
    # Replace commas with underscores
    output_str = output_str.replace(",", "_")
    return output_str
var_str = convert_str_pattern(var)

outfilename = os.path.join(outdir, f"{egene}_{var_str}_results.csv" )


#inputs
input_files_dir = "~/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/analysis/002.v4_calculating_pcs_MAF5/"

results="~/cd4_CellRegMap/001_preprocessing/results_01312025/prep_phenotype_vector/Data/pseudobulk_metacells/"

meta="~/metadata_harmonizing/results/post_WGS_QC/version_080923/revision_103023/cd4_scRNAseq_meta.v3.csv"

sample_mapping_file = f"~/cd4_CellRegMap/001_preprocessing/results_01312025/prep_phenotype_vector/Data/pseudobulk_metacells/pseusobulk_sample_mapping_res{res}.csv"

contexts_file = f"~/cd4_CellRegMap/001_preprocessing/results_01312025/mofa_factors/mofa_trained/cd4_aging_filt_sce_mofa_expectations_model_factors_res{res}.csv"

covariates_file=f"~/cd4_CellRegMap/001_preprocessing/results_01312025/prep_phenotype_vector/Data/pseudobulk_metacells/pseusobulk_covariates_res{res}.csv"
                 
               
print(f"Starting CellRegMap analysis for egene {egene} and variant {var}")             

############################################
########## Sample mapping file #############
############################################
print(f"Reading sample mapping file for cells in res{res}")

## this file will map cells to donors 
## it will also only include donors we have single-cell data for (a subset of all of HipSci donors)
sample_mapping = pd.read_csv(sample_mapping_file, dtype={"WGS_sampleID": str, "scRNA_Sample_Name": str})


## genotype_individual_id are donor IDs, as found in the genotype matrix (G) and GRM covariance (K)
## phenotype_sample_id are cell IDs, as found in the scRNA-seq phenotype vector (y) and cell context covariance (C)
sample_mapping.head()

## extract unique individuals
donors = sample_mapping["WGS_sampleID"].unique()
donors.sort()
print("Number of unique donors: {}".format(len(donors)))



######################################
############### Genotypes ############
######################################
print(f"Preparing genotype matrix for cells in res{res}")

#--------------Get info-------------------------------
chr_num=cd4_eqtl.loc[cd4_eqtl['phenotype_id'] == egene, 'chrom'].values[0]


#open the genotype file for the specified chromosome
plink_file = input_files_dir+f"per_chr_plink_files/chr{chr_num}_ashkenazi.367.AF1.QC.BA.king2.hwe.annot.bed"
## read in genotype file (plink format)
G = read_plink1_bin(plink_file)

print(G)

# adjust naming of samples to have a "g" in front
if isinstance(G, xr.Dataset) or isinstance(G, xr.DataArray):
    # Add "g" prefix to each sample name
    G['sample'] = ['g' + sample for sample in G['sample'].values]
print(G)

#############################
###### SNP selection

########################################
# option: testing only specific eQTLs
# filter file (columns: snp_id, gene)

## (3) select eQTLs for that gene only (from filter file)
var=cd4_eqtl.loc[cd4_eqtl.index == index, 'variant_id'].values[0]
print(var)

## (4) get genotypes
G_sel = G[:,G['snp'].isin(var)]

G_sel.shape
G_sel


#### filtering samples with missing genotypes for the eQTL variant ##########
# Step 1: Identify samples with NaN genotypes for the selected variants
nan_mask = G_sel.isnull().any(dim='variant')

# Step 2: Filter out samples with NaN genotypes
G_sel_filtered = G_sel.sel(sample=~nan_mask)

# Step 3: Filter the sample_mapping matrix for the just the samples that have a genotype 
complete_samples = G_sel_filtered['sample'].values.tolist()
filtered_sample_mapping = sample_mapping[sample_mapping['WGS_sampleID'].isin(complete_samples)]
## extract unique individuals
donors = filtered_sample_mapping["WGS_sampleID"].unique()
donors.sort()


############################################
##### expand from donors to cells ##########

# expand out genotypes from cells to donors (and select relevant donors in the same step)
#filtered_sample_mapping['WGS_sampleID'] = filtered_sample_mapping['WGS_sampleID'].astype(str)
G_expanded = G_sel.sel(sample=filtered_sample_mapping["WGS_sampleID"].values)


print("Genotype Dimensions: ")
print(G_expanded.shape)
print("Genotype Samples: " )
print(G_expanded.sample.values)

G_expanded


# unpack G
GG = G_expanded.values


############################################
################ Kinship matrix ############
############################################
print(f"Preparing Kinship matrix for cells in res{res}")

# Make king IBD matrix: plink2 --bfile {raw_genotype_filtered} --extract {out_pruning_info}.prune.in --make-king square --out {king_ibd_out}
# 
# After running this command, the output *.king and *.king.id can be made into the kinship matrix for QTL. Please make sure you multiply the king values by 2, to get in the normal 0-1 space. The kinship needs to have the king.id as row and column information. https://github.com/single-cell-genetics/limix_qtl/wiki/Inputs#genotype-file

# Read the donor names
names_file = input_files_dir + "CD4_all_chr_ashkenazi.364.AF1.QC.BA.king2.hwe.kinship.king.id"
names = pd.read_csv(names_file, sep="\t", names=["FID", "IID"],  header=0)
names['IID'] = names['IID'].astype(str)
names['IID'] = 'g' + names['IID']

donors_geno = names["IID"]

# Remove the name attribute from the donors Series
donors_geno.name = None

# Print to verify the donors Series
print(donors_geno)


# Read the GRM (genotype relationship matrix) using donor names as both rows and columns
grm_file = input_files_dir + "CD4_all_chr_ashkenazi.364.AF1.QC.BA.king2.hwe.kinship.king"
K = pd.read_csv(grm_file, sep="\t", header=None)
K.head(3) 

# Set the index and column names to the donor names
K.index = donors_geno
K.columns = donors_geno

K.head(3) 
print(K.shape)


##filter for donors with complete genotypes
filtered_donors_geno = donors_geno[donors_geno.isin(complete_samples)]
K_filtered = K.loc[filtered_donors_geno, filtered_donors_geno]
print(K_filtered.shape)


#Negative values between individuals were set to 0, as this indicates that they are less related than random individuals
# Set negative kinship coefficients to zero
K_filtered[K_filtered < 0] = 0
# Scale the kinship coefficients to range between 0 and 1
K_filtered = 2 * K_filtered
K_filtered.head(3)


# Convert to Array
K_filtered = xr.DataArray(K_filtered.values, dims=["sample_0", "sample_1"], coords={"sample_0": K_filtered.columns, "sample_1": K_filtered.index})
K_filtered = K_filtered.sortby("sample_0").sortby("sample_1")

filtered_donors_geno = sorted(set(list(K_filtered.sample_0.values)).intersection(donors))
print("Number of donors after kinship intersection: {}".format(len(filtered_donors_geno)))

## subset to relevant donors
#K = K.sel(sample_0=donors, sample_1=donors)
#assert all(K.sample_0 == donors)
#assert all(K.sample_1 == donors)

## and decompose such as K = hK @ hK.T (using Cholesky decomposition)
hK = cholesky(K_filtered.values)
hK = xr.DataArray(hK, dims=["sample", "col"], coords={"sample": K_filtered.sample_0.values})
assert all(hK.sample.values == K_filtered.sample_0.values)
print(hK)

del K_filtered
print("Sample mapping number of rows BEFORE intersection: {}".format(filtered_sample_mapping.shape[0]))

## subsample sample mapping file to donors in the kinship matrix
#sample_mapping = sample_mapping[sample_mapping["WGS_sampleID"].isin(donors)]
#print("Sample mapping number of rows AFTER intersection: {}".format(sample_mapping.shape[0]))


#-------------- expand kinship matrix from donors to cells ------------------

## use sel from xarray to expand hK (using the sample mapping file)
filtered_sample_mapping['WGS_sampleID'] = filtered_sample_mapping['WGS_sampleID'].astype(str)
hK_expanded = hK.sel(sample=filtered_sample_mapping["WGS_sampleID"].values)
assert all(hK_expanded.sample.values == filtered_sample_mapping["WGS_sampleID"].values)

print("Checking the structure of the expanded kinship matrix for res ", res, " ....")
hK_expanded
print(hK_expanded)
print(hK_expanded.sample.values)

hK = hK_expanded.values


assert all(hK_expanded.sample.values == G_expanded.sample.values)

#####################################
############ Phenotypes #############
#####################################

print(f"Preparing phenotype file for gene {egene} in res {res}")

#------------- Read Metacells Pseudobulk for each resoltuion ------------------
phenodata = results+'Pseudobulk_per_donor_all_cells_all_conditions_Leiden_res'+str(res)+".h5"
# Load AnnData object from HDF5 file
adata = anndata.read_h5ad(phenodata)
print(adata)

#filter samples based on complete genotypes
adata.obs["WGS_sampleID"] = "g" + adata.obs["WGS_sampleID"].astype(str)
filtered_adata = adata[adata.obs['WGS_sampleID'].isin(donors), :]
print(filtered_adata)

phenotype = filtered_adata.X
phenotype = pd.DataFrame(phenotype.toarray(), index=filtered_adata.obs.index, columns=filtered_adata.var.index)
phenotype = phenotype.transpose()
phenotype.columns = phenotype.columns.astype(int) #Because the pesudo cell_ids are numbers, sometimes the software treats them as int or str ...

print(phenotype)


print("Phenotype shape BEFORE selection: {}".format(phenotype.shape))
phenotype = xr.DataArray(phenotype.values, dims=["trait", "cell"], coords={"trait": phenotype.index.values, "cell": phenotype.columns.values})

phenotype = phenotype.sel(cell=filtered_sample_mapping["Cell_ID"].values)

print("Phenotype shape AFTER selection: {}".format(phenotype.shape))
assert all(phenotype.cell.values == filtered_sample_mapping["Cell_ID"].values)

y = phenotype.sel(trait=egene)
# quantile normalise
y = quantile_gaussianize(y)
# reshape
y = y.values.reshape(y.shape[0],1)
print(y.shape)


## extract unique individuals
cells = filtered_sample_mapping["Cell_ID"].unique()


######################################
############ Covariates ##############
######################################
print(f"Preparing covariates for cells in res {res}")

cov = pd.read_csv(covariates_file, index_col = 0)
cov['Sex'] = pd.factorize(cov['Sex'])[0]

# Standardize age (z-score normalization)
scaler = StandardScaler()
cov['age'] = scaler.fit_transform(cov[['age']])

selected_columns = cov[['Sex', 'age']]
selected_columns.columns = range(selected_columns.shape[1])

## subsample covariate matrix file to cells in the pehnotype matrix
selected_columns = selected_columns.loc[selected_columns.index.isin(cells)]

W = xr.DataArray(selected_columns.values, dims=["cell", "cov"], coords={"cell": selected_columns.index.values, "cov": selected_columns.columns.values})
W = W.sel(cell=filtered_sample_mapping["Cell_ID"].values)
W = W.values
#W= np.array(selected_columns)
#W=np.asmatrix(W)

# Create an intercept column of ones
#W_intercept = np.ones((n_cells, 1))
# Combine the intercept column with the covariates DataFrame (optional)
#W_matrix = np.hstack((W_intercept, selected_columns.values))


#batch_columns = cov.filter(regex='^Batch_')

# Concatenate the selected columns with the batch columns
#cov = pd.concat([selected_columns, batch_columns], axis=1)
#######COVARIATE MATRIX 
# just an intercept in this case
#n_cells = phenotype.shape[1]
#W = ones((n_cells, 1))

######################################
########## Cell contexts #############
######################################
print(f"Preparing cellular context for cells in res {res}")

# cells by MOFA factors (5) but only 4 had variants of at least 1%
C = pd.read_csv(contexts_file, index_col = 0)
C = C.iloc[:, :5]

#I decided to just test the the first 5 mofa factors only 
#celltype = cov.filter(regex='^cs_')
#age = cov.filter(regex='age')

## subsample contexts matrix file to cells in the pehnotype matrix
C = C.loc[C.index.isin(cells)]

#add other states to the cell contexts
#C = C.merge(celltype, left_index=True, right_index=True).merge(age, left_index=True, right_index=True)
C = xr.DataArray(C.values, dims=["cell", "pc"], coords={"cell": C.index.values, "pc": C.columns.values})
C = C.sel(cell=filtered_sample_mapping["Cell_ID"].values)

assert all(C.cell.values == filtered_sample_mapping["Cell_ID"].values)


# quantile normalise cell contexts
C_gauss = quantile_gaussianize(C)
print(C_gauss.head)
C_gauss = C_gauss.values



# hK_expanded and G_expanded should be inputed in CellRegMap
# e.g.,:

# from cellregmap import run_interaction
# pv = run_interaction(y=y, G=G_exp, W=W, E=E, hK=hK_exp)[0] 
#note adding the [0] only extracts the pvalues for interactions for the first context 

# run association test using CellRegMap
print(f"Running Interaction Analysis for gene {egene}")

pvals = run_interaction(y=y, W=W, E=C_gauss, G=GG, hK=hK)[0]

pv = pd.DataFrame({"chrom":G_expanded.chrom.values,
               "pv":pvals,
               "variant":G_expanded.snp.values})

print(f"Showing CellRegMap Interaction results for gene {egene}")
pv.head()

print(f"Writing CellRegMap Interaction results for gene {egene} to {outfilename}")
pv.to_csv(outfilename)

print(f"Interaction Analysis completed for gene {egene}")
end_time=time.time()
duration = end_time - start_time
print(f"Duration: {duration} seconds")

print('Process completed!')
print('Good-bye!')
