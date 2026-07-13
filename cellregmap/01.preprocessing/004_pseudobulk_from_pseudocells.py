#-------------------------------------------------------
# Create pseudobulk from sce metacells for CD4 T+ Cells
# Adapted by: Marliette Matos from CellRegMap analysis
# Github: https://github.com/annacuomo/CellRegMap_analyses/blob/main/neuroseq/preprocessing/pseudobulk_from_pseudocells.py
#-------------------------------------------------------


#------------------libraries------------------------------
import sys
import os
import re
import argparse
import numpy as np
import pandas as pd
from anndata import AnnData
import scanpy as sc
import scanpy.external as sce
from scipy import sparse

#-----------------Define Input Arguments------------------------
parser = argparse.ArgumentParser()
parser.add_argument('data_dir', help="absolute path of the directory containing the metacells", type=str)
args = parser.parse_args()

#------------------Prep-variables for pseudobulk ------------
def list_files(directory):
    """List all files in directory and subdirectories"""
    files_list = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.h5'):
                files_list.append(os.path.join(root, file))
    return files_list

# List all files in the directory
files = list_files(args.data_dir)

# Filter for our specific files
pattern = re.compile(r'Metacells_per_donor/(Old|Young)/allcells_(Old|Young)_donors\.h5$')
metacell_files = [f for f in files if pattern.search(f)]

# Extract age categories
age_categories = []
for filepath in metacell_files:
    match = pattern.search(filepath)
    if match:
        age_category = match.group(1)  # Will be either 'Old' or 'Young'
        age_categories.append(age_category)

age_categories = list(set(age_categories))  # Get unique categories
print(f"Found age categories: {age_categories}")
print(f"Found files: {metacell_files}")

# Dictionary for resolution values
res = {
    'leiden_res_34_euclidean': "34",
    'leiden_res_4_euclidean': "4", 
    'leiden_res_45_euclidean': "45",
    'leiden_res_5_euclidean': "5"
}

#------------------Pseudobulk per resolution----------------
for k, v in res.items():
    print(f"\n------- Resolution: {k} --------")
    init = True
    pseudobulk = None
    metadata = None
    gene_names = None
    
    # Process each age group
    for age_group in age_categories:
        print(f"\n------- Processing {age_group} group -------")
        
        # Find file for this age group
        current_file = [f for f in metacell_files if f'/Metacells_per_donor/{age_group}/' in f][0]
        print(f"Processing file: {current_file}")
        
        # Read the file
        myfile = sc.read(current_file)
        print(f"File shape: {myfile.shape}")
        
        # Store gene names if not already stored
        if gene_names is None:
            gene_names = myfile.var_names
        
        # Process each donor
        for i, d in enumerate(myfile.obs.scRNA_Sample_Name.unique(), 1):
            print(f"Donor {i} out of {myfile.obs.scRNA_Sample_Name.nunique()}")
            ct_donor = myfile[myfile.obs.scRNA_Sample_Name == d]
            cells = list(ct_donor.obs.groupby(by=k, observed=False).groups.values())
            
            # Initialize donor-specific data
            donor_pseudobulk = []
            donor_metadata = pd.DataFrame()
            
            # Process each group of cells
            for idx, c in enumerate(cells):
                cell_data = ct_donor[c]
                donor_pseudobulk.append(cell_data.X.mean(axis=0))
                
                # Get metadata for this group
                metadata_cols = ['WGS_sampleID', 'scRNA_Sample_Name', 'cd4_subtypes_l1', 
                               'Sex', 'Group', 'scRNA_batch_Dave.x', k]
                metadata_row = cell_data.obs.loc[:, metadata_cols].drop_duplicates().reset_index(drop=True).iloc[0:1]
                donor_metadata = pd.concat([donor_metadata, metadata_row], ignore_index=True)
            
            # Convert donor_pseudobulk to numpy array
            donor_pseudobulk = np.vstack(donor_pseudobulk)
            
            # Get number of cells per pseudocell
            number_of_cells_list = [len(c) for c in cells]
            
            # Verify lengths match
            assert len(donor_metadata.index) == len(number_of_cells_list), "Length mismatch before concatenation"
            
            # Add cell counts and IDs to metadata
            donor_metadata = pd.concat([
                donor_metadata,
                pd.DataFrame({
                    "Number_of_cells": number_of_cells_list,
                    "cell_ids": cells
                }, index=donor_metadata.index)
            ], axis=1)
            
            # Add to overall data
            if init:
                pseudobulk = donor_pseudobulk
                metadata = donor_metadata
                init = False
            else:
                pseudobulk = np.vstack([pseudobulk, donor_pseudobulk])
                metadata = pd.concat([metadata, donor_metadata], ignore_index=True)

    # After processing all age groups for this resolution
    if pseudobulk is not None and metadata is not None:
        print(f"\nCreating AnnData object for resolution {k}")
        print(f"Pseudobulk shape: {pseudobulk.shape}")
        print(f"Metadata shape: {metadata.shape}")
        
        pseudobulk_res = AnnData(
            X=pseudobulk,
            obs=metadata,
            var=pd.DataFrame(index=gene_names)
        )
        
        # Convert to sparse matrix
        pseudobulk_res.X = sparse.csr_matrix(pseudobulk_res.X)
        
        # Convert cell_ids to string representation
        pseudobulk_res.obs.cell_ids = [
            str(pseudobulk_res.obs.cell_ids.values[i].tolist())
            for i in range(pseudobulk_res.obs.shape[0])
        ]
        
        # Set up output directory
        res_dir = "~/cd4_CellRegMap/001_preprocessing/results_01312025/prep_phenotype_vector/Data/"
        resdir = os.path.join(res_dir, "pseudobulk_metacells")
        os.makedirs(resdir, exist_ok=True)
        
        # Save the pseudobulk data
        output_file = os.path.join(resdir, f"Pseudobulk_per_donor_all_cells_all_conditions_Leiden_res{v}.h5")
        print(f"Saving to: {output_file}")
        pseudobulk_res.write(output_file)
        print(f"Completed resolution {k}")

print("\nPseudobulk of Metacells per Resolution Done")
print("Goodbye")