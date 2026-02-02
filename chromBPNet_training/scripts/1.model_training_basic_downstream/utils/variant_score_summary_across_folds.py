import pandas as pd
import numpy as np
import os
import argparse

## This script is a modified version of the variant_summary_across_folds.py script from the chrombpnet/variant-scorer repository.
## The original script reads variant scores from multiple files and calculates the mean and geometric mean of the scores.
## The script has been modified to allow for different variant file formats and to handle non-numeric score columns.

def get_variant_schema(schema):
    """
    Returns column schema for different variant file formats
    """
    var_SCHEMA = {
        'original': ['chr', 'pos', 'variant_id', 'allele1', 'allele2'],
        'plink': ['chr', 'variant_id', 'ignore1', 'pos', 'allele1', 'allele2'],
        'plink2': ['chr', 'variant_id', 'pos', 'allele1', 'allele2'],
        'bed': ['chr', 'pos', 'end', 'allele1', 'allele2', 'variant_id'],
        'chrombpnet': ['CHR', 'POS0', 'REF', 'ALT', 'META_DATA']
    }
    
    if schema not in var_SCHEMA:
        raise ValueError(f"Unknown schema '{schema}'. Available schemas: {list(var_SCHEMA.keys())}")
    
    return var_SCHEMA[schema]


def geo_mean_overflow(iterable, axis=0):
    """
    Compute geometric mean avoiding overflow
    """
    a = np.array(iterable)
    
    # Handle zero or negative values by replacing with NaN
    a = np.where(a <= 0, np.nan, a)
    
    # Calculate log, ignoring NaNs
    log_a = np.log(a)
    
    # Calculate mean of logs, then exp
    return np.exp(np.nanmean(log_a, axis=axis))


def parse_arguments():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(description='Process variant scores')
    parser.add_argument('--score_dir', type=str, required=True, 
                        help='Directory containing score files')
    parser.add_argument('--score_list', nargs='+', required=True, 
                        help='List of score files')
    parser.add_argument('--out_prefix', type=str, required=True, 
                        help='Output prefix')
    parser.add_argument('--schema', type=str, default='chrombpnet', 
                        help='Variant schema (default: chrombpnet)')
    return parser.parse_args()


def map_columns(df, schema='chrombpnet'):
    """
    Map columns from a specific schema to standard column names
    """
    column_mapping = {
        'chrombpnet': {
            'CHR': 'chr',
            'POS0': 'pos',
            'REF': 'allele1',
            'ALT': 'allele2',
            'META_DATA': 'meta_data'
        }
    }
    
    # Create a copy to avoid modifying the original
    result = df.copy()
    
    # Map columns based on schema
    if schema in column_mapping:
        for src, dest in column_mapping[schema].items():
            if src in result.columns and dest not in result.columns:
                result[dest] = result[src]
    
    # Create variant_id if it doesn't exist
    if 'variant_id' not in result.columns:
        # Check which columns to use for creating variant_id
        chr_col = 'chr' if 'chr' in result.columns else 'CHR'
        pos_col = 'pos' if 'pos' in result.columns else 'POS0'
        ref_col = 'allele1' if 'allele1' in result.columns else 'REF'
        alt_col = 'allele2' if 'allele2' in result.columns else 'ALT'
        
        # Create variant_id
        if all(col in result.columns for col in [chr_col, pos_col, ref_col, alt_col]):
            result['variant_id'] = result.apply(
                lambda row: f"{row[chr_col]}_{row[pos_col]}_{row[ref_col]}_{row[alt_col]}", 
                axis=1
            )
    
    return result


def main():
    # Parse command line arguments
    args = parse_arguments()
    print(args)
    
    variant_score_dir = args.score_dir
    variant_table_list = args.score_list
    output_prefix = args.out_prefix
    schema = args.schema
    
    # Get schema columns
    schema_columns = get_variant_schema(schema)
    print(f"Using schema: {schema} with columns: {schema_columns}")
    
    # Read and process score files
    score_dict = {}
    for i, filename in enumerate(variant_table_list):
        variant_score_file = os.path.join(variant_score_dir, filename)
        
        if not os.path.isfile(variant_score_file):
            print(f"Warning: File {variant_score_file} does not exist. Skipping.")
            continue
            
        try:
            var_score = pd.read_table(variant_score_file)
            print(f"Read file {filename} with {var_score.shape[0]} rows and {var_score.shape[1]} columns")
            
            # Map columns to standard names
            var_score = map_columns(var_score, schema)
            
            score_dict[i] = var_score
        except Exception as e:
            print(f"Error reading file {filename}: {e}")
    
    if not score_dict:
        print("No valid score files found. Exiting.")
        return
    
    # Use the first file as reference
    first_key = list(score_dict.keys())[0]
    first_file = score_dict[first_key]
    
    # Create output dataframe
    variant_scores = pd.DataFrame()
    
    # Copy over schema columns from first score file
    standard_cols = ['chr', 'pos', 'allele1', 'allele2', 'variant_id']
    original_cols = schema_columns
    
    # Ensure standard columns exist in the output
    for col in standard_cols:
        if col in first_file.columns:
            variant_scores[col] = first_file[col]
    
    # Also add original schema columns if they're not already added
    for col in original_cols:
        if col in first_file.columns and col not in variant_scores.columns:
            variant_scores[col] = first_file[col]
    
    # Verify that all files have the same variant positions
    consistent = True
    for i in score_dict:
        if i == first_key:  # Skip the first file
            continue
            
        for col in ['chr', 'pos', 'allele1', 'allele2']:
            if col in variant_scores.columns and col in score_dict[i].columns:
                if not score_dict[i][col].equals(variant_scores[col]):
                    print(f"Warning: Column {col} differs between files")
                    consistent = False
    
    if not consistent:
        print("Warning: Files have inconsistent variants. Results may not be meaningful.")
    
    # Find all score columns present in the files
    all_score_columns = set()
    for i in score_dict:
        for col in score_dict[i].columns:
            # Assume any column not in the schema is a score column
            if col not in schema_columns and col not in ['chr', 'pos', 'allele1', 'allele2', 'variant_id', 'meta_data']:
                all_score_columns.add(col)
    
    print(f"Found {len(all_score_columns)} score columns: {', '.join(sorted(all_score_columns))}")
    
    # Process each score column
    for score in sorted(all_score_columns):
        # Skip columns that are already aggregate results
        if '.mean' in score or '.pval' in score:
            continue
            
        # Find which files contain this score
        files_with_score = [i for i in score_dict if score in score_dict[i].columns]
        
        if files_with_score:
            try:
                # Calculate mean across all files
                values = [score_dict[i][score].values for i in files_with_score]
                
                # Check if we can convert to numeric (some scores might be strings)
                numeric_values = []
                for val_array in values:
                    try:
                        numeric_values.append(pd.to_numeric(val_array))
                    except:
                        print(f"Warning: Column {score} contains non-numeric values in some files. Skipping.")
                        break
                
                if len(numeric_values) == len(values):
                    variant_scores[f"{score}.mean"] = np.mean(numeric_values, axis=0)
                    print(f"Calculated mean for {score}")
                    
                    # Check for p-value columns and process them
                    pval_col = f"{score}.pval"
                    alt_pval_col = f"{score}_pval"
                    
                    pval_files = [i for i in files_with_score if pval_col in score_dict[i].columns]
                    alt_pval_files = [i for i in files_with_score if alt_pval_col in score_dict[i].columns]
                    
                    if pval_files:
                        pval_values = [score_dict[i][pval_col].values for i in pval_files]
                        variant_scores[f"{score}.mean.pval"] = geo_mean_overflow(pval_values)
                        print(f"Calculated geometric mean for {pval_col}")
                    elif alt_pval_files:
                        alt_pval_values = [score_dict[i][alt_pval_col].values for i in alt_pval_files]
                        variant_scores[f"{score}.mean.pval"] = geo_mean_overflow(alt_pval_values)
                        print(f"Calculated geometric mean for {alt_pval_col}")
            except Exception as e:
                print(f"Error processing score {score}: {e}")
    
    print()
    print("First few rows of the output:")
    print(variant_scores.head())
    print(f"Output shape: {variant_scores.shape[0]} rows x {variant_scores.shape[1]} columns")
    print()
    
    # Save results
    out_file = f"{output_prefix}.mean.variant_scores.tsv"
    variant_scores.to_csv(out_file, sep="\t", index=False)
    print(f"Results saved to {out_file}")
    print("DONE")


if __name__ == "__main__":
    main()