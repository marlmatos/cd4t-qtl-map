import pandas as pd
import pybedtools
import argparse
pd.set_option('display.max_columns', 20)


def fetch_variant_annotation_args():
    """
    Parse command line arguments for variant annotation
    """
    parser = argparse.ArgumentParser(description='Annotate variant scores')
    parser.add_argument('--list', type=str, required=True, 
                        help='Path to the variant scores file')
    parser.add_argument('--out_prefix', type=str, required=True, 
                        help='Output prefix')
    parser.add_argument('--peaks', type=str, default=None, 
                        help='Path to peaks file for overlap annotation')
    parser.add_argument('--genes', type=str, default=None, 
                        help='Path to genes file for closest gene annotation')
    parser.add_argument('--schema', type=str, default='chrombpnet', 
                        help='Schema for variant scores file')
    return parser.parse_args()


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


def map_chrombpnet_columns(df):
    """
    Map chrombpnet schema columns to standard column names
    """
    df_copy = df.copy()
    
    # Map column names
    col_mapping = {
        'CHR': 'chr',
        'POS0': 'pos',
        'REF': 'allele1',
        'ALT': 'allele2'
    }
    
    for src, dest in col_mapping.items():
        if src in df_copy.columns:
            df_copy[dest] = df_copy[src]
    
    # Create variant_id if it doesn't exist
    if 'variant_id' not in df_copy.columns:
        df_copy['variant_id'] = df_copy.apply(
            lambda row: f"{row['chr']}_{row['pos']}_{row['allele1']}_{row['allele2']}", 
            axis=1
        )
    
    return df_copy


def main():
    args = fetch_variant_annotation_args()
    print(args)
    variant_scores_file = args.list
    output_prefix = args.out_prefix
    peak_path = args.peaks
    genes_path = args.genes
    schema = args.schema

    # Read variant scores file
    variant_scores = pd.read_table(variant_scores_file)
    print(f"Read variant scores file with {variant_scores.shape[0]} rows and {variant_scores.shape[1]} columns")
    
    # Map columns based on schema
    if schema == "chrombpnet":
        variant_scores = map_chrombpnet_columns(variant_scores)
        print("Mapped chrombpnet columns to standard column names")
    
    # Convert to bed format
    if schema == "bed":
        if 'end' in variant_scores.columns and variant_scores['pos'].equals(variant_scores['end']):
            variant_scores['pos'] = variant_scores['pos'] - 1
        variant_scores_bed_format = variant_scores[['chr','pos','end','allele1','allele2','variant_id']].copy()
        variant_scores_bed_format.sort_values(by=["chr","pos","end"], inplace=True)
    else:
        # Convert to bed format
        variant_scores_bed_format = variant_scores[['chr','pos','allele1','allele2','variant_id']].copy()
        variant_scores_bed_format['pos'] = variant_scores_bed_format.apply(lambda x: int(x.pos)-1, axis=1)
        variant_scores_bed_format['end'] = variant_scores_bed_format.apply(lambda x: int(x.pos)+len(str(x.allele1)), axis=1)
        variant_scores_bed_format = variant_scores_bed_format[['chr','pos','end','allele1','allele2','variant_id']]
        variant_scores_bed_format.sort_values(by=["chr","pos","end"], inplace=True)

    print("\nVariant BED format preview:")
    print(variant_scores_bed_format.head())
    print(f"Variants table shape: {variant_scores_bed_format.shape}")

    # Create BedTool object
    variant_bed = pybedtools.BedTool.from_dataframe(variant_scores_bed_format)

    # Annotate with closest genes if specified
    if args.genes:
        print("\nAnnotating with closest genes...")
        try:
            # Read genes BED file - our custom format has gene name in column 4
            gene_df = pd.read_csv(genes_path, sep='\t', header=None)
            print(f"Read genes file with {gene_df.shape[0]} rows and {gene_df.shape[1]} columns")
            
            # Determine gene info columns based on file format
            gene_name_col = 3  # BED format has name in column 4 (0-indexed: 3)
            
            # Column renaming for pybedtools
            gene_columns = list(range(gene_df.shape[1]))
            gene_df.columns = gene_columns
            
            print("Gene columns:")
            for i in range(min(10, gene_df.shape[1])):
                print(f"Column {i}: {gene_df[i].iloc[0] if i < gene_df.shape[1] else 'N/A'}")
            
            # Convert to BedTool
            gene_bed = pybedtools.BedTool.from_dataframe(gene_df)
            
            # Find closest gene(s) for each variant
            closest_genes_bed = variant_bed.closest(gene_bed, d=True, t='first', k=3)
            
            # Convert results to DataFrame
            closest_gene_df = closest_genes_bed.to_dataframe(header=None)
            
            print("Closest genes preview:")
            print(closest_gene_df.head())
            print(f"Closest genes table shape: {closest_gene_df.shape}")
            
            # Column indices - determine dynamically based on shape
            variant_id_idx = 5  # Based on earlier definition
            gene_name_idx = variant_id_idx + 1 + gene_name_col  # Offset by variant columns + 1
            distance_idx = closest_gene_df.shape[1] - 1  # Last column is distance
            
            print(f"Using variant_id at column {variant_id_idx}, gene name at column {gene_name_idx}")
            
            # Process closest genes information
            closest_genes = {}
            gene_dists = {}
            gene_types = {}
            gene_ids = {}
            
            # Additional annotation columns to extract (adjust indices based on your BED format)
            gene_id_idx = gene_name_idx + 3  # Assuming gene_id is 3 columns after gene_name
            gene_type_idx = gene_name_idx + 4  # Assuming gene_type is 4 columns after gene_name
            
            # Check if these columns exist
            has_gene_id = gene_id_idx < closest_gene_df.shape[1]
            has_gene_type = gene_type_idx < closest_gene_df.shape[1]
            
            for index, row in closest_gene_df.iterrows():
                variant_id = row[variant_id_idx]
                
                if variant_id not in closest_genes:
                    closest_genes[variant_id] = []
                    gene_dists[variant_id] = []
                    gene_ids[variant_id] = []
                    gene_types[variant_id] = []
                
                # Add gene name and distance
                gene_name = row[gene_name_idx]
                distance = row[distance_idx]
                
                closest_genes[variant_id].append(gene_name)
                gene_dists[variant_id].append(distance)
                
                # Add gene ID and type if available
                if has_gene_id:
                    gene_ids[variant_id].append(row[gene_id_idx])
                else:
                    gene_ids[variant_id].append('.')
                    
                if has_gene_type:
                    gene_types[variant_id].append(row[gene_type_idx])
                else:
                    gene_types[variant_id].append('.')
            
            # Create dataframe with closest gene information
            closest_gene_df = pd.DataFrame({'variant_id': list(closest_genes.keys())})
            
            # Add closest genes and distances (up to 3 closest genes)
            for i in range(1, 4):
                closest_gene_df[f'closest_gene_{i}'] = closest_gene_df['variant_id'].apply(
                    lambda x: closest_genes[x][i-1] if i <= len(closest_genes[x]) else '.'
                )
                closest_gene_df[f'gene_distance_{i}'] = closest_gene_df['variant_id'].apply(
                    lambda x: gene_dists[x][i-1] if i <= len(gene_dists[x]) else '.'
                )
                
                # Add gene ID and type if available
                if has_gene_id:
                    closest_gene_df[f'gene_id_{i}'] = closest_gene_df['variant_id'].apply(
                        lambda x: gene_ids[x][i-1] if i <= len(gene_ids[x]) else '.'
                    )
                
                if has_gene_type:
                    closest_gene_df[f'gene_type_{i}'] = closest_gene_df['variant_id'].apply(
                        lambda x: gene_types[x][i-1] if i <= len(gene_types[x]) else '.'
                    )
            
            # Remove duplicates and merge with variant scores
            closest_gene_df.drop_duplicates(inplace=True)
            variant_scores = variant_scores.merge(closest_gene_df, on='variant_id', how='left')
            
            print("Added closest gene annotations")
            
        except Exception as e:
            print(f"Error during gene annotation: {e}")
            import traceback
            traceback.print_exc()

    # Annotate with peak overlap if specified
    if args.peaks:
        print("\nAnnotating with peak overlap...")
        try:
            peak_df = pd.read_table(peak_path, header=None)
            peak_bed = pybedtools.BedTool.from_dataframe(peak_df)
            peak_intersect_bed = variant_bed.intersect(peak_bed, wa=True, u=True)

            peak_intersect_df = peak_intersect_bed.to_dataframe(names=variant_scores_bed_format.columns.tolist())
            
            print("Peak overlap preview:")
            print(peak_intersect_df.head())
            print(f"Peak overlap table shape: {peak_intersect_df.shape}")

            # Mark variants that overlap with peaks
            variant_scores['peak_overlap'] = variant_scores['variant_id'].isin(peak_intersect_df['variant_id'].tolist())
            
            print("Added peak overlap annotations")
        except Exception as e:
            print(f"Error during peak annotation: {e}")

    # Preview final results
    print("\nFinal annotations preview:")
    print(variant_scores.head())
    print(f"Final annotation table shape: {variant_scores.shape}")

    # Save results
    out_file = f"{output_prefix}.annotations.tsv"
    variant_scores.to_csv(out_file, sep="\t", index=False)
    print(f"\nResults saved to {out_file}")
    print("DONE")


if __name__ == "__main__":
    main()