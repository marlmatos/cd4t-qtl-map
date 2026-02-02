import pandas as pd
import argparse

def process_narrow_peaks(narrow_peaks, prefix, output_file):
    # Read the narrow peaks file into a DataFrame
    df = pd.read_csv(narrow_peaks, sep='\t', header=None)
    
    # Drop redundant peaks
    res = df.groupby([0, 1, 2]).first().reset_index()
    
    # Filter peaks that start with "chr"
    res = res[res[0].str.startswith("chr")]
    
    # Rename peak IDs
    res[3] = res[3].apply(lambda x: f'{prefix}_peak_{"".join(filter(str.isdigit, x[10:]))}')
    
    # Drop additional information columns
    res.iloc[:, 6:] = '.'
    
    # Save the processed DataFrame to a new narrowPeak file
    res.to_csv(output_file, sep='\t', header=None, index=None)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process narrow peaks")
    parser.add_argument("narrow_peaks", help="Input narrow peaks file")
    parser.add_argument("prefix", help="Prefix for the output file")
    parser.add_argument("output_file", help="Output narrow peaks file")

    args = parser.parse_args()
    
    process_narrow_peaks(args.narrow_peaks, args.prefix, args.output_file)

