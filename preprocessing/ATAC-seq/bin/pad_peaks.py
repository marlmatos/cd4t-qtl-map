#!/usr/bin/env python3

import pandas as pd
import numpy as np
from typing import Dict, Set
from intervaltree import IntervalTree
import argparse

def process_macs2_summits_hybrid(input_file: str, output_file: str, chrom_lengths_file: str,
                                 extension: int = 300, score_threshold: float = 0.05, 
                                 verbose: bool = True, verbose_freq: int = 1e5) -> pd.DataFrame:
    """
    Process MACS2 narrowpeaks file with summits, padding peaks, filtering based on score,
    removing overlaps, and ensuring peaks are within chromosome boundaries.
    """
    
    # Read the chromosome lengths file
    chrom_lengths = pd.read_csv(chrom_lengths_file, sep='\t', header=None, names=['chrom', 'start', 'end'])
    chrom_lengths = dict(zip(chrom_lengths['chrom'], chrom_lengths['end']))
    
    # Read the narrowpeaks file
    df = pd.read_csv(input_file, sep='\t', header=None,
                     names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'summit'])
    
    # Calculate score threshold based on FDR
    df.score = df.score / (df.score.sum() / 1e6)
    #score_threshold = -10 * np.log10(fdr_threshold)
    #score_threshold = 0.05
    
    # Filter based on score
    df = df[df['score'] >= score_threshold]
    
    # Calculate summit position and new start/end
    df['summit_pos'] = df['start'] + df['summit']
    df['new_start'] = df['summit_pos'] - extension
    df['new_end'] = df['summit_pos'] + extension
    
    # Ensure new_start is not negative and new_end doesn't exceed chromosome length
    df['new_start'] = df['new_start'].clip(lower=0)
    df['new_end'] = df.apply(lambda row: min(row['new_end'], chrom_lengths.get(row['chrom'], row['new_end'])), axis=1)
    
    # Sort by score descending
    df = df.sort_values(by='score', ascending=False).reset_index(drop=True)
    
    summit_list: Dict[str, np.ndarray] = df.groupby('chrom')['summit_pos'].apply(lambda x: np.array(sorted(x))).to_dict()
    summit_peakname_list: Dict[tuple, np.ndarray] = df.groupby(['chrom', 'summit_pos'])['name'].apply(np.array).to_dict()
    
    to_drop: Set[str] = set()
    to_keep: Set[str] = set()
    
    for counter, (_, row) in enumerate(df.iterrows()):
        if verbose and counter % verbose_freq == 0:
            print(f'{counter} records have been checked...')
        
        chrom, start, end, peak_name, score = row[['chrom', 'new_start', 'new_end', 'name', 'score']]
        summit_pos = row['summit_pos']
        
        if peak_name in to_drop:
            continue 
                    
        range_to_check = np.arange(start=start, stop=end)
        
        intersection = np.intersect1d(range_to_check, summit_list[chrom])
        
        to_keep.add(peak_name)

        for overlap in intersection:
            for overlapping_peak_name in summit_peakname_list[(chrom, overlap)]:
                if overlapping_peak_name == peak_name:
                    continue
                to_drop.add(overlapping_peak_name)
    
    # Keep only non-overlapping peaks
    df_no_overlaps = df[df['name'].isin(to_keep)].copy()
    
    # Function to merge overlapping intervals
    def merge_intervals(intervals):
        sorted_intervals = sorted(intervals, key=lambda x: x[0])
        merged = []
        for interval in sorted_intervals:
            if not merged or merged[-1][1] < interval[0]:
                merged.append(interval)
            else:
                merged[-1] = (merged[-1][0], max(merged[-1][1], interval[1]))
        return merged
    
     # Process each chromosome separately to handle any remaining overlaps
    final_peaks = []
    for chrom, group in df_no_overlaps.groupby('chrom'):
        tree = IntervalTree()
        for idx, row in group.iterrows():
            tree[row['new_start']:row['new_end']] = idx

        # Find and merge overlapping intervals
        overlaps = tree.overlap(group['new_start'].min(), group['new_end'].max())
        non_overlapping = set(group.index) - set([item.data for item in overlaps])

        # Convert the set of non-overlapping indices to a list
        non_overlapping = list(non_overlapping)

        # Add non-overlapping peaks directly
        final_peaks.extend(group.loc[non_overlapping].to_dict('records'))

        # Merge overlapping peaks
        overlap_groups = []
        for overlap in overlaps:
            overlap_groups.append((overlap.begin, overlap.end))
        
        merged = merge_intervals(overlap_groups)
        
        for start, end in merged:
            overlapping_peaks = group[(group['new_start'] >= start) & (group['new_end'] <= end)]
            best_peak = overlapping_peaks.loc[overlapping_peaks['score'].idxmax()].copy()
            best_peak['new_start'] = start
            best_peak['new_end'] = end
            final_peaks.append(best_peak.to_dict())
    
    # Create final dataframe
    final_df = pd.DataFrame(final_peaks)
    
    # Reorder columns and rename new_start/new_end to start/end
    final_df = final_df[['chrom', 'new_start', 'new_end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'summit', 'summit_pos']]
    final_df = final_df.rename(columns={'new_start': 'start', 'new_end': 'end'})
    
    # Sort by chromosome and start position
    final_df = final_df.sort_values(['chrom', 'start'])
    
    if verbose:
        print(f"Initial peaks: {len(df)}")
        print(f"Peaks after overlap removal: {len(df_no_overlaps)}")
        print(f"Final peaks after merging: {len(final_df)}")
    
    # Save to file
    final_df.to_csv(output_file, sep='\t', index=False, header=False)
    
    return final_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process MACS2 summits for Nextflow")
    parser.add_argument("--input", required=True, help="Path to the input narrowpeaks file")
    parser.add_argument("--output", required=True, help="Path to save the processed peaks")
    parser.add_argument("--chrom_sizes", required=True, help="Path to the chromosome sizes file")
    parser.add_argument("--extension", type=int, default=300, help="Number of base pairs to extend on each side of the summit")
    parser.add_argument("--score_threshold", type=float, default=0.05, help="score threshold for filtering peaks")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    args = parser.parse_args()

    process_macs2_summits_hybrid(
        input_file=args.input,
        output_file=args.output,
        chrom_lengths_file=args.chrom_sizes,
        extension=args.extension,
        score_threshold=args.score_threshold,
        verbose=args.verbose
    )
    

# Example usage:
# processed_peaks = process_macs2_summits_hybrid('input_narrowpeaks.bed', 'output_processed_peaks.bed', 'chrom_lengths.bed')

# input="/gpfs/commons/groups/lappalainen_lab/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing/merged_library/peak_calling/MACS3/BAMPE/peaks_102024/cd4_atac_summits_peaks.narrowPeak"

# output="/gpfs/commons/groups/lappalainen_lab/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing/STAR/MACS3/anthony_concensus_peaks/marli_processed_peaks.bed"

# chrom_sizes="/gpfs/commons/home/mmatos/resources/genome/hg38.p14.chrom.sizes.fmtd"

# process_macs2_summits_hybrid(input, output, chrom_sizes ,
#                                                extension=300, fdr_threshold=0.01, verbose=True)