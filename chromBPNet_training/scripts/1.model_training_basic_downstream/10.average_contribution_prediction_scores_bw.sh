#!/bin/bash
#SBATCH --job-name=average_prediction_scores         # Job name
#SBATCH --partition=cpu                     # Partition Name
#SBATCH --mem=60G
#SBATCH --time=7-00:00:00                     # total run time limit (HH:MM:SS)
#SBATCH --output=logs/average_prediction_scores_%A_%a.log  # Standard output with array job ID
#SBATCH --error=logs/average_prediction_scores_%A_%a.err   # Error log with array job ID
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mmatos@nygenome.org     # Where to send mail

# Summary: This script averages the prediction scores and contribution scores across the 5 folds of the ChromBPNet model.

# Load module environment
#source /usr/share/lmod/lmod/init/bash

# Activate wiggletools environment
source activate wiggletools
module load deeptools
# full home path
HOME_DIR=/gchm

# directories
CONTRIB_INPUT_DIR=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/contribution_scores_bw
CONTRIB_OUTPUT_DIR=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/contribution_scores_bw/averaged

PRED_INPUT_DIR=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/prediction_scores_bw
PRED_OUTPUT_DIR=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/prediction_scores_bw/averaged
chromsizes=$HOME_DIR/cd4_chrombpnet/data/inputs/hg38.autosomes.chrom.sizes

# Create output directories if they don't exist
mkdir -p $CONTRIB_OUTPUT_DIR
mkdir -p $PRED_OUTPUT_DIR
mkdir -p $CONTRIB_OUTPUT_DIR/temp
mkdir -p $PRED_OUTPUT_DIR/temp

# The original averaging code is commented out as these files have already been created
: '
# Process contribution scores
echo "Processing contribution scores..."
cd $CONTRIB_INPUT_DIR/fold_0
CONTRIB_FILES=$(ls *.bw)

for FILE in $CONTRIB_FILES; do
    echo "Processing contribution file: $FILE..."
    
    # Use wiggletools to average the files and output in wiggle format
    wiggletools mean \
      $CONTRIB_INPUT_DIR/fold_0/$FILE \
      $CONTRIB_INPUT_DIR/fold_1/$FILE \
      $CONTRIB_INPUT_DIR/fold_2/$FILE \
      $CONTRIB_INPUT_DIR/fold_3/$FILE \
      $CONTRIB_INPUT_DIR/fold_4/$FILE \
      > $CONTRIB_OUTPUT_DIR/${FILE%.bw}.wig
    
    echo "Created wiggle file: $CONTRIB_OUTPUT_DIR/${FILE%.bw}.wig"
    
    # Convert this wiggle file to bigWig (instead of bedGraphToBigWig)
    if [ -f "$chromsizes" ]; then
        WIGGLE=$CONTRIB_OUTPUT_DIR/${FILE%.bw}.wig
        BIGWIG="${WIGGLE%.wig}.bw"
        
        # Convert wiggle to bigWig directly
        wigToBigWig $WIGGLE $chromsizes $BIGWIG
        echo "Created bigWig: $BIGWIG"
        
        # Optionally remove wiggle to save space
        # rm $WIGGLE
    fi
done

# Process prediction scores
echo "Processing prediction scores..."
cd $PRED_INPUT_DIR/fold_0
PRED_FILES=$(ls *.bw)

for FILE in $PRED_FILES; do
    echo "Processing prediction file: $FILE..."
    
    # Use wiggletools to average the files and output in wiggle format
    wiggletools mean \
      $PRED_INPUT_DIR/fold_0/$FILE \
      $PRED_INPUT_DIR/fold_1/$FILE \
      $PRED_INPUT_DIR/fold_2/$FILE \
      $PRED_INPUT_DIR/fold_3/$FILE \
      $PRED_INPUT_DIR/fold_4/$FILE \
      > $PRED_OUTPUT_DIR/${FILE%.bw}.wig
    
    echo "Created wiggle file: $PRED_OUTPUT_DIR/${FILE%.bw}.wig"
    
    # Convert this wiggle file to bigWig (compressed) because the output is in wiggle format (uncompressed)
    if [ -f "$chromsizes" ]; then
        WIGGLE=$PRED_OUTPUT_DIR/${FILE%.bw}.wig
        BIGWIG="${WIGGLE%.wig}.bw"
        
        # Convert wiggle to bigWig directly
        wigToBigWig $WIGGLE $chromsizes $BIGWIG
        echo "Created bigWig: $BIGWIG"
        
        # Optionally remove wiggle to save space
        # rm $WIGGLE
    fi
done
'

# =====================================================================================
# NEW CODE: Convert existing wig files to bigWig using wiggletools and deeptools
# =====================================================================================

echo "Converting existing wig files to bigWig format..."

# Process contribution score wig files
echo "Processing contribution wig files..."
cd $CONTRIB_OUTPUT_DIR
WIG_FILES=$(ls *.wig 2>/dev/null)

if [ -n "$WIG_FILES" ]; then
    for WIG_FILE in $WIG_FILES; do
        echo "Converting $WIG_FILE to bedGraph and then to bigWig..."
        
        # Convert wig to bedGraph using wiggletools
        wiggletools write_bg $CONTRIB_OUTPUT_DIR/temp/${WIG_FILE%.wig}.bedGraph $WIG_FILE
        echo "Created bedGraph: $CONTRIB_OUTPUT_DIR/temp/${WIG_FILE%.wig}.bedGraph"
        
        # Sort bedGraph file (required for conversion)
        sort -k1,1 -k2,2n $CONTRIB_OUTPUT_DIR/temp/${WIG_FILE%.wig}.bedGraph > \
            $CONTRIB_OUTPUT_DIR/temp/${WIG_FILE%.wig}.sorted.bedGraph
        
        # Switch to deeptools environment
        source activate deeptools
        
        # Convert bedGraph to bigWig using deeptools
        bedGraphToBigWig $CONTRIB_OUTPUT_DIR/temp/${WIG_FILE%.wig}.sorted.bedGraph \
            $chromsizes \
            $CONTRIB_OUTPUT_DIR/${WIG_FILE%.wig}.bw
            
        echo "Created bigWig: $CONTRIB_OUTPUT_DIR/${WIG_FILE%.wig}.bw"
        
        # Switch back to wiggletools environment for next file
        source activate wiggletools
    done
else
    echo "No wig files found in $CONTRIB_OUTPUT_DIR"
fi

# Process prediction score wig files
echo "Processing prediction wig files..."
cd $PRED_OUTPUT_DIR
WIG_FILES=$(ls *.wig 2>/dev/null)

if [ -n "$WIG_FILES" ]; then
    for WIG_FILE in $WIG_FILES; do
        echo "Converting $WIG_FILE to bedGraph and then to bigWig..."
        
        # Convert wig to bedGraph using wiggletools
        wiggletools write_bg $PRED_OUTPUT_DIR/temp/${WIG_FILE%.wig}.bedGraph $WIG_FILE
        echo "Created bedGraph: $PRED_OUTPUT_DIR/temp/${WIG_FILE%.wig}.bedGraph"
        
        # Sort bedGraph file (required for conversion)
        sort -k1,1 -k2,2n $PRED_OUTPUT_DIR/temp/${WIG_FILE%.wig}.bedGraph > \
            $PRED_OUTPUT_DIR/temp/${WIG_FILE%.wig}.sorted.bedGraph
        
        # Switch to deeptools environment
        source activate deeptools
        
        # Convert bedGraph to bigWig using deeptools
        bedGraphToBigWig $PRED_OUTPUT_DIR/temp/${WIG_FILE%.wig}.sorted.bedGraph \
            $chromsizes \
            $PRED_OUTPUT_DIR/${WIG_FILE%.wig}.bw
            
        echo "Created bigWig: $PRED_OUTPUT_DIR/${WIG_FILE%.wig}.bw"
        
        # Switch back to wiggletools environment for next file
        source activate wiggletools
    done
else
    echo "No wig files found in $PRED_OUTPUT_DIR"
fi

# Clean up temporary files (optional)
# rm -r $CONTRIB_OUTPUT_DIR/temp
# rm -r $PRED_OUTPUT_DIR/temp

echo "All wig files have been converted to bigWig format."