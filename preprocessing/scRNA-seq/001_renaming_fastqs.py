#  rename the lanes in your FASTQ files by extracting the flow cell index from the sample names and adding it to the decimal place in the lane identifier

import os
import shutil

# List of input directories containing your FASTQ files
input_directories = [
    "/gchm/scRNAseq/data/fastq_files/flowcell_1/",
    "/gchm/scRNAseq/data/fastq_files/flowcell_2/",
    "/gchm/scRNAseq/data/fastq_files/flowcell_3/",
    "/gchm/scRNAseq/data/fastq_files/flowcell_4/"
]

# Directory to store the copied files
output_directory = "/gchm/scRNAseq/data/fastq_files/renamed"

# Iterate through each input directory
for input_directory in input_directories:
    # List all files in the current input directory
    files = os.listdir(input_directory)

    # Iterate through each file in the input directory
    for input_file in files:
        # Construct the full path for the input file
        input_path = os.path.join(input_directory, input_file)

        # Extract sample number
        sample_number = int(input_file.split("-")[1][0])

        # Extract lane number
        lane_number = int(input_file.split("_")[2][3:])

        # Create new lane number by appending the sample number
        new_lane_number = f"L00{sample_number}{lane_number}"

        # Modify the filename
        new_filename = input_file.replace(f"L00{lane_number}", new_lane_number)
        new_filename = new_filename.replace(f"-{sample_number}-", "-")

        # Construct the full path for the output file
        output_path = os.path.join(output_directory, new_filename)

        # Create a copy with the modified filename
        shutil.copyfile(input_path, output_path)

        # Print the new filename and its corresponding output path
        print(f"Original: {input_file}")
        print(f"Modified: {new_filename}")
        print(f"Output Path: {output_path}")
        print("\n")
