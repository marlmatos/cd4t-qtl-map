import h5py
import numpy as np
import glob
import hdf5plugin

def average_h5_files(file_pattern, output_path):
    """
    Simple approach to average multiple H5 files based on keys.
    """
    # Find all files
    file_paths = glob.glob(file_pattern)
    print(f"Found {len(file_paths)} files:")
    for file_path in file_paths:
        print(f"  - {file_path}")
    
    if not file_paths:
        print("No files found!")
        return
    
    # Get structure from first file
    with h5py.File(file_paths[0], 'r') as f:
        # Initialize dictionaries to store sums and file counts
        data_sums = {}
        success_counts = {}
        
        # Get all groups and datasets
        for group_name in f.keys():
            for dataset_name in f[group_name].keys():
                key = f"{group_name}/{dataset_name}"
                shape = f[key].shape
                dtype = f[key].dtype
                
                # Use float32 for accumulation of float16 data
                if dtype == np.float16:
                    data_sums[key] = np.zeros(shape, dtype=np.float32)
                else:
                    data_sums[key] = np.zeros(shape, dtype=dtype)
                
                success_counts[key] = 0
    
    # Process each file
    for file_path in file_paths:
        try:
            with h5py.File(file_path, 'r') as f:
                # Process each key
                for key in data_sums.keys():
                    try:
                        data_sums[key] += f[key][()]
                        success_counts[key] += 1
                        print(f"Successfully read {key} from {file_path}")
                    except Exception as e:
                        print(f"Error reading {key} from {file_path}: {e}")
        except Exception as e:
            print(f"Error opening {file_path}: {e}")
    
    # Create output file and save averages
    with h5py.File(output_path, 'w') as f_out:
        for key, data_sum in data_sums.items():
            if success_counts[key] > 0:
                # Calculate average
                avg_data = data_sum / success_counts[key]
                
                # Get group and dataset names
                group_name, dataset_name = key.split('/')
                
                # Create group if needed
                if group_name not in f_out:
                    f_out.create_group(group_name)
                
                # Convert back to original dtype if needed
                original_dtype = h5py.File(file_paths[0], 'r')[key].dtype
                if original_dtype == np.float16:
                    avg_data = avg_data.astype(np.float16)
                
                # Create dataset
                f_out.create_dataset(key, data=avg_data)
                print(f"Saved average of {success_counts[key]} files for {key}")
            else:
                print(f"No successful reads for {key}")
    
    print(f"Averaging complete. Output saved to {output_path}")

# File pattern and output path
file_pattern = "/gchm/cd4_chrombpnet/chrombpnet_model_b7/contribution_scores_bw/fold_*/cd4_tcells.counts_scores.h5"
output_path = "/gchm/cd4_chrombpnet/chrombpnet_model_b7/contribution_scores_bw/averaged_folds_cd4_tcells.counts_scores.h5"

# Run the averaging function
average_h5_files(file_pattern, output_path)