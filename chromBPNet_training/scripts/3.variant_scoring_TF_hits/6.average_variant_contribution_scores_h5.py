import h5py
import numpy as np
import glob
import hdf5plugin

def average_h5_files(file_pattern, output_path):
    """
    Average multiple H5 files obtained from calculating variant contribution scores.
    Preserves non-numeric datasets like 'variant_ids' by copying them from the first file.
    """
    # Find all files
    file_paths = glob.glob(file_pattern)
    print(f"Found {len(file_paths)} files:")
    for file_path in file_paths:
        print(f"  - {file_path}")
    
    if not file_paths:
        print("No files found!")
        return
    
    # Create data structures to track what needs to be averaged
    data_sums = {}
    success_counts = {}
    data_shapes = {}
    data_dtypes = {}
    non_numeric_datasets = set()
    
    # First pass: collect all dataset info from all files
    for file_path in file_paths:
        try:
            with h5py.File(file_path, 'r') as f:
                # Recursively collect datasets
                def collect_datasets(name, obj):
                    if isinstance(obj, h5py.Dataset):
                        if name not in data_sums:
                            shape = obj.shape
                            dtype = obj.dtype
                            
                            # Check if this is a non-numeric dataset that should be preserved
                            if dtype == np.dtype('O') or dtype.kind in ['S', 'U']:
                                non_numeric_datasets.add(name)
                                print(f"Detected non-numeric dataset: {name}, will preserve as-is")
                            else:
                                # Use float32 for accumulation of float16 data
                                if dtype == np.float16:
                                    data_sums[name] = np.zeros(shape, dtype=np.float32)
                                else:
                                    data_sums[name] = np.zeros(shape, dtype=dtype)
                                
                                success_counts[name] = 0
                            
                            data_shapes[name] = shape
                            data_dtypes[name] = dtype
                
                # Visit all items
                f.visititems(collect_datasets)
        except Exception as e:
            print(f"Error scanning {file_path}: {e}")
    
    print(f"Found {len(data_sums)} numeric datasets to average")
    print(f"Found {len(non_numeric_datasets)} non-numeric datasets to preserve")
    
    # Process each file for averaging numeric datasets
    for file_path in file_paths:
        try:
            with h5py.File(file_path, 'r') as f:
                # Process each dataset
                for path, data_sum in data_sums.items():
                    try:
                        if path in f:
                            data_sums[path] += f[path][()]
                            success_counts[path] += 1
                            print(f"Successfully read {path} from {file_path}")
                        else:
                            print(f"Dataset {path} not found in {file_path}")
                    except Exception as e:
                        print(f"Error reading {path} from {file_path}: {e}")
        except Exception as e:
            print(f"Error opening {file_path}: {e}")
    
    # Create output file and save averages
    with h5py.File(output_path, 'w') as f_out:
        # First, copy non-numeric datasets from the first file
        with h5py.File(file_paths[0], 'r') as f_in:
            for path in non_numeric_datasets:
                try:
                    # Ensure parent groups exist
                    parts = path.split('/')
                    current_path = ""
                    for i in range(len(parts) - 1):
                        current_path = current_path + parts[i] if i == 0 else current_path + '/' + parts[i]
                        if current_path not in f_out:
                            f_out.create_group(current_path)
                    
                    # Copy dataset
                    data = f_in[path][()]
                    attrs = dict(f_in[path].attrs.items())
                    
                    # Create dataset
                    dataset = f_out.create_dataset(path, data=data)
                    
                    # Copy attributes
                    for key, value in attrs.items():
                        dataset.attrs[key] = value
                    
                    print(f"Copied non-numeric dataset {path} from first file")
                except Exception as e:
                    print(f"Error copying {path} from first file: {e}")
        
        # Then, save averages for numeric datasets
        for path, data_sum in data_sums.items():
            if success_counts[path] > 0:
                # Calculate average
                avg_data = data_sum / success_counts[path]
                
                # Ensure parent groups exist
                parts = path.split('/')
                current_path = ""
                for i in range(len(parts) - 1):
                    current_path = current_path + parts[i] if i == 0 else current_path + '/' + parts[i]
                    if current_path not in f_out:
                        f_out.create_group(current_path)
                
                # Convert back to original dtype if needed
                if data_dtypes[path] == np.float16:
                    avg_data = avg_data.astype(np.float16)
                
                # Get original attributes
                with h5py.File(file_paths[0], 'r') as f_in:
                    attrs = dict(f_in[path].attrs.items())
                
                # Create dataset
                dataset = f_out.create_dataset(path, data=avg_data)
                
                # Copy attributes
                for key, value in attrs.items():
                    dataset.attrs[key] = value
                
                print(f"Saved average of {success_counts[path]} files for {path}")
            else:
                print(f"No successful reads for {path}")
    
    print(f"Averaging complete. Output saved to {output_path}")
    
#if __name__ == "__main__":
#    # Example usage
#    input_pattern = "*.variant_shap.gradientshap.h5"  # Adjust pattern as needed
#    output_file = "averaged_variant_shap.gradientshap.h5"
    
#    average_deepdish_files(input_pattern, output_file)
#    verify_averaged_file(output_file)    input_pattern = "/path/to/shap/files/*.variant_shap.*.h5"
#    output_file = "my_averaged_results.h5"
#    
#    average_deepdish_files(input_pattern, output_file)
#    verify_averaged_file(output_file)

    
# File pattern and output path
file_pattern = "/gchm/cd4_chrombpnet/chrombpnet_model_b7/variant_contribution_scores/perfold/fold_*/cd4_tcells_AJ_common_variants.shap.variant_shap.counts.h5"
output_path = "/gchm/cd4_chrombpnet/chrombpnet_model_b7/variant_contribution_scores/averaged_cd4_tcells_AJ_common_variants.shap.counts.h5"

# Run the averaging function
average_h5_files(file_pattern, output_path)