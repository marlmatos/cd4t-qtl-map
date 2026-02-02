#python script.py features.tsv barcodes.tsv matrix.mtx
import sys

def count_lines(file_path):
    with open(file_path, 'r') as file:
        return sum(1 for line in file)

def count_nonzero_rows(matrix_path):
    num_nonzeros = 0
    with open(matrix_path, 'r') as matrix_file:
        for line in matrix_file:
            if not line.startswith('%'):
                parts = line.split()
                if int(parts[2]) != 0:
                    num_nonzeros += 1
    return num_nonzeros

def generate_matrix_header(features_path, barcodes_path, matrix_path):
    num_rows = count_lines(features_path)
    num_cols = count_lines(barcodes_path)
    num_nonzeros = count_nonzero_rows(matrix_path)

    with open(matrix_path, 'r') as matrix_file:
        matrix_file_lines = matrix_file.readlines()

    with open(matrix_path, 'w') as matrix_file:
        matrix_file.write("%%MatrixMarket matrix coordinate integer general\n")
        matrix_file.write("%\n")
        matrix_file.write(f"{num_rows} {num_cols} {num_nonzeros}\n")
        matrix_file.writelines(matrix_file_lines)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py features_path barcodes_path matrix_path")
        sys.exit(1)

    features_path = sys.argv[1]
    barcodes_path = sys.argv[2]
    matrix_path = sys.argv[3]

    generate_matrix_header(features_path, barcodes_path, matrix_path)
