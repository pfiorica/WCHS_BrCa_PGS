import subprocess
import os
import math
import pandas as pd

# Input directory with PLINK binary files (one directory per chromosome)
input_directory = "/projects/rpci/wchs/pnfioric/geno_filter_WCHS_Merged_0.01_0.3_AABC_AMBER_FULL/split_files"
output_directory = "/projects/rpci/wchs/pnfioric/geno_filter_WCHS_Merged_0.01_0.3_AABC_AMBER_FULL/split_files/subfiles"

# Maximum number of SNPs per subgroup
max_snps_per_subgroup = 250000
# Number of subgroups per chromosome
desired_subgroups_per_chromosome = 10

# Create the output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Process each chromosome directory
for chromosome in range(1, 23):
    chromosome_directory = os.path.join(input_directory, f"wchs_3048_chr{chromosome}")
    output_chromosome_directory = os.path.join(output_directory, f"wchs_3048_chr{chromosome}")
    os.makedirs(output_chromosome_directory, exist_ok=True)

    chromosome_directory = os.path.join(input_directory, f"wchs_3048_chr{chromosome}")
    output_chromosome_directory = os.path.join(output_directory, f"wchs_3048_chr{chromosome}")
    os.makedirs(output_chromosome_directory, exist_ok=True)
    
    # Define the path to your BIM file
    bim_file_path = f"{chromosome_directory}.bim"  # Replace with your actual file path
    # Read the BIM file into a Pandas DataFrame
    bim_df = pd.read_csv(bim_file_path, sep='\t', header=None, names=['chromosome', 'snp_id', 'genetic_distance', 'base_pair_position', 'allele_1', 'allele_2'])
    
    # Calculate the number of subgroups
    num_subgroups = desired_subgroups_per_chromosome
    
    avg_file_size = round(max(bim_df.base_pair_position)/num_subgroups)
    
    start_bp = min(bim_df.base_pair_position)-1
    end_bp = max(bim_df.base_pair_position)
    
    for subgroup_index in range(num_subgroups):
        first_pos = start_bp + (subgroup_index * avg_file_size)
        end_pos = start_bp + ((subgroup_index + 1) * avg_file_size)
            
        plink_command = f"plink --bfile {chromosome_directory} --chr {chromosome} --from-bp {first_pos} --to-bp {end_pos} --make-bed --out {os.path.join(output_chromosome_directory, f'chr{chromosome}.{subgroup_index + 1}')}"
        subprocess.call(plink_command, shell=True)
