#!/bin/sh
#SBATCH --partition=general-compute --qos=general-compute
#SBATCH --time=19:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --constraint=IB
#SBATCH --mem=16000
##SBATCH --requeue
#SBATCH --job-name="01_split_plink"
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

module load plink

for chr in {1..22}; do \

plink --bfile /projects/rpci/wchs/pnfioric/geno_filter_WCHS_Merged_0.01_0.3_AABC_AMBER_FULL/wchs_geno002_full --chr $chr --keep /projects/rpci/wchs/pnfioric/pgs_wchs/3048_updated_genotypes_wchs.fam --make-bed --out /projects/rpci/wchs/pnfioric/geno_filter_WCHS_Merged_0.01_0.3_AABC_AMBER_FULL/gds_files/wchs_3048_chr${chr}; \

done
