#!/bin/bash -l
#SBATCH --reservation=ubhpc-future
#SBATCH --partition=general-compute --qos=general-compute
#SBATCH --time=19:0:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --constraint=IB
#SBATCH --mem=16000
##SBATCH --requeue
#SBATCH --job-name="02_plink2gds"
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

module load gcc
module load openmpi/4.1.1 
module load r-bundle-bioconductor/3.15-R-4.2.0

cd /projects/rpci/wchs/pnfioric/pgs_wchs/pgs_brca

Rscript plink2gds.R
