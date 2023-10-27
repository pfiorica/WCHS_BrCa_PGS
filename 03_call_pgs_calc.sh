#!/bin/bash
#SBATCH --reservation=ubhpc-future
#SBATCH --partition=general-compute --qos=general-compute
#SBATCH --time=12:00:00
#SBATCH --job-name=PGS_calc_Shieh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16000
##SBATCH --requeue
#SBATCH --output=pgs_job_6_Shieh.out
#SBATCH --error=pgs_job_6_Shieh.err
###SBATCH --array=1,2,4,5,6

module load gcc openmpi r-bundle-bioconductor

Rscript UKB.PGS_PNF1.R 6
