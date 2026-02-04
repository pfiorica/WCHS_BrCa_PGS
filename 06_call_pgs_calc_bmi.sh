#!/bin/bash
#SBATCH --partition=general-compute --qos=general-compute
#SBATCH --time=12:00:00
#SBATCH --job-name=06_BMI_PGS_Calculation
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16000
##SBATCH --requeue
#SBATCH --output=R-%x.pgs_job_%a_.out
#SBATCH --error=R-%x.pgs_job_%a_.err
#SBATCH --array=1-5

module load gcc openmpi r-bundle-bioconductor

Rscript UKB.PGS_PNF1_BMI2.R $SLURM_ARRAY_TASK_ID
