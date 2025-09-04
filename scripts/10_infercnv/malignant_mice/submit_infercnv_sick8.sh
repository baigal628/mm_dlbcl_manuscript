#!/bin/bash
#SBATCH -n 64 #Number of cores
#SBATCH -N 1
##SBATCH --nodelist node12
#SBATCH --mem=64G #Memory per node in MB (see also --mem-per-cpu)
#SBATCH --job-name=sick8

conda activate MAESTROv1.6.0
Rscript sick_inferCNV8.R


