#!/bin/bash
#SBATCH -n 64 #Number of cores
#SBATCH -N 1
##SBATCH --nodelist node12
#SBATCH --mem=84G #Memory per node in MB (see also --mem-per-cpu)
#SBATCH --job-name=14mos10

conda activate MAESTROv1.6.0
Rscript 14mos_inferCNV10.R


