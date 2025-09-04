#!/bin/bash
#SBATCH -n 64 #Number of cores
#SBATCH -N 1
##SBATCH --nodelist node12
#SBATCH --mem=80G #Memory per node in MB (see also --mem-per-cpu)
#SBATCH --job-name=18mos4

source /liulab/galib/miniconda-Py39/bin/activate MAESTROv1.6.0
Rscript run_infer4.R


