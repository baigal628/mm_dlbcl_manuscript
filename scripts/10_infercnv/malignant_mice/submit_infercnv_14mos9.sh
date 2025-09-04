#!/bin/bash
#SBATCH -n 64 #Number of cores
#SBATCH -N 1
##SBATCH --nodelist node12
#SBATCH --mem=200G #Memory per node in MB (see also --mem-per-cpu)
#SBATCH --job-name=14mos9

source /liulab/galib/miniconda-Py39/bin/activate MAESTROv1.6.0
Rscript 14mos_inferCNV9.R


