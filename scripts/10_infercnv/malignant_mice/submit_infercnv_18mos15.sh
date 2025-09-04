#!/bin/bash
#SBATCH -n 64 #Number of cores
#SBATCH -N 1
#SBATCH --mem=80G #Memory per node in MB (see also --mem-per-cpu)
#SBATCH --job-name=18mos15

source /liulab/galib/miniconda-Py39/bin/activate MAESTROv1.6.0
Rscript 18mos_inferCNV15.R


