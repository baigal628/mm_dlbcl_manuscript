#!/bin/bash
#SBATCH --job-name=cellphoneDB_sick
#SBATCH --mem=128G        # total memory need
#SBATCH -n 16 # 1 core
#SBATCH -N 1 # on one node
#SBATCH --error=/liulab/galib/mouse_scRNAseq_margaret/data/cellphoneDB/lsf_%j_%x.err      # error file
#SBATCH --output=/liulab/galib/mouse_scRNAseq_margaret/data/cellphoneDB/lsf_%j_%x.out      # output file
#SBATCH --mail-type=END,FAIL # email notification when job ends/fails
#SBATCH --mail-user=galib@ds.dfci.harvard.edu # email to notify


cellphonedb method statistical_analysis 10X_format/bcl6_sick_meta.txt 10X_format/bcl6_sick --counts-data hgnc_symbol --output-path /liulab/galib/mouse_scRNAseq_margaret/data/cellphoneDB/output/bcl6_sick_out --threads 16 --subsampling --subsampling-log true --subsampling-num-cells 5000
