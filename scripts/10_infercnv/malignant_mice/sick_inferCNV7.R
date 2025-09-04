getwd()
setwd("/liulab/galib/mouse_scRNAseq_margaret/")
library(rBCS)
library(tidyverse)
library(Seurat)
library(here)
library(devtools)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(infercnv)

B_cell<- readRDS("./data/objects/B_cell/B_cell_annotated.rds")
B_cell  %>% dim()

raw_counts_matrix = GetAssayData(B_cell, "counts")
gene_order_file <- read_tsv(file = "./data/infercnv/gene_order_file.txt", col_names = FALSE)

ann_all <- data.frame(row.names = rownames(B_cell@meta.data),group = B_cell$genotype)
                            
mos <- "sick"

tumor_sample <- B_cell@meta.data  %>% filter(age == mos, genotype != "CD70-/-", genotype != "WT")  %>% 
                    pull(pool_id)  %>% unique()

print(length(tumor_sample))

for (i in tumor_sample[7]){

    ref_cell<- B_cell@meta.data  %>% filter(age == "14mos", genotype == "WT")  %>% rownames()
    tumor_cell<- B_cell@meta.data  %>% filter(pool_id == i)  %>% rownames()

    ann<- ann_all[c(ref_cell,tumor_cell),, drop = FALSE]

    ref_group_names<- c('WT')
                              
                              
    infercnv_obj <- CreateInfercnvObject(
        raw_counts_matrix = raw_counts_matrix,
        gene_order_file = "./data/infercnv/gene_order_file.txt",
        annotations_file = ann,
        ref_group_names = ref_group_names)
                              
    infercnv_obj <- infercnv::run(infercnv_obj, cutoff=0.1,
                                  out_dir = file.path(paste0("./data/infercnv/inferCNV_out_new/", i, "_newsick/")),
                                  cluster_by_groups=TRUE, denoise=TRUE, HMM=TRUE,
                                  num_threads = 64)
}


