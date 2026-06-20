
#######################
# N cells per cluster #
#######################


get_cell_number<- function(obj, space, seurat_clusters){
  number_of_cells_per_cluster<- obj@meta.data %>%
    count(sample_id, age, genotype, seurat_clusters, .drop = FALSE)
  
  df<- left_join(number_of_cells_per_cluster, total_number_of_cells_per_sample) %>%
    mutate(percent = n/total_number)
  
  write_tsv(df, here(glue::glue("results/tsv/{space}/{space}_number_of_cells_per_cluster.tsv")))
  
  number_of_cells_per_cluster_wider<- df  %>% 
    pivot_wider(names_from = seurat_clusters, values_from = n)
  
  write_tsv(df, here(glue::glue("results/tsv/{space}/{space}_number_of_cells_per_cluster_wide.tsv")))
  
  p1<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = age, y =n)) +
    geom_boxplot(aes(fill = age), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    facet_grid(seurat_clusters ~ genotype, scale= "free_y")
  
  ggsave(filename = here(glue::glue("results/figures/{space}/final_annotation/{space}_number_of_cells_per_cluster.pdf")),
         plot = p1, width = 8, height = 50,
         limitsize = FALSE)
  
  p2<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = age, y =percent)) +
    scale_y_continuous(labels = scales::percent) +
    geom_boxplot(aes(fill = age), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    facet_grid(seurat_clusters ~ genotype, scale= "free_y")
  
  ggsave(filename = here(glue::glue("results/figures/{space}/final_annotation/{space}_percent_of_cells_per_cluster.pdf")), plot = p2, width = 8, height = 50,
         limitsize = FALSE)
}


#######################
# Get Remain Clusters #
#######################

#>41 cells with >= 3 samples

cluster_remain_41_exclusive<- function(obj, space, seurat_clusters){
  number_of_cells_per_cluster<- obj@meta.data %>% count(sample_id, age, genotype, seurat_clusters, .drop = FALSE)
  
  cluster_remain<- c()
  for (i in 0:(length(number_of_cells_per_cluster$seurat_clusters  %>% unique) -1)){
    n_sample<- number_of_cells_per_cluster  %>% filter(seurat_clusters == i, n > 41)  %>% pull(sample_id)  %>% length()
    if (n_sample >= 3){
      cluster_remain<- c(cluster_remain, i)
    }
  }
  print(space)
  print(cluster_remain)
  
  number_of_cells_per_cluster_wider<- number_of_cells_per_cluster  %>%
    pivot_wider(names_from = seurat_clusters, values_from = n)
  
  number_of_cells_per_cluster_wider_remain<- number_of_cells_per_cluster  %>%
    filter(seurat_clusters %in% cluster_remain)  %>%
    pivot_wider(names_from = seurat_clusters, values_from = n)
  write_tsv(number_of_cells_per_cluster_wider, here("./results/tsv/", paste0(space, "/", space, "_number_of_cells_all_clusters_wide.tsv")))
  write_tsv(number_of_cells_per_cluster_wider_remain, here("./results/tsv/", paste0(space, "/", space, "_number_of_cells_remain_clusters_wide.tsv")))
  
  return(cluster_remain)
}



###################
# Remove Clusters #
###################

remove_clusters<- function(obj, obj_clusters, space, space_name, seurat_clusters, 
                           outdir = "./", fsuffix = "clusters_removed"){
  
  print(paste0("Export data: ", outdir, "data/objects/", space, "/", space, "_",  fsuffix ,".rds"))
  print(paste0("Export file: ", outdir, "results/figures/", space, "/", space, "_",  fsuffix ,".pdf"))
  # Remove clusters
  obj$filter <- ifelse(obj$seurat_clusters %in% obj_clusters, "Keep", "Remove")
  
  obj_removed<- subset(x = obj, subset = filter == "Keep")
  
  obj_removed$seurat_clusters<- factor(x = obj_removed$seurat_clusters,
                                       levels = obj_clusters)
  
  P_umap<- DimPlot(obj_removed, reduction = "umap",
                   label = TRUE, pt.size = 0.2 ) + 
    labs(title = NULL, y = NULL, x = NULL) +
    theme(text = element_text(size = 20))
  
  #return(P_umap)
  #Save rds
  
  ggsave(paste0(outdir, "results/figures/", space, "/", space, "_",  fsuffix ,"_umap.pdf"), P_umap, width = 12, height = 8)
  
  saveRDS(obj_removed, paste0(outdir, "data/objects/", space, "/", space, "_",  fsuffix ,".rds"))
  
  ExportSeurat(obj_removed, paste0(outdir, "data/objects/", space, "/", space, "_", fsuffix ,".bcs"), overwrite=TRUE)
  
}



########################
# Bi-Clustered Heatmap #
########################


plot_bi_clustered_heatmap <- function(obj, ident = "seurat_clusters", slot = "data",
                                      space, outdir = "./", fsuffix = "clusters_removed"){
  set.seed(123)
  clusters <- obj@meta.data %>% pull(ident) %>% unique
  markers <- presto::wilcoxauc(obj, ident , assay = slot,
                               groups_use = clusters)
  write_tsv(markers,
            paste0(outdir, "data/", space, "_", fsuffix, "_",  fsuffix ,"_all_markers.tsv"))
  
  
  top_markers<- markers  %>% filter(padj <=0.05, abs(logFC) >=0.8)
  
  write_tsv(top_markers,
            paste0(outdir, "data/", space, "_",  fsuffix ,"_top_markers(padj<=0.05, abs(logFC) >= 0.8.tsv"))
  
  top_markers_unique <- top_markers %>% pull(feature) %>% unique()
  
  # top_10_markers <- top_markers(markers, n = 10, auc_min = .6,
  #                               pct_in_min = 20, pct_out_max = 30)
  # top10 <- top_10_markers[,-1]  %>% unlist()  %>% unique() %>% na.omit()
  
  
  # top_5_markers <- top_markers(markers, n = 5, auc_min = .6,
  #                              pct_in_min = 20, pct_out_max = 30)
  # 
  # top5 <- top_5_markers[,-1]  %>% unlist()  %>% unique() %>% na.omit()
  
  
  top_10_markers<- top_markers %>% group_by(group) %>% arrange(padj, -abs(logFC), .by_group = TRUE)  %>% slice_head(n = 10)
  
  top10<- top_10_markers  %>% pull(feature)  %>% unique()
  
  write_tsv(top_10_markers,
            paste0(outdir, "data/", space, "_",  fsuffix ,"_top10_markers(padj<=0.05, abs(logFC) >= 0.8.tsv"))
  
  top_5_markers<- top_markers %>% group_by(group) %>% arrange(padj, -abs(logFC), .by_group = TRUE)  %>% slice_head(n = 5)
  
  top5<- top_5_markers  %>% pull(feature)  %>% unique()
  
  write_tsv(top_5_markers,
            paste0(outdir, "data/", space, "_",  fsuffix ,"_top5_markers(padj<=0.05, abs(logFC) >= 0.8.tsv"))
  
  
  ## Built data for plotting
  obj_sub <- obj[, sample.int(ncol(obj), size = 10000)]
  
  all_cells <-  obj_sub@meta.data$seurat_clusters %in% clusters
  
  all_mat<- obj_sub[["RNA"]]@data[top_markers_unique, all_cells] %>% as.matrix()
  top_5_mat<- obj_sub[["RNA"]]@data[top5, all_cells] %>% as.matrix()
  top_10_mat<- obj_sub[["RNA"]]@data[top10, all_cells] %>% as.matrix()
  
  ## scale the rows
  all_mat<- t(scale(t(all_mat)))
  top_5_mat<- t(scale(t(top_5_mat)))
  top_10_mat<- t(scale(t(top_10_mat)))
  
  ## Heatmap annotations
  #cluster_anno<- obj_sub@meta.data$seurat_clusters[all_cells]
  cluster_id <- obj_sub$seurat_clusters[all_cells]
  genotype <-  obj_sub$genotype[all_cells]
  age <- obj_sub$age[all_cells]
  len <- cluster_id  %>% sort()  %>% unique()  %>% length()
  
  genotype <- factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+"))
  age <- factor(age, levels= c("6mos", "14mos", "18mos", "sick"))
  
  reorder_df<- data.frame(cluster_id = cluster_id,
                          genotype = genotype,
                          age = age)
  
  reorder_df<- reorder_df %>% mutate(ind = row_number()) %>%
    arrange(cluster_id, genotype, age)
  
  P_cluster = createPalette(len,  c("#ff0000", "#00ff00", "#0000ff"))
  P4 <- scales::hue_pal()(4)
  P4.2 = brewer.pal(4, "Dark2")
  
  column_ha<- HeatmapAnnotation(cluster_id = cluster_id[reorder_df$ind],
                                genotype = genotype[reorder_df$ind],
                                age = age[reorder_df$ind],
                                col = list(cluster_id = setNames(P_cluster, levels(cluster_id)),
                                           genotype = setNames(P4, levels(genotype)),
                                           age = setNames(P4.2, levels(age))),
                                na_col = "grey")
  column_split_df = data.frame(cluster_id = cluster_id[reorder_df$ind],
                               genotype = genotype[reorder_df$ind])
  
  column_ha_cluster<- HeatmapAnnotation(
    cluster_id = cluster_id,
    col = list(cluster_id = setNames(P_cluster, levels(cluster_id))),
    na_col = "grey",
    show_annotation_name = TRUE)
  
  
  ## Colors
  color_range_all_low = quantile(all_mat, 0.1)
  color_range_top5_low = quantile(top_5_mat, 0.1)
  color_range_top10_low = quantile(top_10_mat, 0.1)
  
  color_range_all_high = quantile(all_mat, 0.95)
  color_range_top5_high = quantile(top_5_mat, 0.95)
  color_range_top10_high = quantile(top_10_mat, 0.95)
  
  # col_fun_all = circlize::colorRamp2(c(color_range_all_low,0,1), viridis(3))
  # col_fun_top5 = circlize::colorRamp2(c(color_range_top5_low, 0, 1), viridis(3))
  # col_fun_top10 = circlize::colorRamp2(c(color_range_top10_low, 0, 1), viridis(3))
  
  col_fun_all = circlize::colorRamp2(c(-1,0,2), c("#440154FF", "#238A8DFF", "#FDE725FF"))
  col_fun_top5 = circlize::colorRamp2(c(-1, 0, 2), c("#440154FF", "#238A8DFF", "#FDE725FF"))
  col_fun_top10 = circlize::colorRamp2(c(-1, 0, 2), c("#440154FF", "#238A8DFF", "#FDE725FF"))
  
  
  p1 = Heatmap(all_mat, name = "Expression",
               column_split = factor(cluster_id),
               top_annotation = column_ha_cluster,
               cluster_columns = TRUE,
               show_column_dend = FALSE,
               show_column_names = FALSE,
               cluster_column_slices = TRUE,
               cluster_rows = TRUE,
               show_row_dend = FALSE,
               column_title_gp = gpar(fontsize = 6),
               column_title_rot = 90,
               row_km = 20,
               column_gap = unit(0.5, "mm"),
               border = NA,
               row_names_gp = gpar(fontsize = 5),
               col = col_fun_all,
               use_raster = TRUE,
               raster_by_magick = FALSE,
               raster_quality = 8)
  
  p2 = Heatmap(top_5_mat, name = "Expression",
               column_split = factor(cluster_id),
               top_annotation = column_ha_cluster,
               cluster_columns = TRUE,
               show_column_dend = FALSE,
               show_column_names = FALSE,
               cluster_column_slices = TRUE,
               cluster_rows = TRUE,
               show_row_dend = FALSE,
               column_title_gp = gpar(fontsize = 5),
               column_gap = unit(0.3, "mm"),
               border = NA,
               row_names_gp = gpar(fontsize = 5),
               col = col_fun_top5,
               use_raster = TRUE,
               raster_by_magick = FALSE,
               raster_quality = 8)
  
  p3 = Heatmap(top_10_mat, name = "Expression",
               column_split = factor(cluster_id),
               top_annotation = column_ha_cluster,
               cluster_columns = TRUE,
               show_column_dend = FALSE,
               show_column_names = FALSE,
               cluster_column_slices = TRUE,
               cluster_rows = TRUE,
               show_row_dend = FALSE,
               column_title_gp = gpar(fontsize = 5),
               column_gap = unit(0.3, "mm"),
               border = NA,
               row_names_gp = gpar(fontsize = 5),
               col = col_fun_top10,
               use_raster = TRUE,
               raster_by_magick = FALSE,
               raster_quality = 8)
  return(list(p1 = p1, p2 = p2, p3 = p3))
}


#########
# Anova #
#########

get_cell_number_between_age<- function(obj, space, vjust = 3, seurat_clusters){
  number_of_cells_per_cluster<- obj@meta.data %>%
    count(sample_id, age, genotype, seurat_clusters, .drop = FALSE)
  
  df<- left_join(number_of_cells_per_cluster, total_number_of_cells_per_sample) %>%
    mutate(percent = n/total_number)
  
  p1<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = age, y =n)) +
    geom_boxplot(aes(fill = age), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    facet_grid(seurat_clusters ~ genotype, scale= "free_y") +
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)), vjust = vjust) +
    theme_bw()
  
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_number_of_cells_per_cluster_between_age_anova.pdf")),
         plot = p1, width = 10, height = 30, limitsize = FALSE)
  
  p2<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = age, y =percent)) +
    scale_y_continuous(labels = scales::percent) +
    geom_boxplot(aes(fill = age), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    facet_grid(seurat_clusters ~ genotype, scale= "free_y") +
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)), vjust = vjust) +
    theme_bw()
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_percent_of_cells_per_cluster_between_age_anova.pdf")),
         plot = p2, width = 10, height = 30, limitsize = FALSE)
}


get_cell_number_between_genotype<- function(obj, space, vjust= 3, seurat_clusters){
  number_of_cells_per_cluster<- obj@meta.data %>%
    count(sample_id, age, genotype, seurat_clusters, .drop = FALSE)
  
  df<- left_join(number_of_cells_per_cluster, total_number_of_cells_per_sample) %>%
    mutate(percent = n/total_number)
  
  p1<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = genotype, y =n)) +
    geom_boxplot(aes(fill = genotype), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    facet_grid(seurat_clusters ~ age, scale= "free_y") +
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)), vjust = vjust) +
    theme_bw()
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_number_of_cells_per_cluster_between_genotype_anova.pdf")),
         plot = p1, width = 10, height = 30, limitsize = FALSE)
  
  p2<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = genotype, y =percent)) +
    scale_y_continuous(labels = scales::percent) +
    geom_boxplot(aes(fill = genotype), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    facet_grid(seurat_clusters ~ age, scale= "free_y") +
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)), vjust = vjust) +
    theme_bw()
  
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_percent_of_cells_per_cluster_between_genotype_anova.pdf")),
         plot = p2, width = 10, height = 30, limitsize = FALSE)
}


##############################
# Custom Comparison with pva #
##############################
###############
# By Genotype #
###############

genotype_comparisons <- list(c("Bcl6tg/+", "CD70-/-;Bcl6tg/+"),
                             c("WT", "Bcl6tg/+"),
                             c("CD70-/-", "CD70-/-;Bcl6tg/+"),
                             c("WT", "CD70-/-;Bcl6tg/+"))

get_cell_number_between_genotype_custom<- function(obj, space, vjust = 0, seurat_clusters){
  number_of_cells_per_cluster<- obj@meta.data %>%
    count(sample_id, age, genotype, seurat_clusters, .drop = FALSE)
  
  df<- left_join(number_of_cells_per_cluster, total_number_of_cells_per_sample) %>%
    mutate(percent = n/total_number)
  
  p1<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = genotype, y = n)) +
    geom_boxplot(aes(fill = genotype), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(seurat_clusters ~ age, scale= "free_y") +
    stat_compare_means(comparisons = genotype_comparisons, method= "wilcox.test", vjust = vjust) +
    labs(x = "Genotypes", y = paste0("Number of cells per cluster "))
  
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_number_of_cells_per_cluster_between_genotype_wilcox.pdf")),
         plot = p1, width = 10, height = 35, limitsize = FALSE)
  
  p2<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = genotype, y = percent)) +
    geom_boxplot(aes(fill = genotype), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(seurat_clusters ~ age, scale= "free_y") +
    stat_compare_means(comparisons = genotype_comparisons, method= "wilcox.test", vjust = vjust) +
    labs(x = "Genotypes", y = paste0("Percent of cells per cluster "))
  
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_percent_of_cells_per_cluster_between_genotype_wilcox.pdf")),
         plot = p2, width = 10, height = 35, limitsize = FALSE)
}


##########
# By Age #
##########
age_comparisons <- list( c("6mos", "14mos"),
                         c("14mos", "18mos"),
                         c("6mos", "18mos"))

get_cell_number_between_age_custom<- function(obj, space, vjust = 0, seurat_clusters){
  number_of_cells_per_cluster<- obj@meta.data %>%
    count(sample_id, age, genotype, seurat_clusters, .drop = FALSE)
  
  df<- left_join(number_of_cells_per_cluster, total_number_of_cells_per_sample) %>%
    mutate(percent = n/total_number)
  
  p1<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = age, y =n)) +
    geom_boxplot(aes(fill = age), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(seurat_clusters ~ genotype, scale= "free_y") +
    stat_compare_means(comparisons = age_comparisons, method= "wilcox.test", vjust = vjust) +
    labs(x = "Age", y = paste0("Number of cells per cluster"))
  
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_number_of_cells_per_cluster_between_age_wilcox.pdf")),
         plot = p1, width = 10, height = 35, limitsize = FALSE)
  
  p2<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = age, y = percent)) +
    geom_boxplot(aes(fill = age), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(seurat_clusters ~ genotype, scale= "free_y") +
    stat_compare_means(comparisons = age_comparisons, method= "wilcox.test", vjust = vjust) +
    labs(x = "Age", y = paste0("Percent of cells per cluster"))
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_percent_of_cells_per_cluster_between_age_wilcox.pdf")),
         plot = p2, width = 10, height = 35, limitsize = FALSE)
}

##############################
# Custom Comparison with sig #
##############################

###############
# By Genotype #
###############

genotype_comparisons <- list(c("Bcl6tg/+", "CD70-/-;Bcl6tg/+"),
                             c("WT", "Bcl6tg/+"),
                             c("CD70-/-", "CD70-/-;Bcl6tg/+"),
                             c("WT", "CD70-/-;Bcl6tg/+"))

get_cell_number_between_genotype_custom_sig<- function(obj, space, vjust = 0.5, seurat_clusters){
  number_of_cells_per_cluster<- obj@meta.data %>%
    count(sample_id, age, genotype, seurat_clusters, .drop = FALSE)
  
  df<- left_join(number_of_cells_per_cluster, total_number_of_cells_per_sample) %>%
    mutate(percent = n/total_number)
  
  p1<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = genotype, y = n)) +
    geom_boxplot(aes(fill = genotype), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(seurat_clusters ~ age, scale= "free_y") +
    stat_compare_means(comparisons = genotype_comparisons, method= "wilcox.test", label = "p.signif",
                       vjust = vjust,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.5),
                                          symbols = c("****", "***", "**", "*", "nsig"))) +
    labs(x = "Genotypes", y = paste0("Number of cells per cluster "))
  
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_number_of_cells_per_cluster_between_genotype_wilcox_symbol.pdf")),
         plot = p1, width = 10, height = 35, limitsize = FALSE)
  
  p2<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = genotype, y = percent)) +
    geom_boxplot(aes(fill = genotype), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(seurat_clusters ~ age, scale= "free_y") +
    stat_compare_means(comparisons = genotype_comparisons, method= "wilcox.test", label = "p.signif",
                       vjust = vjust,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.5),
                                          symbols = c("****", "***", "**", "*", "nsig"))) +
    labs(x = "Genotypes", y = paste0("Percent of cells per cluster "))
  
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_percent_of_cells_per_cluster_between_genotype_wilcox_symbol.pdf")),
         plot = p2, width = 10, height = 35, limitsize = FALSE)
}


##########
# By Age #
##########
age_comparisons <- list( c("6mos", "14mos"),
                         c("14mos", "18mos"),
                         c("6mos", "18mos"))

get_cell_number_between_age_custom_sig<- function(obj, space, vjust = 0.5, seurat_clusters){
  number_of_cells_per_cluster<- obj@meta.data %>%
    count(sample_id, age, genotype, seurat_clusters, .drop = FALSE)
  
  df<- left_join(number_of_cells_per_cluster, total_number_of_cells_per_sample) %>%
    mutate(percent = n/total_number)
  
  p1<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = age, y =n)) +
    geom_boxplot(aes(fill = age), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(seurat_clusters ~ genotype, scale= "free_y") +
    stat_compare_means(comparisons = age_comparisons, method= "wilcox.test", label = "p.signif",
                       vjust = vjust,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.5),
                                          symbols = c("****", "***", "**", "*", "nsig"))) +
    labs(x = "Age", y = paste0("Number of cells per cluster"))
  
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_number_of_cells_per_cluster_between_age_wilcox_symbol.pdf")),
         plot = p1, width = 10, height = 35, limitsize = FALSE)
  
  p2<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = age, y = percent)) +
    geom_boxplot(aes(fill = age), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(seurat_clusters ~ genotype, scale= "free_y") +
    stat_compare_means(comparisons = age_comparisons, method= "wilcox.test", label = "p.signif",
                       vjust = vjust,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.5),
                                          symbols = c("****", "***", "**", "*", "nsig"))) +
    labs(x = "Age", y = paste0("Percent of cells per cluster"))
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_percent_of_cells_per_cluster_between_age_wilcox_symbol.pdf")),
         plot = p2, width = 10, height = 35, limitsize = FALSE)
}

# Extract scaled matrix
GetMatrixFromSeurat<- function(obj, features, group, ...) {
    p<- DotPlot(object = obj, features = features, group.by = group, ...)
    df<- p$data
  
    exp_mat<-df %>% 
        select(-pct.exp, -avg.exp) %>%  
        arrange(id) %>%
        pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
        as.data.frame() 

    row.names(exp_mat) <- exp_mat$features.plot  
    exp_mat <- exp_mat[,-1] %>% as.matrix()
    
    percent_mat<-df %>% 
        select(-avg.exp, -avg.exp.scaled) %>%  
        arrange(id) %>%
        pivot_wider(names_from = id, values_from = pct.exp) %>% 
        as.data.frame()
  
   row.names(percent_mat) <- percent_mat$features.plot  
   percent_mat <- percent_mat[,-1] %>% as.matrix()
   
    
   percent_mat<- percent_mat  %>% replace(is.na(.),0)
   exp_mat <- exp_mat %>% replace(is.na(.),0)
   percent_mat <- percent_mat[complete.cases(exp_mat), ]
   exp_mat <- exp_mat[complete.cases(exp_mat), ] 


   if (!identical(dim(exp_mat), dim(percent_mat))) {
     stop("the dimension of the two matrice should be the same!")
   }
  
   if(! all.equal(colnames(exp_mat), colnames(percent_mat))) {
     stop("column names of the two matrice should be the same!")
   }
  
   if(! all.equal(rownames(exp_mat), rownames(percent_mat))) {
     stop("column names of the two matrice should be the same!")
   }
  
   return(list(exp_mat = exp_mat, percent_mat = percent_mat))
  
}

# Extract unscaled matrix
GetMatrixFromSeurat2<- function (obj, features, group, ...) 
{
  p <- DotPlot(object = obj, features = features, group.by = group, 
               ...)
  df <- p$data
  exp_mat <- df %>% select(-pct.exp, -avg.exp.scaled) %>% arrange(id) %>% 
    pivot_wider(names_from = id, values_from = avg.exp) %>% 
    as.data.frame()
  row.names(exp_mat) <- exp_mat$features.plot
  exp_mat <- exp_mat[, -1] %>% as.matrix()
  percent_mat <- df %>% select(-avg.exp, -avg.exp.scaled) %>% 
    arrange(id) %>% pivot_wider(names_from = id, values_from = pct.exp) %>% 
    as.data.frame()
  row.names(percent_mat) <- percent_mat$features.plot
  percent_mat <- percent_mat[, -1] %>% as.matrix()
  percent_mat <- percent_mat[complete.cases(exp_mat), ]
  exp_mat <- exp_mat[complete.cases(exp_mat), ]
  print(ncol(exp_mat))
  
  if (!identical(dim(exp_mat), dim(percent_mat))) {
    stop("the dimension of the two matrice should be the same!")
  }
  if (!all.equal(colnames(exp_mat), colnames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }
  if (!all.equal(rownames(exp_mat), rownames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }
  return(list(exp_mat = exp_mat, percent_mat = percent_mat))
}

MakeClusterDotPlot<- function(exp_mat, percent_mat, col_fun, legend_title = "expression",
                              column_title = "Clustered Dotplot",
                              row_names_font_size = 3, column_ha, cluster_id = NULL,
                              legend_dot_size = unit(2,"mm"), cluster_rows=TRUE, cluster_columns = TRUE, ...){
  
  if (!is.null(cluster_id)){
    exp_mat <- exp_mat[,cluster_id]
    colnames(exp_mat) <- gsub(".*_", "", cluster_id)
    
    percent_mat <- percent_mat[,cluster_id]
    colnames(percent_mat) <- gsub(".*_", "", cluster_id)
  }
    
    layer_fun = function(j, i, x, y, w, h, fill){
    grid.rect(x = x, y = y, width = w, height = h, 
              gp = gpar(col = NA, fill = NA))
    grid.circle(x=x,y=y,r= sqrt(pindex(percent_mat, i, j)/100) * legend_dot_size,
                gp = gpar(fill = col_fun(pindex(exp_mat, i, j)), col = NA))}
  
  hp<- Heatmap(exp_mat,
               heatmap_legend_param=list(title= legend_title),
               column_title = column_title, 
               col=col_fun,
               cluster_rows = cluster_rows,
               cluster_columns = cluster_columns,
               rect_gp = gpar(type = "none"),
               layer_fun = layer_fun,
               row_names_gp = gpar(fontsize = row_names_font_size),
               border = "black",
               top_annotation = column_ha,
               ...)
  # see https://github.com/jokergoo/ComplexHeatmap/issues/737 for retain k-means
  
  return(hp)
}


#######################
###inferCNV analysis###
#######################


concat_mtx1<- function(poolid_list, obj = "/run.final.infercnv_obj", inpath){
    
    # reference cells are keeped in the very first object
    # (i.e. assuming reference cells are shared among samples in the pool_ids)
    
    inferObj<- readRDS(paste0(inpath, poolid_list[1], obj))
    concat_mtx<- inferObj@expr.data

    for (i in poolid_list[-1]){
        inferObj<- readRDS(paste0(inpath, i, obj))
        obs_cells<- unlist(inferObj@observation_grouped_cell_indices)
        mtx_obs<- inferObj@expr.data[,obs_cells]
        geneList1<- rownames(concat_mtx)
        geneList2<- rownames(mtx_obs)
        shared_gene<- intersect(geneList1, geneList2)
        concat_mtx = cbind(concat_mtx[shared_gene,], mtx_obs[shared_gene,])
        message(paste0("Done adding sample ", i))
    }

    return(concat_mtx)
}


concat_mtx2<- function(poolid_list, obj = "/run.final.infercnv_obj", inpath){
    
    #WT from 14mos is used as ref
    inferObj<- readRDS(paste0(inpath, poolid_list[1], obj))
    obs_cells<- unlist(inferObj@observation_grouped_cell_indices)
    concat_mtx<- inferObj@expr.data[,obs_cells]

    for (i in poolid_list[-1]){
        inferObj<- readRDS(paste0(inpath, i, obj))
        obs_cells<- unlist(inferObj@observation_grouped_cell_indices)
        mtx_obs<- inferObj@expr.data[,obs_cells]
        geneList1<- rownames(concat_mtx)
        geneList2<- rownames(mtx_obs)
        shared_gene<- intersect(geneList1, geneList2)
        concat_mtx = cbind(concat_mtx[shared_gene,], mtx_obs[shared_gene,])
        message(paste0("Done adding sample ", i))
    }

    return(concat_mtx)
}

###########################
###cellphoneDB analysis####
###########################
build_cpdb_mtx<- function(cols, scaled = FALSE, pval = 0.01, 
                          sig.cellpairs= 3,
                          sig.timepoint = 2,
                          group_name = "",
                          ph_6mos_pvals, ph_6mos_means,
                          ph_14mos_pvals, ph_14mos_means, 
                          ph_18mos_pvals, ph_18mos_means,
                          ph_sick_pvals, ph_sick_means){
    
    ph_6mos_pvals_sub<- ph_6mos_pvals[,cols]
    interactions_6mos = ph_6mos_pvals[rowSums(ph_6mos_pvals_sub < pval) >= sig.cellpairs, 'interacting_pair']

    ph_14mos_pvals_sub<- ph_14mos_pvals[,cols]
    interactions_14mos = ph_14mos_pvals[rowSums(ph_14mos_pvals_sub < pval) >= sig.cellpairs, 'interacting_pair']

    ph_18mos_pvals_sub<- ph_18mos_pvals[,cols]
    interactions_18mos = ph_18mos_pvals[rowSums(ph_18mos_pvals_sub < pval) >= sig.cellpairs, 'interacting_pair']

    ph_sick_pvals_sub<- ph_sick_pvals[,cols]
    interactions_sick = ph_sick_pvals[rowSums(ph_sick_pvals_sub < pval) >= sig.cellpairs, 'interacting_pair']

    interaction_df = table(c(interactions_6mos, interactions_14mos, interactions_18mos, interactions_sick))
    interactions = names(interaction_df)[interaction_df>= sig.timepoint]
    # interactions<- c(interactions_6, interactions_14, interactions_18, interactions_sick)  %>% sort()  %>% unique()
    
    message("total number of intersections considered: ", interactions  %>% length())
    
    message("6 mos...")
    ph_6mos_pvals_dedup = ph_6mos_pvals[!ph_6mos_pvals$interacting_pair  %>% duplicated(), ]
    ph_6mos_means_dedup = ph_6mos_means[!ph_6mos_means$interacting_pair  %>% duplicated(), ]
    
    rownames(ph_6mos_pvals_dedup) = ph_6mos_pvals_dedup$interacting_pair
    matchID= match(interactions, rownames(ph_6mos_pvals_dedup))

    ph_6mos_mtx= ph_6mos_means_dedup[matchID, cols]
    rownames(ph_6mos_mtx) = interactions
    ph_6mos_mtx[is.na(ph_6mos_mtx)] = 0
    ph_6mos_mtx  %>% dim()

    ph_6mos_mtxp= ph_6mos_pvals_dedup[matchID, cols]
    rownames(ph_6mos_mtxp) = interactions
    ph_6mos_mtxp[is.na(ph_6mos_mtxp)] = 1
    ph_6mos_mtxp  %>% dim()
    
    message("14 mos...")
    ph_14mos_pvals_dedup = ph_14mos_pvals[!ph_14mos_pvals$interacting_pair  %>% duplicated(), ]
    ph_14mos_means_dedup = ph_14mos_means[!ph_14mos_means$interacting_pair  %>% duplicated(), ]
    
    rownames(ph_14mos_pvals_dedup) = ph_14mos_pvals_dedup$interacting_pair
    matchID= match(interactions, rownames(ph_14mos_pvals_dedup))

    ph_14mos_mtx= ph_14mos_means_dedup[matchID, cols]
    rownames(ph_14mos_mtx) = interactions
    ph_14mos_mtx[is.na(ph_14mos_mtx)] = 0
    ph_14mos_mtx  %>% dim()

    ph_14mos_mtxp= ph_14mos_pvals_dedup[matchID, cols]
    rownames(ph_14mos_mtxp) = interactions
    ph_14mos_mtxp[is.na(ph_14mos_mtxp)] = 1
    ph_14mos_mtxp  %>% dim()
    
    message("18 mos...")
    ph_18mos_pvals_dedup = ph_18mos_pvals[!ph_18mos_pvals$interacting_pair  %>% duplicated(), ]
    ph_18mos_means_dedup = ph_18mos_means[!ph_18mos_means$interacting_pair  %>% duplicated(), ]
    
    rownames(ph_18mos_pvals_dedup) = ph_18mos_pvals_dedup$interacting_pair
    matchID= match(interactions, rownames(ph_18mos_pvals_dedup))

    ph_18mos_mtx= ph_18mos_means_dedup[matchID, cols]
    rownames(ph_18mos_mtx) = interactions
    ph_18mos_mtx[is.na(ph_18mos_mtx)] = 0
    ph_18mos_mtx  %>% dim()

    ph_18mos_mtxp= ph_18mos_pvals_dedup[matchID, cols]
    rownames(ph_18mos_mtxp) = interactions
    ph_18mos_mtxp[is.na(ph_18mos_mtxp)] = 1
    ph_18mos_mtxp  %>% dim()
    
    message("sick...")
    ph_sick_pvals_dedup = ph_sick_pvals[!ph_sick_pvals$interacting_pair  %>% duplicated(), ]
    ph_sick_means_dedup = ph_sick_means[!ph_sick_means$interacting_pair  %>% duplicated(), ]
    
    rownames(ph_sick_pvals_dedup) = ph_sick_pvals_dedup$interacting_pair
    matchID= match(interactions, rownames(ph_sick_pvals_dedup))
    
    ph_sick_mtx= ph_sick_means_dedup[matchID, cols]
    rownames(ph_sick_mtx) = interactions
    ph_sick_mtx[is.na(ph_sick_mtx)] = 0
    ph_sick_mtx  %>% dim()

    ph_sick_mtxp= ph_sick_pvals_dedup[matchID, cols]
    rownames(ph_sick_mtxp) = interactions
    ph_sick_mtxp[is.na(ph_sick_mtxp)] = 1
    ph_sick_mtxp  %>% dim()
    
    message("Generating count matrix...")
    combined_mtx = cbind(ph_6mos_mtx, ph_14mos_mtx, ph_18mos_mtx, ph_sick_mtx)
    combined_mtx_scaled = t(apply(combined_mtx, 1, scale))
    
    group_n = paste0(c("06mos", "14mos", "18mos", "sick"), group_name)
    
    colnames(combined_mtx_scaled) = paste0(rep(colnames(ph_6mos_mtx), times =4), "_", rep(group_n, each=ncol(ph_6mos_mtx)))
    colnames(combined_mtx) = paste0(rep(colnames(ph_6mos_mtx), times =4), "_", rep(group_n, each=ncol(ph_6mos_mtx)))
    
    ordered_name = colnames(combined_mtx_scaled)  %>% sort()
    
    combined_mtx_scaled = combined_mtx_scaled[,ordered_name]
    combined_mtx = combined_mtx[,ordered_name]
    
    message("Generating pvalue matrix...")
    pval_mtx = cbind(ph_6mos_mtxp, ph_14mos_mtxp, ph_18mos_mtxp, ph_sick_mtxp)
    pval_mtx[pval_mtx == 0] = 1e-5
    pval_mtx = -log10(pval_mtx)
    
    colnames(pval_mtx) = paste0(rep(colnames(ph_6mos_mtxp), times =4), "_", rep(group_n, each=ncol(ph_6mos_mtxp)))
    pval_mtx = pval_mtx[,ordered_name]
    pval_mtx<- as.matrix(pval_mtx)
    
   
    if (scaled == TRUE){
        message("output scaled data")
        obj = list("exp_mtx" = combined_mtx_scaled, "pval_mtx" = pval_mtx)
    }
    else {
        message("output unscaled data")
        obj = list("exp_mtx" = combined_mtx, "pval_mtx" = pval_mtx)
    }
    
    return(obj)
}

build_cpdb_mtx_wt<- function(cols, scaled = FALSE, pval = 0.01, 
                          sig.cellpairs= 3,
                          sig.timepoint = 2,
                          group_name = "",
                          ph_6mos_pvals, ph_6mos_means,
                          ph_14mos_pvals, ph_14mos_means, 
                          ph_18mos_pvals, ph_18mos_means){
    
    ph_6mos_pvals_sub<- ph_6mos_pvals[,cols]
    interactions_6mos = ph_6mos_pvals[rowSums(ph_6mos_pvals_sub < pval) >= sig.cellpairs, 'interacting_pair']

    ph_14mos_pvals_sub<- ph_14mos_pvals[,cols]
    interactions_14mos = ph_14mos_pvals[rowSums(ph_14mos_pvals_sub < pval) >= sig.cellpairs, 'interacting_pair']

    ph_18mos_pvals_sub<- ph_18mos_pvals[,cols]
    interactions_18mos = ph_18mos_pvals[rowSums(ph_18mos_pvals_sub < pval) >= sig.cellpairs, 'interacting_pair']

    interaction_df = table(c(interactions_6mos, interactions_14mos, interactions_18mos))
    interactions = names(interaction_df)[interaction_df>= sig.timepoint]
    
    message("total number of intersections considered: ", interactions  %>% length())
    
    message("6 mos...")
    ph_6mos_pvals_dedup = ph_6mos_pvals[!ph_6mos_pvals$interacting_pair  %>% duplicated(), ]
    ph_6mos_means_dedup = ph_6mos_means[!ph_6mos_means$interacting_pair  %>% duplicated(), ]
    
    rownames(ph_6mos_pvals_dedup) = ph_6mos_pvals_dedup$interacting_pair
    matchID= match(interactions, rownames(ph_6mos_pvals_dedup))

    ph_6mos_mtx= ph_6mos_means_dedup[matchID, cols]
    rownames(ph_6mos_mtx) = interactions
    ph_6mos_mtx[is.na(ph_6mos_mtx)] = 0
    ph_6mos_mtx  %>% dim()

    ph_6mos_mtxp= ph_6mos_pvals_dedup[matchID, cols]
    rownames(ph_6mos_mtxp) = interactions
    ph_6mos_mtxp[is.na(ph_6mos_mtxp)] = 1
    ph_6mos_mtxp  %>% dim()
    
    message("14 mos...")
    ph_14mos_pvals_dedup = ph_14mos_pvals[!ph_14mos_pvals$interacting_pair  %>% duplicated(), ]
    ph_14mos_means_dedup = ph_14mos_means[!ph_14mos_means$interacting_pair  %>% duplicated(), ]
    
    rownames(ph_14mos_pvals_dedup) = ph_14mos_pvals_dedup$interacting_pair
    matchID= match(interactions, rownames(ph_14mos_pvals_dedup))

    ph_14mos_mtx= ph_14mos_means_dedup[matchID, cols]
    rownames(ph_14mos_mtx) = interactions
    ph_14mos_mtx[is.na(ph_14mos_mtx)] = 0
    ph_14mos_mtx  %>% dim()

    ph_14mos_mtxp= ph_14mos_pvals_dedup[matchID, cols]
    rownames(ph_14mos_mtxp) = interactions
    ph_14mos_mtxp[is.na(ph_14mos_mtxp)] = 1
    ph_14mos_mtxp  %>% dim()
    
    message("18 mos...")
    ph_18mos_pvals_dedup = ph_18mos_pvals[!ph_18mos_pvals$interacting_pair  %>% duplicated(), ]
    ph_18mos_means_dedup = ph_18mos_means[!ph_18mos_means$interacting_pair  %>% duplicated(), ]
    
    rownames(ph_18mos_pvals_dedup) = ph_18mos_pvals_dedup$interacting_pair
    matchID= match(interactions, rownames(ph_18mos_pvals_dedup))

    ph_18mos_mtx= ph_18mos_means_dedup[matchID, cols]
    rownames(ph_18mos_mtx) = interactions
    ph_18mos_mtx[is.na(ph_18mos_mtx)] = 0
    ph_18mos_mtx  %>% dim()

    ph_18mos_mtxp= ph_18mos_pvals_dedup[matchID, cols]
    rownames(ph_18mos_mtxp) = interactions
    ph_18mos_mtxp[is.na(ph_18mos_mtxp)] = 1
    ph_18mos_mtxp  %>% dim()

    
    message("Generating count matrix...")
    combined_mtx = cbind(ph_6mos_mtx, ph_14mos_mtx, ph_18mos_mtx)
    combined_mtx_scaled = t(apply(combined_mtx, 1, scale))
    
    group_n = paste0(c("06mos", "14mos", "18mos"), group_name)
    
    colnames(combined_mtx_scaled) = paste0(rep(colnames(ph_6mos_mtx), times =3), "_", rep(group_n, each=ncol(ph_6mos_mtx)))
    colnames(combined_mtx) = paste0(rep(colnames(ph_6mos_mtx), times =3), "_", rep(group_n, each=ncol(ph_6mos_mtx)))
    
    ordered_name = colnames(combined_mtx_scaled)  %>% sort()
    
    combined_mtx_scaled = combined_mtx_scaled[,ordered_name]
    combined_mtx = combined_mtx[,ordered_name]
    
    message("Generating pvalue matrix...")
    pval_mtx = cbind(ph_6mos_mtxp, ph_14mos_mtxp, ph_18mos_mtxp)
    pval_mtx[pval_mtx == 0] = 1e-5
    pval_mtx = -log10(pval_mtx)
    
    colnames(pval_mtx) = paste0(rep(colnames(ph_6mos_mtxp), times =3), "_", rep(group_n, each=ncol(ph_6mos_mtxp)))
    pval_mtx = pval_mtx[,ordered_name]
    pval_mtx<- as.matrix(pval_mtx)
    
   
    if (scaled == TRUE){
        message("output scaled data")
        obj = list("exp_mtx" = combined_mtx_scaled, "pval_mtx" = pval_mtx)
    }
    else {
        message("output unscaled data")
        obj = list("exp_mtx" = combined_mtx, "pval_mtx" = pval_mtx)
    }
    
    return(obj)
}


build_cpdb_mtx_noSick<- function(cols, scaled = FALSE, pval = 0.01, 
                          sig.cellpairs= 3,
                          sig.timepoint = 2,
                          group_name = "",
                          ph_6mos_pvals, ph_6mos_means,
                          ph_14mos_pvals, ph_14mos_means, 
                          ph_18mos_pvals, ph_18mos_means){
    
    ph_6mos_pvals_sub<- ph_6mos_pvals[,cols]
    interactions_6mos = ph_6mos_pvals[rowSums(ph_6mos_pvals_sub < pval) >= sig.cellpairs, 'interacting_pair']

    ph_14mos_pvals_sub<- ph_14mos_pvals[,cols]
    interactions_14mos = ph_14mos_pvals[rowSums(ph_14mos_pvals_sub < pval) >= sig.cellpairs, 'interacting_pair']

    ph_18mos_pvals_sub<- ph_18mos_pvals[,cols]
    interactions_18mos = ph_18mos_pvals[rowSums(ph_18mos_pvals_sub < pval) >= sig.cellpairs, 'interacting_pair']

    interaction_df = table(c(interactions_6mos, interactions_14mos, interactions_18mos))
    interactions = names(interaction_df)[interaction_df>= sig.timepoint]
    
    message("total number of intersections considered: ", interactions  %>% length())
    
    message("6 mos...")
    ph_6mos_pvals_dedup = ph_6mos_pvals[!ph_6mos_pvals$interacting_pair  %>% duplicated(), ]
    ph_6mos_means_dedup = ph_6mos_means[!ph_6mos_means$interacting_pair  %>% duplicated(), ]
    
    rownames(ph_6mos_pvals_dedup) = ph_6mos_pvals_dedup$interacting_pair
    matchID= match(interactions, rownames(ph_6mos_pvals_dedup))

    ph_6mos_mtx= ph_6mos_means_dedup[matchID, cols]
    rownames(ph_6mos_mtx) = interactions
    ph_6mos_mtx[is.na(ph_6mos_mtx)] = 0
    ph_6mos_mtx  %>% dim()

    ph_6mos_mtxp= ph_6mos_pvals_dedup[matchID, cols]
    rownames(ph_6mos_mtxp) = interactions
    ph_6mos_mtxp[is.na(ph_6mos_mtxp)] = 1
    ph_6mos_mtxp  %>% dim()
    
    message("14 mos...")
    ph_14mos_pvals_dedup = ph_14mos_pvals[!ph_14mos_pvals$interacting_pair  %>% duplicated(), ]
    ph_14mos_means_dedup = ph_14mos_means[!ph_14mos_means$interacting_pair  %>% duplicated(), ]
    
    rownames(ph_14mos_pvals_dedup) = ph_14mos_pvals_dedup$interacting_pair
    matchID= match(interactions, rownames(ph_14mos_pvals_dedup))

    ph_14mos_mtx= ph_14mos_means_dedup[matchID, cols]
    rownames(ph_14mos_mtx) = interactions
    ph_14mos_mtx[is.na(ph_14mos_mtx)] = 0
    ph_14mos_mtx  %>% dim()

    ph_14mos_mtxp= ph_14mos_pvals_dedup[matchID, cols]
    rownames(ph_14mos_mtxp) = interactions
    ph_14mos_mtxp[is.na(ph_14mos_mtxp)] = 1
    ph_14mos_mtxp  %>% dim()
    
    message("18 mos...")
    ph_18mos_pvals_dedup = ph_18mos_pvals[!ph_18mos_pvals$interacting_pair  %>% duplicated(), ]
    ph_18mos_means_dedup = ph_18mos_means[!ph_18mos_means$interacting_pair  %>% duplicated(), ]
    
    rownames(ph_18mos_pvals_dedup) = ph_18mos_pvals_dedup$interacting_pair
    matchID= match(interactions, rownames(ph_18mos_pvals_dedup))

    ph_18mos_mtx= ph_18mos_means_dedup[matchID, cols]
    rownames(ph_18mos_mtx) = interactions
    ph_18mos_mtx[is.na(ph_18mos_mtx)] = 0
    ph_18mos_mtx  %>% dim()

    ph_18mos_mtxp= ph_18mos_pvals_dedup[matchID, cols]
    rownames(ph_18mos_mtxp) = interactions
    ph_18mos_mtxp[is.na(ph_18mos_mtxp)] = 1
    ph_18mos_mtxp  %>% dim()
    
    message("Generating count matrix...")
    combined_mtx = cbind(ph_6mos_mtx, ph_14mos_mtx, ph_18mos_mtx)
    combined_mtx_scaled = t(apply(combined_mtx, 1, scale))
    
    group_n = paste0(c("06mos", "14mos", "18mos"), group_name)
    
    colnames(combined_mtx_scaled) = paste0(rep(colnames(ph_6mos_mtx), times =3), "_", rep(group_n, each=ncol(ph_6mos_mtx)))
    colnames(combined_mtx) = paste0(rep(colnames(ph_6mos_mtx), times =3), "_", rep(group_n, each=ncol(ph_6mos_mtx)))
    
    ordered_name = colnames(combined_mtx_scaled)  %>% sort()
    
    combined_mtx_scaled = combined_mtx_scaled[,ordered_name]
    combined_mtx = combined_mtx[,ordered_name]
    
    message("Generating pvalue matrix...")
    pval_mtx = cbind(ph_6mos_mtxp, ph_14mos_mtxp, ph_18mos_mtxp)
    pval_mtx[pval_mtx == 0] = 1e-5
    pval_mtx = -log10(pval_mtx)
    
    colnames(pval_mtx) = paste0(rep(colnames(ph_6mos_mtxp), times =3), "_", rep(group_n, each=ncol(ph_6mos_mtxp)))
    pval_mtx = pval_mtx[,ordered_name]
    pval_mtx<- as.matrix(pval_mtx)
    
   
    if (scaled == TRUE){
        message("output scaled data")
        obj = list("exp_mtx" = combined_mtx_scaled, "pval_mtx" = pval_mtx)
    }
    else {
        message("output unscaled data")
        obj = list("exp_mtx" = combined_mtx, "pval_mtx" = pval_mtx)
    }
    
    return(obj)
}

merge_mtx<- function(wt_mtx, bcl6_mtx, bcl6_double_mtx){
    
    wt_interactions<- rownames(wt_mtx$exp_mtx)
    bcl6_interactions <- rownames(bcl6_mtx$exp_mtx)
    bcl6_double_interactions<- rownames(bcl6_double_mtx$exp_mtx)
    
    interactions<- c(wt_interactions, bcl6_interactions, bcl6_double_interactions)  %>% sort()  %>% unique()
    
    message("wt...")
    wt_matchID= match(interactions, rownames(wt_mtx$exp_mtx))
    
    wt_mtx$exp_mtx = wt_mtx$exp_mtx[wt_matchID,]
    rownames(wt_mtx$exp_mtx) = interactions
    wt_mtx$exp_mtx[is.na(wt_mtx$exp_mtx)] = 0
    
    wt_mtx$pval_mtx = wt_mtx$pval_mtx[wt_matchID,]
    rownames(wt_mtx$pval_mtx) = interactions
    wt_mtx$pval_mtx[is.na(wt_mtx$pval_mtx)] = 0
    
    message("bcl6...")
    bcl6_matchID= match(interactions, rownames(bcl6_mtx$exp_mtx))
    
    bcl6_mtx$exp_mtx = bcl6_mtx$exp_mtx[bcl6_matchID,]
    rownames(bcl6_mtx$exp_mtx) = interactions
    bcl6_mtx$exp_mtx[is.na(bcl6_mtx$exp_mtx)] = 0
    
    bcl6_mtx$pval_mtx = bcl6_mtx$pval_mtx[bcl6_matchID,]
    rownames(bcl6_mtx$pval_mtx) = interactions
    bcl6_mtx$pval_mtx[is.na(bcl6_mtx$pval_mtx)] = 0
    
    message("bcl6 double...")
    bcl6_double_matchID= match(interactions, rownames(bcl6_double_mtx$exp_mtx))
    
    bcl6_double_mtx$exp_mtx = bcl6_double_mtx$exp_mtx[bcl6_double_matchID,]
    rownames(bcl6_double_mtx$exp_mtx) = interactions
    bcl6_double_mtx$exp_mtx[is.na(bcl6_double_mtx$exp_mtx)] = 0
    
    bcl6_double_mtx$pval_mtx = bcl6_double_mtx$pval_mtx[bcl6_double_matchID,]
    rownames(bcl6_double_mtx$pval_mtx) = interactions
    print(range(bcl6_double_mtx$pval_mtx))
    bcl6_double_mtx$pval_mtx[is.na(bcl6_double_mtx$pval_mtx)] = 0
    
    message(nrow(bcl6_mtx$exp_mtx) == nrow(bcl6_double_mtx$exp_mtx))
    
    
    exp_mtx = cbind(wt_mtx$exp_mtx, bcl6_mtx$exp_mtx, bcl6_double_mtx$exp_mtx)
    message("scaling the count matrix")
    exp_mtx_scaled = t(apply(exp_mtx, 1, scale))
    colnames(exp_mtx_scaled) = colnames(exp_mtx)

    pval_mtx = cbind(wt_mtx$pval_mtx, bcl6_mtx$pval_mtx, bcl6_double_mtx$pval_mtx)
    
    
    clusterID = gsub(pattern = "_06mos.*|_14mos.*|_18mos.*|_sick.*", replacement = "",x = colnames(exp_mtx_scaled))
    genotype = paste0("bcl", gsub(pattern = ".*bcl|.*wt", replacement = "", x = colnames(exp_mtx_scaled)))
    genotype[genotype=='bcl'] = 'wt'
    age = gsub(pattern = "_bcl.*|_wt", replacement = "",x = colnames(exp_mtx_scaled))
    
    ordered_colname = data.frame(colID = colnames(exp_mtx_scaled),  clusterID = clusterID, genotype = genotype, age = age) %>% 
        arrange(clusterID, age)  %>% pull(colID)

    exp_mtx_scaled <- exp_mtx_scaled[,ordered_colname]
    pval_mtx <- pval_mtx[,ordered_colname]
    
    combined_mtx = list("exp_mtx" = exp_mtx_scaled, "pval_mtx" = pval_mtx)
    
    return(combined_mtx)

}

plt_cpdb_dotplot<- function(exp_mtx, pval_mtx, dt.size = 1.5, 
                            column_title = "merged"){
    
    col_len = ncol(exp_mtx)/9
    
    message("quantile: ", c(quantile(exp_mtx, 0.1), quantile(exp_mtx, 0.5), quantile(exp_mtx, 0.95)))

    col_fun = circlize::colorRamp2(c(-1.5, 0, 1.5), c("#440154FF", "#238A8DFF", "#FDE725FF"))

    fact = gsub("06mos.*|14mos.*|18mos.*|sick.*", "", colnames(exp_mtx))
    dend = cluster_between_groups(exp_mtx, fact)
    
    legend_dot_size = unit(dt.size,"mm")
    
    lgd = Legend(labels = c(0, 1, 2, 3, 4, 5),
                 title = "-logPvalue",
                 type = "points",
                 graphics = list(
                    function(x, y, w, h) grid.circle(x, y, r= sqrt(0/5) * legend_dot_size,
                                                    gp = gpar(fill = "black")),
                    function(x, y, w, h) grid.circle(x, y, r= sqrt(1/5) * legend_dot_size,
                                                    gp = gpar(fill = "black")),
                    function(x, y, w, h) grid.circle(x, y, r= sqrt(2/5) * legend_dot_size,
                                                    gp = gpar(fill = "black")),
                    function(x, y, w, h) grid.circle(x, y, r= sqrt(3/5) * legend_dot_size,
                                                    gp = gpar(fill = "black")),
                    function(x, y, w, h) grid.circle(x, y, r= sqrt(4/5) * legend_dot_size,
                                                    gp = gpar(fill = "black")),
                    function(x, y, w, h) grid.circle(x, y, r= sqrt(5/5) * legend_dot_size,
                                                    gp = gpar(fill = "black"))            
                            )
                     )

    layer_fun_dot = function(j, i, x, y, w, h, fill){
        grid.rect(x = x, y = y, width = w, height = h, 
                  gp = gpar(col = NA, fill = NA))

        grid.circle(x=x, y=y, r=sqrt(pindex(pval_mtx, i, j)/5) * legend_dot_size,
                    gp = gpar(fill = col_fun(pindex(exp_mtx, i, j)), col = NA))
    }
    
    clusterID =  gsub("_06mos.*|_14mos.*|_18mos.*|_sick.*", "", colnames(exp_mtx))
    genotype = paste0("bcl", gsub(pattern = ".*bcl|.*wt", replacement = "", x = colnames(exp_mtx)))
    genotype[genotype=='bcl'] = 'wt'
    age = rep(rep(c("06mos", "14mos", "18mos"), times = c(3,3,3)), ncol(exp_mtx)/9)
    
    print(table(genotype))
    print(table(age))
    anno_df = data.frame(clusterID = clusterID, genotype = genotype, age = age)
    
    # clusterID_cols
    clusterID_cols= kelly(n = length(unique(clusterID)) +2)[-c(1,2)] %>% unname()
    clusterID_cols = setNames(clusterID_cols, unique(clusterID))
    # age_cols
    age_cols = brewer.pal(4, "Dark2")[c(1:3)]
    age_cols = setNames(age_cols, c("06mos", "14mos", "18mos"))
    
    genotype_cols = brewer.pal(5, "Set1")[c(1,3,5)]
    genotype_cols = setNames(genotype_cols, c("wt", "bcl6", "bcl6_double"))
    
    col_list = list(clusterID = clusterID_cols, genotype = genotype_cols, age = age_cols)

    column_ha <- HeatmapAnnotation(
        clusterID = clusterID, age = age, genotype = genotype, col = col_list,
        annotation_legend_param = list(clusterID = list(title = "annotation"),
                                       age = list(title = "age"),
                                       genotype = list(title = "genotype")), show_annotation_name = TRUE)
    
    dotplot<- Heatmap(exp_mtx, 
        name = "meanExpr",
        column_title = column_title,
        column_names_side = "top",
        column_title_side = "bottom",
        show_column_names = FALSE,
        column_names_gp = gpar(fontsize = 10),
        column_title_gp = gpar(fontsize = 12),
        cluster_columns = dend,
        column_split = col_len,
        rect_gp = gpar(type = "none"),
        layer_fun = layer_fun_dot,
        col = col_fun,
        row_km = 1,
        row_names_gp = gpar(fontsize = 8),
        top_annotation = column_ha,
        border = "black",
        show_row_dend =FALSE,
        show_column_dend = FALSE)
        
    return(dotplot)
    
}


plt_cpdb_dotplot_with_sick<- function(exp_mtx, pval_mtx, dt.size = 1.5, 
                            column_title = "merged"){
    
    col_len = ncol(exp_mtx)/11
    
    message("quantile: ", c(quantile(exp_mtx, 0.1), quantile(exp_mtx, 0.5), quantile(exp_mtx, 0.95)))

    col_fun = circlize::colorRamp2(c(-1.5, 0, 1.5), c("#440154FF", "#238A8DFF", "#FDE725FF"))

    fact = gsub("06mos.*|14mos.*|18mos.*|sick.*", "", colnames(exp_mtx))
    dend = cluster_between_groups(exp_mtx, fact)
    
    legend_dot_size = unit(dt.size,"mm")
    
    lgd = Legend(labels = c(0, 1, 2, 3, 4, 5),
                 title = "-logPvalue",
                 type = "points",
                 graphics = list(
                    function(x, y, w, h) grid.circle(x, y, r= sqrt(0/5) * legend_dot_size,
                                                    gp = gpar(fill = "black")),
                    function(x, y, w, h) grid.circle(x, y, r= sqrt(1/5) * legend_dot_size,
                                                    gp = gpar(fill = "black")),
                    function(x, y, w, h) grid.circle(x, y, r= sqrt(2/5) * legend_dot_size,
                                                    gp = gpar(fill = "black")),
                    function(x, y, w, h) grid.circle(x, y, r= sqrt(3/5) * legend_dot_size,
                                                    gp = gpar(fill = "black")),
                    function(x, y, w, h) grid.circle(x, y, r= sqrt(4/5) * legend_dot_size,
                                                    gp = gpar(fill = "black")),
                    function(x, y, w, h) grid.circle(x, y, r= sqrt(5/5) * legend_dot_size,
                                                    gp = gpar(fill = "black"))            
                            )
                     )

    layer_fun_dot = function(j, i, x, y, w, h, fill){
        grid.rect(x = x, y = y, width = w, height = h, 
                  gp = gpar(col = NA, fill = NA))

        grid.circle(x=x, y=y, r=sqrt(pindex(pval_mtx, i, j)/5) * legend_dot_size,
                    gp = gpar(fill = col_fun(pindex(exp_mtx, i, j)), col = NA))
    }
    
    clusterID =  gsub("_06mos.*|_14mos.*|_18mos.*|_sick.*", "", colnames(exp_mtx))
    genotype = paste0("bcl", gsub(pattern = ".*bcl|.*wt", replacement = "", x = colnames(exp_mtx)))
    genotype[genotype=='bcl'] = 'wt'
    age = rep(rep(c("06mos", "14mos", "18mos", "sick"), times = c(3,3,3,2)), ncol(exp_mtx)/11)
    
    message(table(genotype))
    message(table(age))
    anno_df = data.frame(clusterID = clusterID, genotype = genotype, age = age)
    
    # clusterID_cols
    clusterID_cols= kelly(n = length(unique(clusterID)) +2)[-c(1,2)] %>% unname()
    clusterID_cols = setNames(clusterID_cols, unique(clusterID))
    # age_cols
    age_cols = brewer.pal(4, "Dark2")
    age_cols = setNames(age_cols, c("06mos", "14mos", "18mos", "sick"))
    
    genotype_cols = brewer.pal(5, "Set1")[c(1,3,5)]
    genotype_cols = setNames(genotype_cols, c("wt", "bcl6", "bcl6_double"))
    
    col_list = list(clusterID = clusterID_cols, genotype = genotype_cols, age = age_cols)

    column_ha <- HeatmapAnnotation(
        clusterID = clusterID, age = age, genotype = genotype, col = col_list,
        annotation_legend_param = list(clusterID = list(title = "annotation"),
                                       age = list(title = "age"),
                                       genotype = list(title = "genotype")), show_annotation_name = TRUE)
    
    dotplot<- Heatmap(exp_mtx, 
        name = "meanExpr",
        column_title = column_title,
        column_names_side = "top",
        column_title_side = "bottom",
        show_column_names = FALSE,
        column_names_gp = gpar(fontsize = 10),
        column_title_gp = gpar(fontsize = 12),
        cluster_columns = dend,
        column_split = col_len,
        rect_gp = gpar(type = "none"),
        layer_fun = layer_fun_dot,
        col = col_fun,
        row_km = 6,
        row_names_gp = gpar(fontsize = 8),
        top_annotation = column_ha,
        border = "black",
        show_row_dend =FALSE,
        show_column_dend = FALSE)
        
    return(dotplot)
    
}

legend_dot_size = unit(1.5,"mm")
lgd = Legend(labels = c(0, 1, 2, 3, 4, 5),
         title = "-logPvalue",
         type = "points",
         graphics = list(
            function(x, y, w, h) grid.circle(x, y, r= sqrt(0/5) * legend_dot_size,
                                            gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x, y, r= sqrt(1/5) * legend_dot_size,
                                            gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x, y, r= sqrt(2/5) * legend_dot_size,
                                            gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x, y, r= sqrt(3/5) * legend_dot_size,
                                            gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x, y, r= sqrt(4/5) * legend_dot_size,
                                            gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x, y, r= sqrt(5/5) * legend_dot_size,
                                            gp = gpar(fill = "black"))            
                    )
             )


EnrichedView = function(enrichment,
                        rank_by = "pvalue",
                        mode = 1,
                        subset = NULL,
                        top = 0, bottom = 0,
                        x = "LogFDR",
                        charLength = 40,
                        filename = NULL,
                        width = 7, height = 4, ...){

  if(is(enrichment, "enrichResult")) enrichment = slot(enrichment, "result")
  if(is(enrichment, "gseaResult")) enrichment = slot(enrichment, "result")
  # No enriched pathways
  if(is.null(enrichment) || nrow(enrichment)==0){
    p1 = noEnrichPlot("No enriched terms")
    if(!is.null(filename)){
      ggsave(plot=p1,filename=filename, units = "in", width=width, height=height, ...)
    }
    return(p1)
  }
  flag = ifelse("NES"%in%colnames(enrichment), TRUE, FALSE)
  if(!flag) colnames(enrichment)[colnames(enrichment)=="Count"] = "NES"
  ## Rank enriched pathways ##
  enrichment$logP = round(-log10(enrichment$pvalue), 1)
  enrichment$logFDR = round(-log10(enrichment$p.adjust), 1)
  enrichment = enrichment[!is.na(enrichment$ID), ]
  if(tolower(rank_by) == "pvalue"){
    enrichment = enrichment[order(enrichment$pvalue, -abs(enrichment$NES)), ]
  }else if(tolower(rank_by) == "nes"){
    enrichment = enrichment[order(-abs(enrichment$NES), enrichment$pvalue), ]
  }

  ## Normalize term description ##
  terms = as.character(enrichment$Description)
  terms = lapply(terms, function(x,k){
    x = as.character(x)
    if(nchar(x)>k){x=substr(x,start=1,stop=k)}
    return(x)}, charLength)
  enrichment$Description = do.call(rbind, terms)
  enrichment = enrichment[!duplicated(enrichment$Description),]

  # Select pathways to show ##
  if(bottom>0){
    tmp = enrichment[enrichment$NES<0, ]
    subset = c(subset, tmp$ID[1:min(nrow(tmp), bottom)])
  }
  if(top>0){
    tmp = enrichment[enrichment$NES>0, ]
    subset = c(subset, tmp$ID[1:min(nrow(tmp), top)])
  }
  if(!is.null(subset)){
    idx = enrichment$ID %in% subset
    if(sum(idx)==0) return(noEnrichPlot("No eligible terms!!!"))
    enrichment = enrichment[idx, ]
  }
  # enrichment$ID = factor(enrichment$ID, levels=c(pid_neg, pid_pos))
  # enrichment = enrichment[order(enrichment$ID), ]
  enrichment$Description = factor(enrichment$Description,
                                  levels=unique(enrichment$Description))
  enrichment$ID = factor(enrichment$ID, levels=unique(enrichment$ID))

  ## Prepare data for plotting ##
  if(x=="NES"){
    enrichment$x = enrichment$NES
    enrichment$size = enrichment$logFDR
  }else if(x=="LogP"){
    enrichment$x = enrichment$logP
    enrichment$size = abs(enrichment$NES)
  }else{
    enrichment$x = enrichment$logFDR
    enrichment$size = abs(enrichment$NES)
  }

  # Visualize the top enriched terms
  enrichment$col = "UG"
  enrichment$col[enrichment$NES<0] = "DG"

  if(mode==1){
    ## Plot the figure ##
    p1 = ggplot(data=enrichment, aes_string(x="x", y="Description", size="size"))
    p1 = p1 + geom_point(aes(color = col))
    p1 = p1 + scale_color_manual(values = c("DG"="#377eb8", "UG"="#e41a1c"))
    p1 = p1 + theme(panel.grid.major=element_line(colour="gray90"),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank())
    p1 = p1 + theme(legend.position="right")
    p1 = p1 + theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))
    p1 = p1 + theme_bw(base_size = 14)
    p1 = p1 + theme(plot.title = element_text(hjust = 0.5))
    p1 = p1 + guides(color = "none")
  }else if(mode == 2){
    idx = (max(enrichment$x)-enrichment$x) > enrichment$x
    enrichment$hjust = 1.1
    enrichment$hjust[idx] = -0.1
    p1 = ggplot(enrichment, aes_string("x", "ID", label = "Description"))
    p1 = p1 + geom_point(aes_string(color = "col", size = "size"))
    p1 = p1 + xlim(0, NA)
    p1 = p1 + geom_text(aes_string(hjust = "hjust"))
    p1 = p1 + theme_bw(base_size = 14) + theme(plot.title = element_text(hjust = 0.5))
  }
  if(x=="NES"){
    p1 = p1 + labs(x = "Enrichment score", y = NULL, color = NULL,
                   size = expression(-log[10]*FDR))
    if(!flag) p1 = p1 + labs(x = "Count")
  }else{
    p1 = p1 + labs(x = expression(-log[10]*FDR), y = NULL, color = NULL, size = "NES")
    if(!flag) p1 = p1 + labs(size = "Count")
  }

  if(!is.null(filename)){
    ggsave(plot=p1, filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p1)
}

noEnrichPlot = function(main = "No enriched terms"){
  p1 = ggplot()
  p1 = p1 + geom_text(aes(x=0, y=0, label="No enriched terms"), size=6)
  p1 = p1 + labs(title=main)
  p1 = p1 + theme(plot.title = element_text(size=12))
  p1 = p1 + theme_bw(base_size = 14)
  p1 = p1 + theme(plot.title = element_text(hjust = 0.5))
  p1
}

# run_de<- function(obj, outdir, prefix){
    
    
#     n_median = obj@meta.data  %>% count(pool_id)  %>% pull(n)  %>% median()
    
#     cells_to_keep <- obj@meta.data %>% 
#                         group_by(pool_id) %>% 
#                         slice_sample(n = min(n(), n_median)) %>%
#                         pull(cell = rownames(.))

    
#     obj <- subset(obj, cells = cells_to_keep)
    
#     df <- obj@meta.data %>% count(pool_id)

#     ggplot(df, aes(x = reorder(pool_id, n), y = n)) +
#         geom_col() + coord_flip() + 
#         labs(x = "sample", y = "Number of cells") + theme_bw()
    
    
#     compare <- obj$genotype  %>% unique( )
#     print(compare)
#     markers <- presto::wilcoxauc(obj, "genotype" , assay = "data", groups_use = compare)

#     write_tsv(markers, paste0(outdir, prefix, "_all_markers.tsv"))

#     hist(markers$padj, breaks = 20)
#     hist(markers$logFC, breaks = 20)

#     top_markers<- markers  %>% filter(padj <=0.05, abs(logFC) >=0.1)  %>% filter() %>% arrange(padj, -abs(logFC), .by_group = TRUE)

#     write_tsv(top_markers,
#             paste0(outdir, prefix ,"_significant_markers_padj_0.05_abslogFC_0.1.tsv"))


#     up_in_double<- top_markers %>% 
#         filter(group == "CD70-/-;Bcl6tg/+", logFC > 0) %>% 
#         arrange(padj)  %>% 
#         slice_head(n = 10)

#     down_in_double<- top_markers %>% 
#         filter(group == "CD70-/-;Bcl6tg/+", logFC < 0) %>% 
#         arrange(padj)  %>% 
#         slice_head(n = 10)


#     write_tsv(up_in_double,
#             paste0(outdir, prefix ,"_top10_up_regulated_in_double.tsv"))

#     write_tsv(down_in_double,
#             paste0(outdir, prefix ,"_top10_down_regulated_in_double.tsv"))

#     genelist = c(up_in_double$feature, down_in_double$feature)  %>% unique()

#     p<- DotPlot(object = obj, features = genelist, group.by = "genotype", scale = T) +
#             theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = prefix, y = "genotype")
#     p 
#     ggsave(paste0(outdir, prefix ,"_top10_up_and_down_regulated_genes_in_double.pdf"), p, width = 7, height = 4)
    
#     p<- VlnPlot(obj, features = genelist, group.by = "genotype", ncol = 10, pt.size = 0)

#     ggsave(paste0(outdir, prefix ,"_top10_up_and_down_regulated_genes_in_double_vlnplot.pdf"), p, width = 10, height = 5)
# }

# run_de_pseudobulk <- function(obj, outdir, prefix) {

#     set.seed(123)

#     raw_counts <- AggregateExpression(obj, group.by = "orig.ident", slot = "counts")
#     cell_counts <- table(obj$orig.ident)
#     mean_counts <- sweep(raw_counts$RNA, 2, cell_counts[colnames(raw_counts$RNA)], "/")

#     keep_genes <- rowSums(raw_counts$RNA >= 3) >= 3
#     mean_counts <- mean_counts[keep_genes, ]

#     log_mean_counts <- log1p(mean_counts)

#     genotype_group <- obj@meta.data %>% distinct(orig.ident, genotype) %>% deframe()
#     genotype_group <- genotype_group[colnames(mean_counts)]

#     results <- sapply(rownames(mean_counts), function(gene) {
#         g1_log <- log_mean_counts[gene, genotype_group == "Bcl6tg/+"]
#         g2_log <- log_mean_counts[gene, genotype_group == "CD70-/-;Bcl6tg/+"]

#         g1_raw <- mean_counts[gene, genotype_group == "Bcl6tg/+"]
#         g2_raw <- mean_counts[gene, genotype_group == "CD70-/-;Bcl6tg/+"]

#         if (length(unique(c(g1_log, g2_log))) <= 1) {
#             return(c(p = 1, logFC = 0))
#         }

#         wt <- wilcox.test(g2_log, g1_log, exact = FALSE)
#         c(p = wt$p.value, logFC = log2(mean(g2_raw) + 1) - log2(mean(g1_raw) + 1))
#     })

#     markers <- as.data.frame(t(results))
#     markers$feature <- rownames(markers)
#     markers$padj <- p.adjust(markers$p, method = "BH")
#     markers <- markers %>% arrange(padj)

#     write_tsv(markers, paste0(outdir, prefix, "_all_markers.tsv"))

#     # Volcano plot
#     markers$neg_log10_padj <- -log10(markers$padj)
#     markers$sig <- case_when(
#         markers$padj <= 0.1 & markers$logFC >= 0.5  ~ "Up",
#         markers$padj <= 0.1 & markers$logFC <= -0.5 ~ "Down",
#         TRUE ~ "NS"
#     )

#     top_genes <- markers %>%
#         filter(sig != "NS") %>%
#         arrange(padj) %>%
#         slice_head(n = 20)

#     p_volcano <- ggplot(markers, aes(x = logFC, y = neg_log10_padj, color = sig)) +
#         geom_point(alpha = 0.5, size = 1) +
#         scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "grey70")) +
#         geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
#         geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "black") +
#         ggrepel::geom_text_repel(
#             data = top_genes,
#             aes(label = feature),
#             size = 3, max.overlaps = 20, color = "black"
#         ) +
#         labs(
#             title = paste0(prefix, ":\nVolcano plot padj<= 0.1 abs(logFC) >= 0.5"),
#             x = "log2 Fold Change",
#             y = "-log10(adjusted p-value)",
#             color = NULL
#         ) +
#         theme_classic(base_size = 14)

#     ggsave(paste0(outdir, prefix, "_volcano_plot.pdf"), p_volcano, width = 7, height = 5)

#     # Histograms
#     p_padj <- ggplot(markers, aes(x = padj)) +
#         geom_histogram(bins = 30, fill = "#4DBBD5", color = "black") +
#         geom_vline(xintercept = 0.1, linetype = "dashed", color = "red") +
#         labs(
#             title = paste0(prefix, ":\nAdjusted p-value distribution"),
#             x = "Adjusted p-value",
#             y = "Number of genes"
#         ) + theme_classic(base_size = 14)

#     p_logfc <- ggplot(markers, aes(x = logFC)) +
#         geom_histogram(bins = 30, fill = "#E64B35", color = "black") +
#         geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +
#         labs(
#             title = paste0(prefix, ":\nlog2 Fold Change distribution"),
#             x = "log2 Fold Change",
#             y = "Number of genes"
#         ) + theme_classic(base_size = 14)

#     ggsave(paste0(outdir, prefix, "_padj_histogram.pdf"), p_padj, width = 5, height = 4)
#     ggsave(paste0(outdir, prefix, "_logFC_histogram.pdf"), p_logfc, width = 5, height = 4)

#     # Significant markers
#     markers_sig <- markers %>% filter(padj <= 0.1, abs(logFC) >= 0.5) %>% arrange(padj)
#     write_tsv(markers_sig, paste0(outdir, prefix, "_significant_markers_padj_0.1_abslogFC_0.5.tsv"))

#     up_in_double <- markers_sig %>%
#         filter(logFC > 0) %>%
#         arrange(padj) %>%
#         slice_head(n = 10)

#     down_in_double <- markers_sig %>%
#         filter(logFC < 0) %>%
#         arrange(padj) %>%
#         slice_head(n = 10)

#     genelist <- c(up_in_double$feature, down_in_double$feature) %>% unique()

#     if (length(genelist) == 0) {
#         message("No significant DE genes for ", prefix)
#     } else {
#         write_tsv(up_in_double, paste0(outdir, prefix, "_top10_up_regulated_in_double.tsv"))
#         write_tsv(down_in_double, paste0(outdir, prefix, "_top10_down_regulated_in_double.tsv"))

#         p <- DotPlot(object = obj, features = genelist, group.by = "genotype", scale = TRUE) +
#             theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#             labs(title = prefix, y = "genotype")

#         ggsave(paste0(outdir, prefix, "_top10_up_and_down_regulated_genes_in_double.pdf"), p, width = 7, height = 4)

#         p <- VlnPlot(obj, features = genelist, group.by = "genotype", ncol = 10, pt.size = 0)
#         ggsave(paste0(outdir, prefix, "_top10_up_and_down_regulated_genes_in_double_vlnplot.pdf"), p, width = 10, height = 5)

#         p <- VlnPlot(obj, features = genelist, group.by = "orig.ident", split.by = "genotype",
#                       pt.size = 0, stack = TRUE, sort = TRUE, flip = TRUE)
#         ggsave(paste0(outdir, prefix, "_top10_up_and_down_regulated_genes_in_double_vlnplot_by_sample.pdf"), p, width = 14, height = 8)
#     }
# }