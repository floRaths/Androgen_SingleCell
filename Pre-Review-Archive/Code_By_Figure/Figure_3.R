


panel_3b <- function(){

query <- c("CIDEA", "LPL", "PPARG")

#P3 <- 
Sobj_integrated@assays$RNA@data[query, ] %>% as.matrix() %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("ID") %>% 
  left_join(select(rownames_to_column(Sobj_integrated@meta.data, "ID"), ID, CellType, Type, Subcluster)) %>% 
  #filter(CellType == "Fibroblast") %>%
  #filter(Subcluster %in% levels(Sobj_integrated)[1:3]) %>% 
  pivot_longer(query, names_to = "Gene") %>% filter(Gene != "AR") %>% 
  mutate(value = value + 1) %>% 
  
  ggplot(aes(x = Type, y = value, fill = Type)) + 
  geom_violin(size = 1, color = "grey20") + 
  
  scale_fill_manual(values = c("#9d3396", "#eb933e")) + 
  
  scale_y_log10() + 
  #scale_y_sqrt() + 
  facet_wrap(facets = "Gene") + 
  
  ylab("expression") + 
  
  theme(text = element_text(family = "Lato", size = 20), 
        legend.position = "top", 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(color = "grey30", size = 2, fill = NA),
        panel.grid = element_line(color = "grey90"))




P2 <- 
adipocyte_features_2020_10_06_v1 %>% #filter(case_label %notin% c("CF-787-821", "TM-7714")) %>% 
  group_by(case_label, type) %>% 
  summarise(mean = mean(area)) %>% 
  ggplot(aes(x = type, y = mean, fill = type)) + 
  
  geom_boxplot(#aes(color = type), 
    color = "grey20",
    width = 0.5, alpha = 1, 
    lwd = 1,
    position=position_dodge(0.7), outlier.alpha = 0.8) +
  
  
  geom_point(color = "grey60", size = 4, alpha = 0.65, 
             position = position_jitterdodge(jitter.width = 0.2,
                                             jitter.height = 0,
                                             dodge.width = 0.9)) +
  
  
  scale_fill_manual(values = c("#9d3396", "#eb933e")) +
  scale_color_manual(values = c(cols, "grey")) +
  
  stat_compare_means(paired = F, method = "t.test", 
                     label.x.npc = 0.25) +
  
  ggtitle("avg. adipocyte area") +
  
  #coord_flip() + 
  
  theme(text = element_text(family = "Lato", size = 20), 
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(color = "grey30", size = 2, fill = NA),
        panel.grid = element_line(color = "grey90"))




P1 + P2 + patchwork::plot_layout(ncol = 1, heights = c(2.5,1))

}
p %>% save_x(data = ., name = "Lipolysis_Markers_And_Area", 1.5, 5, 5, svg = T)



panel_S3a <- function(){

celltype = "Adipocyte"

Sobj_integrated <- AddMetaData(Sobj_integrated, metadata = as.data.frame(Sobj_integrated@reductions$umap@cell.embeddings))
df_rna = Sobj_integrated@meta.data

set.seed(467)
rows_rna  <- sample(nrow(df_rna)) #, 10000)


P2 <- 
  df_rna[rows_rna, ] %>% #filter(CellType == "Myeloid") %>% 
  
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = factor(Sample, levels = levels))) + 
  
  geom_point(data = filter(select(df_rna[rows_rna, ], -Type)), color = "grey88",
             alpha = 1, 
             size = 6, 
             stroke = 1) +
  
  geom_point(alpha = 1, 
             size = 3, 
             stroke = 0) +
  
  facet_wrap(~Type, ncol = 1) +
  
  scale_color_manual(values = cols) + 
  
  look + 
  theme(strip.background = element_blank(), panel.background = element_blank(),
        strip.text = element_blank())




}
P2 <- panel_S3a()
P2 %>% save_x(data = ., name = paste0(celltype, "Subcluster_Umap"), 0.8, 6, 10, svg = F)






panel_S3c <- function(){

### find TFs in Adipocytes with AR motif
target_tfs <- NarrowPeaks_DistanceToSummit200_AllCellTypes_Nov16 %>% 
  filter(CellType == "Adipocyte", Motif == "AR_689", GeneRna %in% TFs$Gene) %>% 
  pull(GeneRna) %>% unique()

df <- readRDS("Output/Adipocyte_DE_Genes_NoFC_Cutoff.rds") %>% filter(Cluster== "Adipocyte", Gene %in% TFs$Gene)

df %>% 
  ggplot(aes(x = avg_logFC, y = -log10(p_val_adj))) + 
  
  geom_vline(xintercept =  0, color = "grey70", size = 1) +
  geom_hline(yintercept =  -log10(0.05), color = "grey70", size = 1, linetype = "dashed") +
  
  geom_point(data = filter(df, Gene %notin% target_tfs), size = 4, alpha = 0.5, color = "grey65") + 
  geom_point(data = filter(df, Gene %in% target_tfs),    size = 5, alpha = 1, shape = 17, color = "#9d3396") + 
  geom_text_repel(data = filter(df, Gene %in% target_tfs, abs(avg_logFC) > 0.3), aes(label = Gene), size = 5, alpha = 1) + 
  
  xlim(c(-1, 1)) +
  
  theme(text = element_text(family = "Lato", size = 25),
        #axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        #axis.ticks.y = element_blank(), 
        panel.grid = element_blank(),
        legend.position = "bottom",
        panel.background = element_blank(),
        panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
        )

}

p <- panel_S3c()

p %>% save_x(data = ., name = paste0("Adipocyte_AR_Tfs_Vulcano"), 1, 10, 10, svg = T)



panel_S3d <- function() {
  
  avg_fib <- readRDS("Output/Previous_Analysis/Fibroblast_Subcluster_avgExpression.rds")
  marks_fib  <- readRDS("Output/Previous_Analysis/Fibroblast_Subcluster_Markers.rds")
  
  genes1 <- 
    marks_fib %>% 
    filter(!str_detect(gene, "-AS1|-AS2|\\.|LINC0")) %>% 
    group_by(cluster) %>% 
    top_n(25, pct.1) %>% 
    top_n(5, avg_log2FC) %>% 
    arrange(cluster) %>% pull(gene)
  
  
  
  
  mat1 <- 
    avg_fib %>% 
    rownames_to_column("gene") %>% 
    filter(gene %in% c(genes1)) %>% 
    column_to_rownames("gene")
  
  
  
  pheatmap(t(mat1[genes1,levels(marks_fib$cluster)]), 
           #clustering_distance_rows = "correlation",
           #clustering_distance_cols = "correlation", 
           #clustering_method = "ward.D2",
           cluster_rows = F,
           cluster_cols = F,
           scale = "column", 
           color = colorRampPalette(brewer.pal(11, "BrBG"))(100)
           #color = viridis::magma(n = 100, end = 0.9)
  )
  
}
p <- panel_S3d()
p %>% save_x(data = ., name = paste0("Fibroblast_Marker_Calling"), 1, 12, 10, svg = T)

#panel_S3d <- function(){

Prop <- 
  as.data.frame(prop.table(table(Sobj_integrated$Subcluster, 
                                 Sobj_integrated$Sample), 
                           margin = 2)) %>% 
  dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq") %>% 
  mutate(Type = substring(.$Sample, 1, 2))# %>% arrange(factor(Cluster, levels = levels), desc(Freq))



Prop %>% #filter(!is.na(Freq), Cluster %in% c("T-Effector", "CD4_T", "CD8_T", "NK")) %>% 
  ggplot(aes(x = Cluster, y = Freq, fill = Type)) + 
  
  #geom_violin(alpha = 0.65, trim = F, scale = "width" , draw_quantiles = T) +
  
  geom_boxplot(color = "grey40", width = 0.5, alpha = 1, lwd = 1.25,
               position=position_dodge(0.75), outlier.alpha = 0) +
  
  geom_point(aes(color = Sample), 
             position = position_jitterdodge(jitter.width = 0.5, 
                                             jitter.height = 0,
                                             dodge.width = 0.4), 
             alpha = 1,
             size  = 4) +
  scale_y_sqrt() + 
  
  scale_color_manual(values =  cols) +
  
  ggtitle("Fibroblast Proportions") + 
  scale_fill_manual(values = c("#9d3396", "#eb933e")) +
  theme(text = element_text(family = "Lato", size = 25),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        #axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.background = element_blank(), 
        panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
  ) 

}
#p <- panel_S3d()
#p %>% save_x(data = ., name = paste0("Fibroblast_Subcluster_Proportions"), 1, 12, 10, svg = T)





panel_S3f <- function() {
  data %>% select(Sample, CellType, Subcluster, Type, query) %>% 
    #filter(CellType %in% celltype) %>% 
    #filter(Subcluster %in% c("CD4_T", "CD8_T")) %>% 
    #scores %>% select(Subcluster, query) %>% 
    pivot_longer(query, names_to = "Module") %>% 
    #group_by(Sample, Subcluster, Module) %>% 
    #mutate(sample_avg = mean(value)) %>% 
    
    ggplot(aes(x = factor(Subcluster, levels = levels), y = value, fill = Type)) + 
    
    geom_split_violin(scale = "width", 
                      alpha = 0.75, 
                      trim = F,
                      lwd = 1.25, 
                      color = "grey30") + 
    
    
    geom_boxplot(aes(color = Type), 
                 fill = "grey30", 
                 width = 0.1, alpha = 1, 
                 lwd = 1.25,
                 position=position_dodge(0.7), outlier.alpha = 0) +
    
    scale_fill_viridis_d(option = "C", begin = 0.4, end = 0.8) +
    scale_color_manual(values = c("grey30", "grey30")) +
    
    stat_summary(fun = "median", 
                 geom = "point", size = 1.5,
                 position = position_dodge(0.7),
                 color = "floralwhite") +
    
    facet_wrap(vars(Module), scales = "free_y") +
    
    #scale_y_sqrt() +
    
    theme(text = element_text(family = "Lato", size = 25),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          #axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none",
          panel.background = element_blank(), 
          panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
    )
  
}


query= c("BIOCARTA_ECM_PATHWAY", "REACTOME_LAMININ_INTERACTIONS", "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX", "KEGG_ECM_RECEPTOR_INTERACTION")

p <- panel_S3f()

p %>% save_x(data = ., name = paste0("Fibroblast_ECM_Pathways"), 1.5, 12, 10, svg = T)





Lipolysis_Markers %>% 
  left_join(filter(Treatment_Response, Cluster == "Adipocyte")) %>% 
  ggplot(aes(x = avg_logFC, y = -log10(p_val_adj), color = `lipo function`)) + 
  geom_vline(xintercept = 0, color = "grey70") +
  
  geom_point(size = 3) + 
  geom_text_repel(aes(label = Gene), size = 4)  +
  ggtitle("lipolysis related genes in adipocytes") +
  
  scale_color_d3() +
  theme(text = element_text(family = "Lato", size = 25),
                                          #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                                          #axis.text.y = element_blank(),
                                          #axis.title.x = element_blank(),
                                          #axis.title.y = element_blank(),
                                          panel.grid = element_blank(),
                                          legend.position = "bottom",
                                          panel.background = element_blank(), 
                                          panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
  )






p <- 
Lipolysis_Markers %>% rename("LipoFunction" = 2) %>% 
  left_join(filter(Treatment_Response, Cluster == "Adipocyte")) %>% 
  ggplot(aes(x = avg_logFC, y = -log10(p_val_adj), color = LipoFunction)) + 
  
  geom_vline(xintercept = 0, color = "grey70") +
  
  geom_point(size = 3) + 
  geom_text_repel(aes(label = Gene), size = 4)  +
  ggtitle("lipolysis related genes in adipocytes") +
  
  scale_color_d3() +
  theme(text = element_text(family = "Lato", size = 25),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        #axis.text.y = element_blank(),
        #axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        panel.background = element_blank(), 
        panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
  )




Lipolysis_Markers2 %>% distinct() %>% 
  left_join(filter(Treatment_Response, Cluster == "Adipocyte")) %>% 
  ggplot(aes(x = avg_logFC, y = -log10(p_val_adj), color = LipoFunction)) + 
  
  geom_vline(xintercept = 0, color = "grey70") +
  
  geom_point(size = 3) + 
  geom_text_repel(aes(label = Gene), size = 4)  +
  ggtitle("lipolysis related genes in adipocytes") +
  
  scale_color_d3() +
  theme(text = element_text(family = "Lato", size = 25),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        #axis.text.y = element_blank(),
        #axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        panel.background = element_blank(), 
        panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
  )




p %>% save_x(data = ., name = paste0("Adipocyte_PTEN_ETC"), 1, 16, 9, svg = T)

p %>% save_x(data = ., name = paste0("Adipocyte_Lipolysis_Volcano"), 1, 16, 9, svg = T)
