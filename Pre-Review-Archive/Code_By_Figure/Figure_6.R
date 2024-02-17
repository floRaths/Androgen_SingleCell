
panel_S6b <- function(){

avg_lym    <- readRDS("Output/Lymphoid_Subcluster_avgExpression.rds")
marks_lym  <- readRDS("Output/Lymphoid_Subcluster_Markers.rds")
avg_mye    <- readRDS("Output/Myeloid_Subcluster_avgExpression.rds")
marks_mye  <- readRDS("Output/Myeloid_Subcluster_Markers.rds")

genes1 <- 
  marks_lym %>% 
  filter(!is.na(avg_logFC), p_val_adj <= 0.05, pct.1 > 0.25) %>% 
  filter(!str_detect(gene, "-AS1|-AS2|\\.|LINC0")) %>% 
  select(-p_val, -p_val_adj) %>% 
  group_by(cluster) %>% 
  top_n(3, avg_logFC) %>% 
  arrange(cluster) %>% pull(gene)



mat1 <- 
  avg_lym %>% 
  rownames_to_column("gene") %>% 
  filter(gene %in% c(genes1, "NCAM1", "CD8A", "CD3E", "CD4", "TOX", "IL7R", "GZMB", "LAG3")) %>% 
  #mutate_at(colnames(avg_lym), percent_rank) %>% 
  column_to_rownames("gene")



#mat1 <- 
#  avg_lym %>% 
#  rownames_to_column("gene") %>% 
#  filter(gene %in% c(genes1, "NCAM1", "CD8A", "CD3E", "CD4", "TOX", "IL7R", "GZMB", "LAG3")) %>% 
#  column_to_rownames("gene") %>% t() %>%  as.data.frame() %>%  rownames_to_column("Name") %>%
#    mutate_at(c(genes1, "NCAM1", "CD8A", "CD3E", "CD4", "TOX", "IL7R", "GZMB", "LAG3"), percent_rank) %>% 
#    column_to_rownames("Name") %>% t()
  

p1 <- 
pheatmap(mat1, 
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation", 
         clustering_method = "ward.D2",
         scale = "row", 
         color = colorRampPalette(brewer.pal(11, "BrBG"))(100)
         #color = viridis::magma(n = 100, end = 0.9)
         )




genes2 <- 
  marks_mye %>% 
  filter(!is.na(avg_logFC), p_val_adj <= 0.05, pct.1 > 0.25) %>% 
  filter(!str_detect(gene, "-AS1|-AS2|\\.|LINC0")) %>% 
  select(-p_val, -p_val_adj) %>% 
  group_by(cluster) %>% 
  top_n(3, avg_logFC) %>% 
  arrange(cluster) %>% pull(gene)



mat2 <- avg_mye %>% 
  rownames_to_column("gene") %>% 
  filter(gene %in% c(genes2, "MRC1", "CD163", "ITGAX", "ITGAM", "C3", "PALD1", "HLA-DRB1", "CD36", "SNTB1")) %>% 
  column_to_rownames("gene")


#mat2 <- avg_mye %>% 
#  rownames_to_column("gene") %>% 
#  filter(gene %in% c(genes2, "MRC1", "CD163", "ITGAX", "ITGAM", "C3", "PALD1")) %>% 
#  mutate_at(colnames(avg_mye), percent_rank) %>% 
#  column_to_rownames("gene")


p2 <- pheatmap(mat2, 
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation", 
         clustering_method = "ward.D2",
         color = colorRampPalette(brewer.pal(11, "BrBG"))(100)
         )


as.ggplot(p1) + as.ggplot(p2)

}

p <- panel_S6b()

p %>% save_x(data = ., name = "Immune_Marker_Calling_Heatmap", 1, 16, 10, svg = T)
