library(venneuler)
library(VennDiagram)

AER <- 
adjacencies %>% filter(TF %in% c("AR", "ESR1", "PGR", "ANKRD30A", "CUX2", 
                                 "TRPS1", "CUX2", "NR4A1", "PPARG",  
                                 "GPAM", "IKZF1", "FOSB", 
                                 "TP63", "CREB5", "RPS4X", "MYLK", "MLXIPL",
                                 "TBX1", "ATF3", "PRRX1", "PRDM16",
                                 "NFE2L3", "RBPJ", "FOXO3", "YWHAZ",  "THRB", "SMARCA4", "SNRNP70", "ZBTB7A")) %>% 
  arrange(-importance)# %>% group_by(TF) %>% top_n(100, importance)












conftargets <- 
adjacencies %>% 
  group_by(TF) %>% 
  mutate(percentile_rank = ntile(importance, 100)) %>% 
  filter(percentile_rank >= 90) %>% 
  pull(target) %>% 
  unique()







mat <- 
AER %>% filter(target %in% conftargets) %>% 
  pivot_wider(names_from = TF, values_from = importance) %>% column_to_rownames("target")


mat <- 
percentile %>% filter(percentile_rank > 90, TF %in% c("CUX2", "ZNF689", "ESR1", "PGR", 
                                                      "AR", "TRPS1", "ANKRD30A","ANKRD30B", "HNF4G", "MYB", "TFAP2A",
                                                      "TBX3", "BATF",  "HES1", "XBP1", "BHLHE40", 
                                                      "GATA3", "KDM4B")) %>% 
  dplyr::select(-percentile_rank) %>% 
  pivot_wider(names_from = TF, values_from = importance) %>% 
  column_to_rownames("target") 





celltype <- "Fibroblast"


expressed_TFs <- 
  Average_TF_Expression %>% ungroup() %>% 
  filter(CellType_All == celltype, expression == "expressed", perc >= 30) %>% 
  arrange(Type, -perc) %>% pull(TF) %>% unique() %>% sort()

TF_candidates <- 
  Marker_Bind %>% filter(CellType == celltype) %>% 
  pull(Gene) %>% 
  unique() %>% 
  intersect(expressed_TFs)


adjacencies <- read_tsv(paste0("~/Documents/SCENIC/pySCENIC_Full/scenicdata/CT_All_4kvar/Results/", celltype, "/expr_mat.adjacencies.tsv"))
adjacencies <- read_tsv("~/Documents/SCENIC/pySCENIC_Full/scenicdata/TFs_VarFeatures/Results_4k/expr_mat.adjacencies.tsv")

percentile <- 
  adjacencies %>% 
  group_by(TF) %>% 
  mutate(percentile_rank = ntile(importance, 100))


fc <- 0.25

targets <- 
percentile %>% 
  filter(TF %in% TF_candidates) %>% 
  group_by(TF) %>% 
  filter(percentile_rank > 95) %>% 
  #group_by(target) %>% 
  #top_n(5, importance) %>% 
  arrange(target) %>% 
  pull(target) %>% unique() %>% 
  intersect(filter(Marker_Bind, avg_logFC >= fc | avg_logFC <= -fc, CellType == celltype)$Gene) %>% 
  unique()




#TFs <- c("AR", "ESR1", "PGR", "ANKRD30A", "CUX2", "TRPS1", "BATF", "MYB", "ZNF689", "ESRRG", "GATA3", "KDM4B", "TBX3", "CERS4", "NFIA")

#x <- percentile %>% filter(TF == "ANKRD30A", percentile_rank > 99)
#TFs <- percentile %>% filter(target %in% x$target) %>% group_by(target) %>% top_n(3, importance) %>% pull(TF) %>% unique() %>% sort()

#targets <- filter(Marker_Bind, avg_logFC >= 0.25 | avg_logFC <= -0.25, CellType == celltype)$Gene

mat <- 
  percentile %>% filter(TF %in% TF_candidates, target %in% targets) %>% 
  dplyr::select(-4) %>%  
  pivot_wider(names_from = TF, values_from = importance) %>% 
  column_to_rownames("target")

mat[is.na(mat)] <- 0

image <- 
  pheatmap(t(mat), scale = "column", fontsize = 10, 
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           color = viridis::magma(n = 100),
           annotation_col = column_to_rownames(dplyr::select(filter(Marker_Bind, CellType == celltype), 1, 3), "Gene"),
           #annotation_colors = ann_colors, 
           treeheight_row = 0, treeheight_col = 0, 
           cutree_rows = 1, 
           cutree_cols = 1,
           fontsize_row = 10,
           fontsize_col = 8, 
           show_colnames = T
            
           )

  
Ann <-   
column_to_rownames(mean_logFC, "Gene")

  
  
  
query <-   
percentile %>% filter(percentile_rank > 90, TF %in% c("TRPS1")) %>% pull(target)
  
reactomeBar(name = "Test", query, 15)
  







  
  
  
  VENN.LIST <- list(filter(AER, TF == "AR")$target, 
                    filter(AER, TF == "ESR1")$target, 
                    filter(AER, TF == "PGR")$target,
                    filter(AER, TF == "ANKRD30A")$target
  )
  
  venn.plot <- venn.diagram(VENN.LIST , NULL, 
                            fill=c("darkmagenta", "darkblue", "Orange"), 
                            alpha=c(0.5, 0.5, 0.5), 
                            cex = 2, 
                            cat.fontface=4, 
                            category.names=c("AR", "ESR", "PGR"), 
                            main="Hormone Receptor Overlap")
  
  # To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
  grid.newpage();
  grid.draw(venn.plot)
  
  
  
  
  


















  