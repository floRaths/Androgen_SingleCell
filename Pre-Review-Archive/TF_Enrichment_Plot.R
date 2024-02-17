Motif_Summary <- Motif_Summary %>% filter(dtbs == "CISBP") %>% select(Motif, TF) %>% separate(Motif, into = c("Dat", "Motif"), sep = ".cisBP.") %>% left_join(Treatment_Response_Subcluster, by = c("TF" = "Gene"))

Motif_Summary_Complete_Unfiltered <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/ATAC_Integration/Motif_Summary_Complete_Unfiltered.rds")
Enrichment <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/ATAC_Integration/Motif_Enrichment_RDS_Files/EnrichmentOfExpressedTfsOrMotifsInDifferentCellTypes_onlyCisBp.rds")

cols <- c("#9d3396", "#eb933e")


celltype <- "Fibroblast"

motif_exp <- Motif_Summary_Complete_Unfiltered %>% 
  filter(dtbs == "CISBP", CellType == celltype, perc_CF >= 5 | perc_TM >= 5, avg_CF >= 0.5 | avg_TM >= 0.5) %>% 
  separate(Motif, into = c("A", "Motif"), sep = "otif.cisBP.") %>% pull(Motif)


motif_exp <- Motif_Summary_Complete_Unfiltered %>% 
  filter(dtbs == "CISBP", CellType == celltype, !is.na(avg_logFC), perc_CF >= 5 | perc_TM >= 5, avg_CF >= 0.5 | avg_TM >= 0.5) %>% 
  separate(Motif, into = c("A", "Motif"), sep = "otif.cisBP.") %>% pull(Motif)


df <- 
  Enrichment %>% as_tibble() %>% unique() %>% 
  separate(Comparison, into = c("CellType", "Type", "C", "D", "E", "F", "G"), sep = " ", extra =  "merge") %>% 
  select(-C, -D, -E, -G, -"F") %>% filter(CellType == celltype, TF %in% motif_exp) %>% 
  mutate(Enrichment_Div = ifelse(Type == "Trans", Enrichment, Enrichment*-1)) %>% 
  mutate(xlim = plyr::round_any(max(Enrichment), 0.1, ceiling))


motifs <- df %>% group_by(Type) %>% top_n(-10, rank) %>% pull(TF) %>% unique()
xlim   <- df %>% filter(TF %in% motifs) %>% pull(xlim) %>% unique()
  
  
#image <- 
df %>% filter(TF %in% motifs) %>% 

  ggplot(aes(x = reorder(TF, Enrichment_Div), y = Enrichment_Div, fill = Type)) + 
  geom_col(width = 0.85, alpha = 1) + 
  coord_flip() +
  scale_fill_manual(values = cols) +
  
  ggtitle(label = paste0(celltype, " Motif Enrichment")) + 
  
  #ylim(-xlim, xlim) +
  
  geom_hline(yintercept = 0, lwd = 3, color = "white") +
  
  theme(text = element_text(family = "Lato", face = "bold", size = 15),
        legend.position = "bottom", 
        
        title = element_text(size = 16, colour = "grey20"),
        
        axis.text.x = element_text(size = 15),
        
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y  = element_blank(),
        panel.grid.major.y  = element_line(color = "grey70", 
                                           lineend = "round",
                                           linetype = "dotted", 
                                           size = 0.6),
        
        panel.background = element_rect(fill = "white", 
                                        color = "grey70", 
                                        size = 2),
      
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

        
save(name = paste0(celltype, " Motif Enrichment"), 0.8, 6, 8, svg = F)



df <- 
Treatment_Response_Subcluster %>% 
  filter(Gene %in% tfs, CellType == celltype) %>% 
  select(Gene, avg_logFC, Cluster) %>% 
  pivot_wider(names_from = Cluster, values_from = avg_logFC)

df[match(tfs, df$Gene),]  


#mat[is.na(mat)] <- 0


#image <- 
pheatmap(mat, 
         #scale = "row", 
         fontsize = 10, cluster_cols = F,
         cluster_rows = F,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "euclidean",
         #annotation_row = Ann,
         #annotation_colors = ann_colors, 
         treeheight_row = 0, treeheight_col = 0, 
         #cutree_rows = 1, 
         #cutree_cols = 2,
         fontsize_col = 15,
         fontsize_row = 10, 
         #labels_col = NULL, 
         color = viridis::magma(n = 100, begin = 0, end = 0.9), 
         main = "Pathway logFC"
)








levels <- c("LUM_HR-pos", "LUM_HR-neg", "Basal", "Fibroblast", "Adipocyte", "Blood_EC", "Lymphatic_EC", "Vasc_Accessory", "Myeloid", "Lymphoid")
Marker_list    <- vector("list", length = length(levels))




for (i in 1:length(levels)) {
  

celltype <- levels[i]

motif_exp <- Motif_Summary_Complete_Unfiltered %>% 
  filter(dtbs == "CISBP", CellType == celltype, perc_CF >= 5 | perc_TM >= 5, avg_CF >= 0.5 | avg_TM >= 0.5) %>% 
  separate(Motif, into = c("A", "Motif"), sep = "otif.cisBP.") %>% pull(Motif)


df <- 
  Enrichment %>% as_tibble() %>% unique() %>% 
  separate(Comparison, into = c("CellType", "Type", "C", "D", "E", "F", "G"), sep = " ", extra =  "merge") %>% 
  select(-C, -D, -E, -G, -"F") %>% filter(CellType == celltype, TF %in% motif_exp) %>% 
  mutate(Enrichment_Div = ifelse(Type == "Trans", Enrichment, Enrichment*-1)) %>% 
  mutate(xlim = plyr::round_any(max(Enrichment), 0.1, ceiling))


motifs_rank <- df %>% group_by(Type) %>% top_n(-10, rank) %>% pull(TF) %>% unique()




motif_exp <- Motif_Summary_Complete_Unfiltered %>% 
  filter(dtbs == "CISBP", CellType == celltype, !is.na(avg_logFC), perc_CF >= 5 | perc_TM >= 5, avg_CF >= 0.5 | avg_TM >= 0.5) %>% 
  separate(Motif, into = c("A", "Motif"), sep = "otif.cisBP.") %>% pull(Motif)

df <- 
  Enrichment %>% as_tibble() %>% unique() %>% 
  separate(Comparison, into = c("CellType", "Type", "C", "D", "E", "F", "G"), sep = " ", extra =  "merge") %>% 
  select(-C, -D, -E, -G, -"F") %>% filter(CellType == celltype, TF %in% motif_exp) %>% 
  mutate(Enrichment_Div = ifelse(Type == "Trans", Enrichment, Enrichment*-1)) %>% 
  mutate(xlim = plyr::round_any(max(Enrichment), 0.1, ceiling))


motifs_FC <- df %>% group_by(Type) %>% top_n(-10, rank) %>% pull(TF) %>% unique()


Marker_list[[i]] <- c(motifs_rank, motifs_FC) %>% unique() %>% as.data.frame() %>% rename(Motif = 1) %>% mutate(CellType = celltype)
}


Bind <- do.call(rbind.data.frame, Marker_list)
