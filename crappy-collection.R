mat <- 
Sobj_integrated@meta.data %>% 
  rownames_to_column("X1") %>% 
  select(X1, Type, CellType) %>% 
  left_join(select(bind, X1, names)) %>% 
  pivot_longer(names) %>% 
  group_by(name, CellType, Type) %>% 
  summarise(avg = mean(value)) %>% 
  unite(CellType, CellType, Type) %>% 
  pivot_wider(names_from = CellType, values_from = avg) %>% mutate(name = tolower(name), name = str_replace_all(name, "_", " ")) %>% 
  column_to_rownames("name")




wilcox_Msig <- 
Sobj_integrated@meta.data %>% 
  rownames_to_column("X1") %>% 
  select(X1, Type, CellType) %>% 
  left_join(select(bind, X1, names)) %>% 
  as_tibble() %>% pivot_longer(names, names_to = "Module") %>% 
  group_by(CellType, Module) %>% 
  
  do(w = wilcox.test(value~Type, data = ., paired = F)) %>% 
  summarise(CellType, Module, Wilcox = w$p.value)







Sobj_integrated@meta.data %>% 
  rownames_to_column("X1") %>% 
  select(X1, Type, CellType) %>% 
  left_join(select(bind, X1, names)) %>% 
  pivot_longer(names) %>% 
  group_by(name, CellType, Type) %>% 
  summarise(avg = mean(value)) %>% 
  pivot_wider(names_from = Type, values_from = avg) %>% 
  mutate(diff = TM - CF) %>% 
  left_join(wilcox_Msig, by = c("CellType", "name" = "Module")) %>% 
  filter(Wilcox <= 0.05, abs(diff) > 0.05)



p <- 
pheatmap(mat, #scale = "row",
         color = viridis::magma(n = 100, begin = 0, end = 1), 
         #color = my.colors,
         #breaks = my.breaks,
         clustering_distance_rows = "correlation", 
         cluster_cols = F, 
         fontsize_col = 10,
         fontsize_row = 10, 
         treeheight_row = 10, 
         treeheight_col = 0, 
         main = "Pathway avg_Score")


p %>% save_x(data = ., name = paste0("Cancer_Pathways2"), 1, 16, 12, svg = T)  



####################################################
####################################################
####################################################



df1 <- 
  Bind %>% 
  filter(GeneRna %in% gmt$HALLMARK_INTERFERON_GAMMA_RESPONSE) %>% 
  select(GeneRna, Motif) %>% 
  distinct() %>% 
  group_by(Motif) %>% 
  summarise(n = n()) %>% arrange(-n) %>% 
  mutate(Module = "IFNG", perc = (n/172)*100)




df2 <- 
  Bind %>% 
  filter(GeneRna %notin% gmt$HALLMARK_INTERFERON_GAMMA_RESPONSE) %>% 
  select(GeneRna, Motif) %>% 
  distinct() %>% 
  group_by(Motif) %>% 
  summarise(n = n()) %>% arrange(-n) %>% 
  mutate(Module = "Others", perc = (n/12473)*100)


bind_rows(df1, df2) %>% select(-n) %>% 
  pivot_wider(names_from = Module, values_from = perc) %>% 
  mutate(diff = IFNG-Others, avg = mean(diff), adj = Others + avg)




#%>% filter(perc > 30) %>% ggplot(aes(y = reorder(Motif, perc), x = perc, color = Module)) + geom_point()



look <- theme(axis.text = element_blank(), axis.title = element_blank(),
              axis.ticks = element_blank())

p <- AR + ESR1 + PGR

p <- AR

p %>% save_x(data = ., name = "UMAPS_Legend", 1, 6, 6, svg = T)



####################################################
####################################################
####################################################



AtacSeqEmbeddingBasedOnPeakMatrix_20201007 <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/ATAC_Integration/AtacSeqEmbeddingBasedOnPeakMatrix_20201007.rds")



plot_prop <- function() {
  
  df <- Sobj_integrated@meta.data %>% select(CellType, Sample, Type, Subcluster) #%>% filter(CellType %in% c("LUM_HR-pos", "LUM_HR-neg", "Basal"))
  
  df$CellType <- as.character(df$CellType)
  df$Subcluster <- as.character(df$Subcluster)
  
  
  Prop <- 
    as.data.frame(prop.table(table(df$Subcluster, 
                                   df$Sample), 
                             margin = 2)) %>% 
    dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq") %>% 
    mutate(Type = substring(.$Sample, 1, 2))# %>% arrange(factor(Cluster, levels = levels), desc(Freq))
  
  #P2 <- 
  Prop %>% 
    ggplot(aes(x = reorder(Cluster, Freq), y = Freq, fill = Type)) + 
    geom_boxplot() + 
    #geom_point(aes(color = Sample), position = position_dodge(width = 0.5)) + 
    ggtitle("RNA Data") + 
    scale_fill_manual(values = c("#9d3396", "#eb933e")) +
    theme(text = element_text(family = "Lato", size = 25),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          #axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          panel.background = element_blank(), 
          panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
          panel.grid = element_line(color = "grey90")) 
  
  
  
  
  df <- AtacSeqEmbeddingBasedOnPeakMatrix_20201007 %>% select(CellType = "predictedCellTypeGroup", Sample = "Sample", Type = "SampleType") %>% filter(CellType %in% c("LUM_HR-pos", "LUM_HR-neg", "Basal"))
  
  df$CellType <- as.character(df$CellType)
  
  Prop <- 
    as.data.frame(prop.table(table(df$CellType, 
                                   df$Sample), 
                             margin = 2)) %>% 
    dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq") %>% 
    mutate(Type = substring(.$Sample, 1, 2))# %>% arrange(factor(Cluster, levels = levels), desc(Freq))
  
  P1 <- 
    Prop %>% 
    ggplot(aes(x = Cluster, y = Freq, fill = Type)) + 
    geom_boxplot() + 
    ggtitle("ATAC Data") + 
    #geom_point(aes(color = Sample)) + 
    scale_fill_manual(values = c("#9d3396", "#eb933e")) +
    theme(text = element_text(family = "Lato", size = 25),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          #axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none",
          panel.background = element_blank(), 
          panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
          panel.grid = element_line(color = "grey90")) 
  
  
  
  P2 + ylim(c(0,0.65)) + P1 + ylim(c(0,0.65))
  
  
}

plot_prop()







source_data_MaxEcc0_7_MinGlands5 %>% group_by(Type, sample) %>% summarise(mean = n()) %>% 
  ggplot(aes(x = Type, y = mean)) + geom_boxplot()



tissue_areas <- read_csv("utilities/tissue_areas.csv")
source_data_means_MaxEcc0_7_MinGlands5 <- read_csv("utilities/source_data_means_MaxEcc0.7_MinGlands5.csv")

vars <- source_data_means_MaxEcc0_7_MinGlands5 %>% colnames()

look <- theme(text = element_text(family = "Lato", size = 25),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              #axis.text.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.position = "bottom",
              panel.background = element_blank(), 
              panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
              panel.grid = element_line(color = "grey90")) 



P1 <- 
  source_data_means_MaxEcc0_7_MinGlands5 %>% 
  ggplot(aes(x = Type, y = log_convex_area, fill = Type)) + 
  geom_boxplot(width = 0.5, size = 1.5) + 
  geom_point(position = position_jitter(width = 0.2), size = 4, alpha = 0.5) + 
  
  ggtitle("log_convex_area") + 
  scale_fill_manual(values = c("#9d3396", "#eb933e")) + look

P2 <- 
  source_data_means_MaxEcc0_7_MinGlands5 %>% 
  ggplot(aes(x = Type, y = glands_per_area, fill = Type)) + 
  geom_boxplot(outlier.alpha = 0, width = 0.5, size = 1.5) + 
  geom_point(position = position_jitter(width = 0.2), size = 4, alpha = 0.5) + 
  
  ggtitle("glands_per_area") + 
  scale_fill_manual(values = c("#9d3396", "#eb933e")) + look

P3 <- 
  source_data_means_MaxEcc0_7_MinGlands5 %>% 
  ggplot(aes(x = Type, y = percent_interior_cells, fill = Type)) + 
  geom_boxplot(width = 0.5, size = 1.5) + 
  geom_point(position = position_jitter(width = 0.2), size = 4, alpha = 0.5) + 
  
  ggtitle("percent_interior_cells") + 
  scale_fill_manual(values = c("#9d3396", "#eb933e")) + look

P4 <- 
  source_data_means_MaxEcc0_7_MinGlands5 %>% 
  ggplot(aes(x = Type, y = log_nuclei, fill = Type)) + 
  geom_boxplot(width = 0.5, size = 1.5) + 
  geom_point(position = position_jitter(width = 0.2), size = 4, alpha = 0.5) + 
  
  ggtitle("log_nuclei") + 
  scale_fill_manual(values = c("#9d3396", "#eb933e")) + look



p <- P1 + P2 + P3 + P4 + patchwork::plot_layout(nrow = 1)

p %>% save_x(data = ., name = paste0("Epithelial_Acini_Features"), 1, 16, 9, svg = T)    





####################################################
####################################################
####################################################




mat_query <- c("REACTOME_SMOOTH_MUSCLE_CONTRACTION", "REACTOME_TYPE_I_HEMIDESMOSOME_ASSEMBLY", "REACTOME_CELL_JUNCTION_ORGANIZATION",
               #"PID_RHOA_PATHWAY", 
               #"REACTOME_MYOGENESIS", 
               #"REACTOME_NUCLEAR_RECEPTOR_TRANSCRIPTION_PATHWAY", 
               #"BIOCARTA_EIF_PATHWAY", 
               #"REACTOME_RESPONSE_TO_METAL_IONS", 
               #"HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
               #"KEGG_RIBOSOME", 
               #"REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
               "REACTOME_ATTENUATION_PHASE", 
               "WP_HEMATOPOIETIC_STEM_CELL_DIFFERENTIATION" 
               #"PID_HDAC_CLASSII_PATHWAY"
)

basal_motifs <- 
  R2_values_motifZscoresVersusCellEnrichmentInPathways %>% colnames() %>% as_tibble() %>% 
  separate(value, into = c("A", "B"), remove = F)

basal_filtered <- 
  Treatment_Response %>% filter(Gene %in% basal_motifs$A, Cluster == "Basal", pct.1 > 0.2 | pct.2 > 0.2, abs(avg_logFC) > 0.2)


pathways <- 
  R2_values_motifZscoresVersusCellEnrichmentInPathways %>% rownames() %>% as_tibble() %>% 
  separate(value, into = c("A", "B"), sep = " Cis| Trans", remove = F) %>% 
  filter(A %in% mat_query) %>% pull(value)

anno <- R2_values_motifZscoresVersusCellEnrichmentInPathways %>% rownames() %>% as_tibble() %>% 
  separate(value, into = c("A", "B"), sep = " Cis| Trans", remove = F) %>% 
  filter(A %in% mat_query) %>% select(-A) %>% column_to_rownames("value")



p <- pheatmap(R2_values_motifZscoresVersusCellEnrichmentInPathways[pathways, filter(basal_motifs, A %in% basal_filtered$Gene)$value], 
              color = viridis::viridis(n = 100), annotation_row = anno)

p %>% save_x(data = ., name = paste0("NR4A1_Scatter"), 1, 10, 7, svg = T)    





mat <- 
  motifZscoresVersusCellEnrichmentInPathways %>% 
  unite(X, c(Pathway, Subcluster)) %>% separate(Motif, into = c("Gene", "B"), remove = F) %>% 
  left_join(filter(Treatment_Response, Cluster == "Basal")) %>% 
  filter(Gene %in% basal_motifs$A, Cluster == "Basal", pct.1 > 0.2 | pct.2 > 0.2, abs(avg_logFC) > 0.2) %>% 
  
  select(Motif, X , R2) %>% 
  distinct() %>% 
  pivot_wider(names_from = Motif, values_from = R2) %>% column_to_rownames("X")

pheatmap(mat)


