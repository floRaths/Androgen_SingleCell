panel_5a <- function(){

look = theme(panel.background = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

A <- 
FeaturePlot(Sobj_integrated, reduction = "umap",
            features = c("Artery_1"),
            pt.size = 1.5, 
            order = T, 
            #max.cutoff = 4.2,
            cols = viridis::magma(n = 100),
            ncol = 1
) + look
B <- 
FeaturePlot(Sobj_integrated, reduction = "umap",
            features = c("Veins_1"),
            pt.size = 1.5, 
            order = T, 
            #max.cutoff = 4.2,
            cols = viridis::magma(n = 100),
            ncol = 1
) + look
C <- 
FeaturePlot(Sobj_integrated, reduction = "umap",
            features = c("Capillary_1"),
            pt.size = 1.5, 
            order = T, 
            #max.cutoff = 4.2,
            cols = viridis::magma(n = 100),
            ncol = 1
) + look


p <- A + B + C

}

p <- panel_5a()

p %>% save_x(data = ., name = paste0("BloodEC_Suptype_Module_Callin"), 1, 8, 16, svg = F)



panel_5b <- function(){

levels(Sobj_integrated) <- c("Vein", "Capillary", "Artery", "Immature", "Lymphatic_EC", "Lymphatic_EC_(SV2C+)", "Pericyte", "Vasc_SMC_(CLSTN2+)", "Vasc_SMC_(CACNB2+)")

query <- c("TLL1", "LNX1", "PCSK5", "RPL13", "RELN", "RADIL", "PDGFRB", "ACTA2", "THSD7B", "CLSTN2", "CACNB2")

p <- 
VlnPlot(Sobj_integrated, log = F, 
        features = query,
        #group.by = "Subcluster", 
        pt.size = 0, 
        ncol = 1
        )
}

p <- panel_5b()

p %>% save_x(data = ., name = paste0("BloodEC_Suptype_Specific_Markers"), 1, 6, 30, svg = T)



panel_5e <- function(){

percentile <- 
  read_tsv(paste0("Output/Reboot/SCENIC/expr_mat.adjacencies_all_cells.tsv")) %>% 
  group_by(TF) %>% 
  mutate(percentile_rank = ntile(importance, 100))

#percentile <- 
#  read_tsv(paste0("Output/Reboot/SCENIC/CellType_SCENIC/Output/expr_mat.adjacencies_Blood_EC.tsv")) %>% 
#  group_by(TF) %>% 
#  mutate(percentile_rank = ntile(importance, 100))


library(gprofiler2)
library(viridis)

clust = "PPARG"

query <- percentile %>% filter(TF == clust, percentile_rank > 95) %>% pull(target)

Results <- gost(query, organism = "hsapiens")


#p4 <- 
ggplot(head(filter(Results$result, source == "GO:BP"), 6), 
       aes(fct_reorder(term_name, p_value, .desc = T), intersection_size, fill = p_value)) + 
  geom_col() +
  
  scale_fill_viridis(begin = 0.4, end = 0.7, direction = -1) +
  
  ylab(label = "Number of Genes in Pathway") +
  
  ggtitle(paste(unique(head(filter(Results$result, source == "GO:BP"), 15)$source), "-", clust, "Module")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50), position = "top") +
  coord_flip() + 
  
  theme(text = element_text(family = "Lato", size = 20), 
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "grey30", size = 2, fill = NA),
        panel.grid = element_blank()
        )

}

p <- panel_5e()

p %>% save_x(data = ., name = paste0("PPARG_Module_GOBP"), 1, 10, 6, svg = T)



panel_5f <- function() {
  
  data %>% select(Sample, CellType, Subcluster, Type, query) %>% 
    pivot_longer(query, names_to = "Module") %>% 
    ggplot(aes(x = Subcluster, y = value, fill = Type)) + 
    
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
    
    theme(text = element_text(family = "Lato", size = 25),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none",
          panel.background = element_blank(), 
          panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
    )
  
}

query <- c("KEGG_PPAR_SIGNALING_PATHWAY")

p <- panel_5f()

p %>% save_x(data = ., name = paste0("PPARG_Signaling_Kegg"), 1, 10, 6, svg = T)




celltype = "Lymphatic_EC"

data_diff_CT <- readRDS(paste0("Output/MsigDB_Scores/CellType_Response_", celltype, ".rds"))

query <- c("REACTOME_VEGFR2_MEDIATED_CELL_PROLIFERATION", "WP_NOTCH_SIGNALING", "HALLMARK_PROTEIN_SECRETION", "REACTOME_ECM_PROTEOGLYCANS")


p <- 
data %>% 
  select(Sample, CellType, Subcluster, Type, query) %>% 
  filter(CellType %in% celltype) %>% 
  pivot_longer(query, names_to = "Module") %>% 
  
  ggplot(aes(x = CellType, y = value, fill = Type)) + 
  
  geom_split_violin(scale = "width", 
                    alpha = 1, 
                    trim = F,
                    lwd = 1.25, 
                    color = "grey30") + 
  
  geom_boxplot(aes(color = Type), 
               fill = "grey30", 
               width = 0.1, alpha = 1, 
               lwd = 1.25,
               position=position_dodge(0.7), outlier.alpha = 0) +
  
  #scale_fill_viridis_d(option = "C", begin = 0.4, end = 0.8) +
  scale_fill_manual(values = c("#A6499B", "#FAA42F")) +
  scale_color_manual(values = c("grey30", "grey30")) +
  
  stat_summary(fun = "median", 
               geom = "point", size = 3,
               position = position_dodge(0.7),
               color = "floralwhite") +
  
  facet_wrap(vars(Module), scales = "free_x") +
  #scale_y_sqrt() +
  theme(text = element_text(family = "Lato", size = 25),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), 
        panel.grid = element_blank(),
        legend.position = "bottom",
        panel.background = element_blank(), 
        panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
  ) + coord_flip()


p %>% save_x(data = ., name = paste0("LymphEC_Pathway_Violins"), 1, 16, 9, svg = T)  



plot_split <- function() {
  data %>% select(Sample, CellType, Subcluster, Type, query) %>% 
    filter(CellType %in% celltype) %>% 
    #filter(Subcluster %in% c("CD4_T", "CD8_T")) %>% 
    #scores %>% select(Subcluster, query) %>% 
    pivot_longer(query, names_to = "Module") %>% 
    #group_by(Sample, Subcluster, Module) %>% 
    #mutate(sample_avg = mean(value)) %>% 
    
    ggplot(aes(x = Subcluster, y = value, fill = Type)) + 
    
    geom_split_violin(scale = "width", 
                      alpha = 1, 
                      trim = F,
                      lwd = 1.25, 
                      color = "grey30") + 
    
    
    geom_boxplot(aes(color = Type), 
                 fill = "grey30", 
                 width = 0.1, alpha = 1, 
                 lwd = 1.25,
                 position=position_dodge(0.7), outlier.alpha = 0) +
    
    scale_fill_manual(values = c("#A6499B", "#FAA42F")) +
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


query= c("KEGG_VASCULAR_SMOOTH_MUSCLE_CONTRACTION", "HALLMARK_ANGIOGENESIS")

p <- plot_split()

p %>% save_x(data = ., name = paste0("Vasc_Acc_Pathway_Violins"), 1, 16, 9, svg = T)  





blood_ec_pathways <- function(){

processed <- readRDS(paste0("Output/MsigDB_Scores/Subcluster_Response_", celltype, ".rds"))

query <- c("REACTOME_HDACS_DEACETYLATE_HISTONES", "REACTOME_NOREPINEPHRINE_NEUROTRANSMITTER_RELEASE_CYCLE",
           "PID_ALPHA_SYNUCLEIN_PATHWAY", "PID_ALK2_PATHWAY",
           "REACTOME_RETROGRADE_NEUROTROPHIN_SIGNALLING", "REACTOME_SIGNALING_BY_BMP",
           "BIOCARTA_LEPTIN_PATHWAY", "KEGG_RIBOSOME",
           "REACTOME_INTERFERON_GAMMA_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
           "REACTOME_SIGNALING_BY_INTERLEUKINS", "WP_ESTROGEN_RECEPTOR_PATHWAY",
           "REACTOME_REGULATION_OF_TLR_BY_ENDOGENOUS_LIGAND", "KEGG_PPAR_SIGNALING_PATHWAY", "HALLMARK_ADIPOGENESIS",
           "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION")

mat <- 
  #data_diff %>% left_join(wilcox_Msig) %>% 
  processed %>% 
  distinct() %>% 
  filter(Module %in% query, CellType == celltype) %>% 
  pivot_longer(c(CF, TM), names_to = "Type") %>% 
  arrange(Type) %>% unite(Subcluster, Subcluster, Type) %>% ungroup( ) %>% 
  select(-Wilcox, -FC, -CellType) %>% 
  pivot_wider(names_from = Module, values_from = value) %>% 
  column_to_rownames("Subcluster") %>% t() %>% as.data.frame() %>% 
  rownames_to_column("Module") %>% 
  mutate(Module = tolower(Module), Module = str_replace_all(Module, "_", " ")) %>% 
  column_to_rownames("Module")

p <- 
pheatmap(mat, scale = "row",
         color = viridis::magma(n = 100, begin = 0, end = 0.9), 
         #color = my.colors,
         #breaks = my.breaks,
         clustering_distance_rows = "correlation", 
         cluster_cols = F, 
         #clustering_distance_cols = "correlation",
         
         fontsize_col = 20,
         fontsize_row = 20, 
         treeheight_row = 10, treeheight_col = 0, 
         #annotation_col = Ann,
         main = "Pathway avg_Score")

}

p <- blood_ec_pathways()

p %>% save_x(data = ., name = paste0("Blood_EC_Pathway_Heatmap"), 1, 12, 12, svg = T)  





sample_exclusion <- function(){
  
  margin = 0.5
  
  look = theme(panel.background = element_rect(fill  = "white", 
                                               color = "grey60", 
                                               size  = 1), 
               panel.grid = element_blank(), plot.margin = unit(c(margin,margin,margin,margin), "cm"),
               axis.title = element_blank(), 
               axis.text  = element_blank(),
               axis.ticks = element_blank(),
               legend.position = "none")
  
  celltype = "Blood_EC"
  
  Sobj_integrated <- AddMetaData(Sobj_integrated, metadata = as.data.frame(Sobj_integrated@reductions$umap@cell.embeddings))
  df_rna = Sobj_integrated@meta.data
  
  set.seed(467)
  rows_rna  <- sample(nrow(df_rna)) #, 10000)
  
  
  #P1 <- 
    df_rna[rows_rna, ] %>% #filter(CellType == "Myeloid") %>% 
    
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = factor(Sample, levels = levels_samp))) + 
    geom_point(alpha = 0.35, 
               size = 3, 
               stroke = 0) +
    
    #facet_wrap(~ CellType, nrow = 1) +
    
    scale_color_manual(values = cols_samp) + 
    #scale_color_manual(values = colors[1:4]) + 
    #scale_color_manual(values = colors[7:10]) + 
    
    look + 
    theme(strip.background = element_blank(), panel.background = element_blank(),
          strip.text = element_blank())
  
}

p <- sample_exclusion()

p %>% save_x(data = ., name = paste0("Blood_EC_sample_exclusion"), 0.75, 12, 12, svg = F)  








