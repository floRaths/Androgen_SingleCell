GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}




celltype = "Myeloid"

scores <- read_csv(paste0("Output/MsigDB_Scores/gene_sets_GO_scores/", celltype, "_scores.csv"))
#scores <- read_csv(paste0("Output/MsigDB_Scores/all_cells_scores.csv"))

gene_counts <- readRDS("Output/MsigDB_Scores/Module_Gene_Counts.rds")
modules_keep <- gene_counts %>% filter(between(value, 5, 1000)) %>% pull(Module)
names <- scores %>% colnames() %>% intersect(modules_keep)


data <- Sobj_integrated@meta.data %>% 
  rownames_to_column("X1") %>% 
  select(X1, Sample, Type, CellType, Subcluster) %>% 
  mutate_all(as.character) %>% 
  inner_join(select(scores, X1, names), by = c("X1"))


# Create FC test across subclusters
##################################################

#contains(c("REACTOME", "BIOCARTA", "PID", "KEGG", "HALLMARK"))

cluster <- data %>% pull(Subcluster) %>% unique()
Marker_list    <- vector("list", length = length(cluster))

for (i in seq_along(Marker_list)) {
  
  Marker_list[[i]] <- 
    
    data %>% 
    select(Type, Subcluster, CellType, names) %>% 
    as_tibble() %>% pivot_longer(names, names_to = "Module") %>% 
    mutate(Test = Subcluster, Group = ifelse(Test == cluster[i], "group1", "group2")) %>% 
    
    group_by(Group, Module) %>% mutate(group_avg = mean(value)) %>% 
    select(Module, Group, group_avg) %>% 
    distinct() %>% 
    pivot_wider(names_from = Group, values_from = group_avg) %>% 
    mutate(FC = log2(group1/group2)) %>% 
    add_column(CLuster = cluster[i])
  
  #group_by(Module) %>% do(w = wilcox.test(value~Group, data = ., paired = F)) %>% 
  #summarise(Module, Wilcox = w$p.value) %>% add_column(Subcluster = cluster[i])
}


Bind <- do.call(rbind.data.frame, Marker_list)

Bind %>% saveRDS(paste0("Output/MsigDB_Scores/Subcluster_Markers_", celltype, ".rds"))


plot_vln <- function() {
  data %>% select(Sample, CellType, Subcluster, Type, query) %>% 
    #scores %>% select(Subcluster, query) %>% 
    pivot_longer(query, names_to = "Module") %>% 
    #group_by(Sample, Subcluster, Module) %>% 
    #mutate(sample_avg = mean(value)) %>% 
    
    ggplot(aes(x = Subcluster, y = value, fill = Subcluster)) + 
    
    geom_violin(scale = "width", 
                alpha = 0.75, 
                trim = F,
                lwd = 1.25, 
                color = "grey30") + 
    
    
    geom_boxplot(color = "grey30", 
                 fill = "grey30", 
                 width = 0.1, alpha = 1, 
                 lwd = 1.25,
                 position=position_dodge(0.7), outlier.alpha = 0) +
    
    #scale_fill_viridis_d(option = "C", begin = 0.4, end = 0.8) +
    #scale_color_manual(values = c("grey30", "grey30")) +
    
    stat_summary(fun = "median", 
                 geom = "point", size = 1.5,
                 position = position_dodge(0.7),
                 color = "floralwhite") +
    
    facet_wrap(vars(Module), scales = "free_y") +
    
    #scale_y_sqrt() +
    
    theme(text = element_text(family = "Lato"),
          axis.text = element_text(face = "bold"),
          panel.background = element_rect(fill = "gray96"),
          strip.text = element_text(size = 13, face = "bold", hjust = 0), 
          legend.key.size = unit(1,"cm"), legend.justification = c(1,1))
}

query <- Bind %>% ungroup() %>% filter(CLuster == "Lymphatic_EC", group1 > 0.2, str_detect(Module, "REACTOME|BIOCARTA|WP_|PID|KEGG|HALLMARK")) %>% top_n(20, FC) %>% arrange(-group1) %>% pull(Module)



plot_vln()



# Basal
#mat_query <- c("REACTOME_SMOOTH_MUSCLE_CONTRACTION", "REACTOME_TYPE_I_HEMIDESMOSOME_ASSEMBLY", "REACTOME_CELL_JUNCTION_ORGANIZATION",
#               "PID_RHOA_PATHWAY", "REACTOME_MYOGENESIS", "REACTOME_NUCLEAR_RECEPTOR_TRANSCRIPTION_PATHWAY", "BIOCARTA_EIF_PATHWAY", 
#               "REACTOME_RESPONSE_TO_METAL_IONS", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "KEGG_RIBOSOME", "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION")



mat <- 
  data %>% select(Subcluster, mat_query) %>% 
  pivot_longer(mat_query, names_to = "Module") %>% 
  group_by(Subcluster, Module) %>% 
  mutate(avg = mean(value)) %>% select(-value) %>% 
  ungroup() %>%  distinct() %>% 
  pivot_wider(names_from = Subcluster, values_from = avg) %>% 
  column_to_rownames("Module")

image <- 
  pheatmap(mat, scale = "row",
           color = viridis::magma(n = 100, begin = 0, end = 0.9), 
           clustering_distance_rows = "correlation", #cluster_cols = F, 
           clustering_distance_cols = "correlation",
           
           fontsize_col = 15,
           fontsize_row = 12, 
           treeheight_row = 10, treeheight_col = 0, 
           #annotation_col = Ann,
           main = "Pathway avg_Score")


save(name = paste0("Msig_Subcluster_Markers_", celltype), 0.8, 10, 5)






# Test per subcluster
##################################################


wilcox_Msig <- 
  data %>% 
  select(Type, Subcluster, CellType, names) %>% 
  as_tibble() %>% pivot_longer(names, names_to = "Module") %>% 
  group_by(Subcluster, Module) %>% 
  
  do(w = wilcox.test(value~Type, data = ., paired = F)) %>% 
  summarise(Subcluster, Module, Wilcox = w$p.value)



sig_modules <- wilcox_Msig %>% filter(Wilcox <= 0.05) %>% pull(Module) %>% unique()

data_avg <- 
  data %>% 
  select(Type, CellType, Subcluster, sig_modules) %>% 
  as_tibble() %>% group_by(Subcluster, Type) %>% 
  mutate_at(vars(sig_modules), mean) %>% ungroup %>% 
  distinct()


#data_avg <- 
#  data %>% 
#  select(Type, CellType, Subcluster, names) %>% 
#  as_tibble() %>% group_by(Type) %>% 
#  mutate_at(vars(names), mean) %>% ungroup %>% 
#  distinct()


data_diff <- 
  data_avg %>% 
  pivot_longer(sig_modules, names_to = "Module") %>% 
  pivot_wider(names_from = Type, values_from = value) %>% 
  group_by(Subcluster, Module) %>% 
  mutate(FC = log2(TM/CF))


data_diff %>% left_join(wilcox_Msig) %>% saveRDS(paste0("Output/MsigDB_Scores/gene_sets_GO_scores/Subcluster_Response_", celltype, ".rds"))

processed <- readRDS(paste0("Output/MsigDB_Scores/gene_sets_GO_scores/Subcluster_Response_", celltype, ".rds"))



plot_split <- function() {
  data %>% select(Sample, CellType, Subcluster, Type, query) %>% filter(CellType %in% celltype) %>% 
    #scores %>% select(Subcluster, query) %>% 
    pivot_longer(query, names_to = "Module") %>% 
    #group_by(Sample, Subcluster, Module) %>% 
    #mutate(sample_avg = mean(value)) %>% 
    
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
    
    #scale_y_sqrt() +
    
    theme(text = element_text(family = "Lato"),
          axis.text.x = element_text(face = "bold", size = 15, angle = -25, hjust = 0, vjust = 1),
          axis.title.x = element_blank(),
          panel.background = element_rect(fill = "gray96"),
          strip.text = element_text(size = 13, face = "bold", hjust = 0), 
          legend.key.size = unit(1,"cm"), legend.justification = c(1,1))
  
}

#data_diff %>% left_join(wilcox_Msig) %>% 
processed %>% 
  filter(Wilcox <= 0.05, FC < -0.1, CF > 0.1,
         #str_detect(Module, "REACTOME|BIOCARTA|WP_|PID|KEGG|HALLMARK"), #Subcluster == 0,
         str_detect(Module, "PHAGO")
  ) %>% 
  #filter(TM > 0) %>% 
  ungroup() %>% top_n(-40, Wilcox) %>% arrange(FC)


#query <- c("REACTOME_TYPE_I_HEMIDESMOSOME_ASSEMBLY",
#           "BIOCARTA_GATA3_PATHWAY", "REACTOME_GAP_JUNCTION_DEGRADATION", "REACTOME_STRIATED_MUSCLE_CONTRACTION",
#           "BIOCARTA_TNFR2_PATHWAY", "REACTOME_CASPASE_MEDIATED_CLEAVAGE_OF_CYTOSKELETAL_PROTEINS",
#           "BIOCARTA_EGFR_SMRTE_PATHWAY", "BIOCARTA_SPPA_PATHWAY", "REACTOME_SHC1_EVENTS_IN_ERBB4_SIGNALING", "REACTOME_GOLGI_TO_ER_RETROGRADE_TRANSPORT", 
#           "HALLMARK_PROTEIN_SECRETION", "PID_PI3KCI_AKT_PATHWAY", "PID_FOXO_PATHWAY", "PID_TGFBR_PATHWAY", "HALLMARK_TGF_BETA_SIGNALING", "KEGG_RIBOSOME",
#           "KEGG_TIGHT_JUNCTION", "KEGG_ADHERENS_JUNCTION", "HALLMARK_MITOTIC_SPINDLE", "KEGG_FOCAL_ADHESION", "REACTOME_ESR_MEDIATED_SIGNALING",
#           "REACTOME_DNA_REPLICATION", "REACTOME_KILLING_MECHANISMS")

#query= c("BIOCARTA_ECM_PATHWAY", "REACTOME_LAMININ_INTERACTIONS", "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX")

query = c("GO_AUTOPHAGOSOME_ORGANIZATION")

plot_split()



mat <- 
  #data_diff %>% left_join(wilcox_Msig) %>% 
  Fibroblast_processed %>% 
  distinct() %>% 
  filter(Module %in% c(query), CellType == celltype) %>% 
  pivot_longer(c(CF, TM), names_to = "Type") %>% 
  arrange(Type) %>% unite(Subcluster, Subcluster, Type) %>% ungroup( ) %>% 
  select(-Wilcox, -FC, -CellType) %>% 
  pivot_wider(names_from = Module, values_from = value) %>% column_to_rownames("Subcluster")


my.breaks <- c(seq((min(mat)), -0.00001, length.out = 50), seq(0.00001, (max(mat)), length.out = 50))
my.colors <- c(colorRampPalette(colors = c("#002c2d", "#00494b", "#076769", "#208288", "#3ea8a6", "#76c7be", "#bbe4d1", "#ffffe0"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("#ffffe0", "#fed693", "#f5ad52", "#dd8629", "#bd651a", "#9c4511", "#79260b", "#580000"))(length(my.breaks)/2))



#image <- 
pheatmap(t(mat), scale = "row",
         #color = viridis::magma(n = 100, begin = 0, end = 0.9), 
         #color = my.colors,
         #breaks = my.breaks,
         clustering_distance_rows = "correlation", cluster_cols = F, 
         #clustering_distance_cols = "correlation",
         
         fontsize_col = 15,
         fontsize_row = 12, 
         treeheight_row = 10, treeheight_col = 0, 
         #annotation_col = Ann,
         main = "Pathway avg_Score")


save(name = paste0("Msig_Changes_in_", celltype), 1, 12, 9, svg = F)





# Test per CellType
##################################################

names <- names_all[1:3500]
names <- names_all[3500:length(names_all)]

wilcox_Msig_CT2 <- 
  data %>% filter(CellType == "LUM_HR-pos") %>% 
  select(Type, CellType, names) %>% 
  as_tibble() %>% pivot_longer(names, names_to = "Module") %>% 
  group_by(CellType, Module) %>% 
  
  do(w = wilcox.test(value~Type, data = ., paired = F)) %>% 
  summarise(CellType, Module, Wilcox = w$p.value)



sig_modules_CT <- wilcox_Msig_CT %>% filter(Wilcox <= 0.05) %>% pull(Module) %>% unique()

data_avg_CT <- 
  data %>% 
  select(Type, CellType, Subcluster, sig_modules_CT) %>% 
  as_tibble() %>% group_by(CellType, Type) %>% 
  mutate_at(vars(sig_modules_CT), mean) %>% ungroup %>% 
  distinct()


data_diff_CT <- 
  data_avg_CT %>% 
  pivot_longer(sig_modules_CT, names_to = "Module") %>% 
  pivot_wider(names_from = Type, values_from = value) %>% 
  group_by(CellType, Module) %>% 
  mutate(FC = log2(TM/CF))


data_diff_CT %>% left_join(wilcox_Msig_CT) %>% saveRDS(paste0("Output/MsigDB_Scores/gene_sets_GO_scores/CellType_Response_", celltype, ".rds"))
data_diff_CT 



data_diff_CT %>% left_join(wilcox_Msig_CT) %>% select(-Subcluster) %>% distinct() %>% 
  filter(Wilcox <= 0.05, abs(FC) > 0.2, TM > 0.05,
         #str_detect(Module, "REACTOME|BIOCARTA|WP_|PID|KEGG|HALLMARK")
         str_detect(Module, "STEROID")
  ) %>% 
  ungroup() %>% top_n(-40, Wilcox) %>% arrange(-FC)

### Fibroblasts
#query <- c("REACTOME_LAMININ_INTERACTIONS", "BIOCARTA_PTEN_PATHWAY", 
#           "KEGG_SPLICEOSOME", "BIOCARTA_P53HYPOXIA_PATHWAY")#, 
#           "WP_NAD_BIOSYNTHETIC_PATHWAYS", "REACTOME_OTHER_SEMAPHORIN_INTERACTIONS", 
#           "STEGER_ADIPOGENESIS_UP", "PID_INTEGRIN4_PATHWAY", "HENDRICKS_SMARCA4_TARGETS_UP")


#query <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "REACTOME_BINDING_AND_UPTAKE_OF_LIGANDS_BY_SCAVENGER_RECEPTORS", "WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS")


data$Subcluster <- factor(data$Subcluster, levels = c("LAMB1+_Matrix-F", "FBN1+_Matrix-F", "Vasc-F", "Lipo-F", "Chondrocyte"))

#Adipocyte
#query <- c("BIOCARTA_IGF1MTOR_PATHWAY", "BIOCARTA_PTEN_PATHWAY", "REACTOME_NUCLEAR_SIGNALING_BY_ERBB4", "PID_ERBB_NETWORK_PATHWAY", "REACTOME_GLYCOGEN_METABOLISM", "WP_GLYCOGEN_SYNTHESIS_AND_DEGRADATION")

query <- c("GO_REGULATION_OF_STEROID_BIOSYNTHETIC_PROCESS")

#image <- 
data %>% select(Sample, CellType, Subcluster, Type, query) %>% filter(CellType %in% celltype) %>% 
  #scores %>% select(Subcluster, query) %>% 
  pivot_longer(query, names_to = "Module") %>% 
  #group_by(Sample, Subcluster, Module) %>% 
  #mutate(sample_avg = mean(value)) %>% 
  
  ggplot(aes(x = CellType, y = value, fill = Type)) + 
  
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
  
  theme(text = element_text(family = "Lato"),
        axis.text.x = element_text(face = "bold", size = 15, angle = -25, hjust = 0, vjust = 1),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "gray96"),
        strip.text = element_text(size = 13, face = "bold", hjust = 0), 
        legend.key.size = unit(1,"cm"), legend.justification = c(1,1))




save(name = paste0(celltype, " MsigDB Scores"), 0.8, 16, 12, svg = T)












Ann <- data_avg %>% select(Subcluster) %>% 
  filter(Subcluster %notin% c("lun_2", "lun_0", "lun_3", "lun_5", "bas_1")) %>% mutate(N = Subcluster) %>% column_to_rownames("N")


#mat_query <- c("REACTOME_CGMP_EFFECTS", "BIOCARTA_LEPTIN_PATHWAY", "REACTOME_CARNITINE_METABOLISM", "REACTOME_ADENYLATE_CYCLASE_ACTIVATING_PATHWAY",
##  "WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_PUBERTY_STAGE_2_OF_4", "STEIN_ESTROGEN_RESPONSE_NOT_VIA_ESRRA", "REACTOME_MAPK3_ERK1_ACTIVATION",
#  "WP_TRYPTOPHAN_CATABOLISM_LEADING_TO_NAD_PRODUCTION", "REACTOME_KETONE_BODY_METABOLISM", "WP_CHOLESTEROL_BIOSYNTHESIS_PATHWAY",
#  "REACTOME_ATTENUATION_PHASE", "REACTOME_HSF1_ACTIVATION", "KEGG_RIBOSOME", "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
#  "REACTOME_DNA_STRAND_ELONGATION", "WP_DNA_REPLICATION", "TIAN_TNF_SIGNALING_NOT_VIA_NFKB", "WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_INVOLUTION_STAGE_4_OF_4", 
#  "REACTOME_ESTROGEN_DEPENDENT_NUCLEAR_EVENTS_DOWNSTREAM_OF_ESR_MEMBRANE_SIGNALING")


mat_query <- c("REACTOME_SMOOTH_MUSCLE_CONTRACTION", "REACTOME_RHO_GTPASES_ACTIVATE_PAKS", "KEGG_RIBOSOME", "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION", "REACTOME_METALLOTHIONEINS_BIND_METALS",
               "REACTOME_RESPONSE_TO_METAL_IONS")



mat_query <- c("PID_INTEGRIN4_PATHWAY", "REACTOME_LAMININ_INTERACTIONS", "BIOCARTA_ACH_PATHWAY", "BIOCARTA_CTCF_PATHWAY",
               "BIOCARTA_PTEN_PATHWAY", "REACTOME_SIGNALING_BY_ROBO_RECEPTORS", "REACTOME_ATTENUATION_PHASE", "REACTOME_SEMA3A_PAK_DEPENDENT_AXON_REPULSION")


#mat <- data_avg %>% select(Subcluster, mat_query) %>% 
#  #filter(Subcluster %notin% c("lun_2", "lun_0", "lun_3", "lun_5", "bas_1")) %>% 
#  column_to_rownames("Subcluster")



mat_query <- c(up, down)


mat <- 
  Fibroblast_processed %>% filter(Module %in% mat_query) %>% pivot_longer(c(CF, TM), names_to = "Type")  %>% arrange(Type) %>% unite(Subcluster, Subcluster, Type) %>% select(-Wilcox, -FC) %>% pivot_wider(names_from = Module, values_from = value) %>% column_to_rownames("Subcluster")









mods1 <- data_diff %>% left_join(wilcox) %>% 
  filter(Wilcox <= 0.05, FC < -0.2, 
         str_detect(Module, "REACTOME|BIOCARTA|WP_|PID|KEGG|HALLMARK|NFIC|SMARCA4")
  ) %>% filter(CF > 0) %>% group_by(Subcluster) %>% top_n(-8, Wilcox) %>% pull(Module) %>% unique()

mods2 <- data_diff %>% left_join(wilcox) %>% 
  filter(Wilcox <= 0.05, FC > 0.2, 
         str_detect(Module, "REACTOME|BIOCARTA|WP_|PID|KEGG|HALLMARK|NFIC|SMARCA4")
  ) %>% filter(TM > 0) %>% group_by(Subcluster) %>% top_n(-8, Wilcox) %>% pull(Module) %>% unique()



mat <- data_diff %>% 
  filter(Module %in% c(mods1, mods2)) %>% 
  select(Subcluster, Module, CF, TM) %>% 
  pivot_wider(names_from = Subcluster, values_from = FC) %>% 
  column_to_rownames("Module")


mat <- data_diff %>% 
  filter(Module %in% c(mods1, mods2)) %>% 
  select(Subcluster, Module, CF, TM) %>% 
  pivot_longer(c(CF, TM)) %>% unite(Subcluster, Subcluster, name) %>% 
  pivot_wider(names_from = Subcluster, values_from = value) %>% 
  column_to_rownames("Module")


#mat <- data_diff %>% left_join(wilcox) %>% 
#  filter(Wilcox <= 0.05, abs(FC) > 0.2) %>% 
#  select(Subcluster, Module, FC) %>% 
#  pivot_wider(names_from = Subcluster, values_from = FC) %>% 
#  column_to_rownames("Module")


mat[is.na(mat)] <- 0


#image <- 
pheatmap(mat, 
         scale = "row", 
         fontsize = 10, cluster_cols = F,
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





query <- "REACTOME_SEMA3A_PAK_DEPENDENT_AXON_REPULSION"


add <- data %>% select(X1, query) %>% column_to_rownames("X1")


Sobj_integrated <- AddMetaData(Sobj_integrated, add)


VlnPlot(Sobj_integrated, 
        idents = c("Lipo-F", "LAMB1+_Matrix-F", "FBN1+_Matrix-F", "Vasc-F"),
        features = query,
        #features = c("HENDRICKS_SMARCA4_TARGETS_UP"),
        group.by = "Subcluster", 
        split.by = "Type",
        pt.size = 0, 
        #y.max = 0.5, 
        #ncol = 5,  
        cols = viridis::magma(n = 2, begin = 0.3, end = 0.8)
)





names <- names(gmt)
Marker_list    <- vector("list", length = length(names))

for (i in 1:length(names)) {
  
  Marker_list[[i]]  <- length(gmt[[names[i]]]) %>% as_tibble() %>% add_column(Module = names[i])
  
}

Bind <- do.call(rbind.data.frame, Marker_list)


