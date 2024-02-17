panel_4d <- function() {

query = c("ENAH", "ITGA2", "ITGB8")

Sample_averages %>% filter(Gene %in% query, Cluster %in% c("LUM_HR-neg")) %>% 
  
  ggplot(aes(x = factor(Cluster, levels = levels_samp), y = avg_Expr, fill = Type)) +
  
  geom_boxplot(color = "grey30", 
               width = 0.5, 
               alpha = 1, 
               position=position_dodge(0.9), 
               lwd = 1.5,
               outlier.alpha = 0) +
  
  geom_point(aes(color = factor(Sample, levels = levels_samp)),
             position = position_jitterdodge(jitter.width = 0.000005, 
                                             jitter.height = 0,
                                             dodge.width = 0.5), 
             alpha = 1,
             size  = 4) +
  
  
  scale_fill_manual(values = c("#A6499B", "#FAA42F")) +
  scale_color_manual(values = cols_samp) +
  
  theme(text = element_text(family = "Lato", size = 65),
        strip.background = element_blank(),
        axis.ticks = element_line(size = 1, colour = "grey30", lineend = "round"),
        axis.text.x = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey30", size = 2),
        panel.grid = element_blank()
        
  ) +
  
  facet_wrap(vars(Gene), scales = "free_y", ncol = 3)

}

p <- panel_4d()

p %>% save_x(data = ., name = paste0("Panel_4d_ENAH_ITG"), 1, 12, 5, svg = T)



Modules_List <- readRDS("utilities/Modules_List.rds")
Sobj_integrated <- AddModuleScore  (Sobj_integrated, list((top_n(Modules_List$LUMA, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "LUM")
Sobj_integrated <- AddModuleScore  (Sobj_integrated, list((top_n(Modules_List$LUPR, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "LUP")
Sobj_integrated <- AddModuleScore  (Sobj_integrated, list((top_n(Modules_List$MASC, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "MAS")
Sobj_integrated <- AddModuleScore  (Sobj_integrated, list((top_n(Modules_List$STRM, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "STR")


panel_S4a <- function(size, nrow) {

nrow = nrow
size = size
  
margin = 0
look = theme(panel.background = element_rect(fill  = "white", 
                                             color = "grey60", 
                                             size  = 0.25), 
             panel.grid = element_blank(), plot.margin = unit(c(margin,margin,margin,margin), "cm"),
             axis.line = element_blank(),
             axis.title = element_blank(), 
             axis.text  = element_blank(),
             axis.ticks = element_blank(),
             legend.position = "none")



Sobj_integrated <- AddMetaData(Sobj_integrated, metadata = as.data.frame(Sobj_integrated@reductions$umap@cell.embeddings))
df_rna = Sobj_integrated@meta.data

set.seed(467)
rows_rna  <- sample(nrow(df_rna))


p <- 
  df_rna[rows_rna, ] %>% 
  
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = CellType)) + 
  geom_point(alpha = 0.35, 
             size = c(size+1), 
             stroke = 0) +
  
  scale_color_manual(values = c("#44aa99", "#88ccee", "#117733")) + 
  look


p1 <- 
FeaturePlot(Sobj_integrated, reduction = "umap",
            features = c("LUM1"),
            pt.size = size, 
            order = T, 
            cols = viridis::magma(n = 100)
) + look

p2 <- 
FeaturePlot(Sobj_integrated, reduction = "umap",
            features = c("LUP1"),
            pt.size = size, 
            order = T, 
            cols = viridis::magma(n = 100)
) + look

p3 <- 
FeaturePlot(Sobj_integrated, reduction = "umap",
            features = c("MAS1"),
            pt.size = size, 
            order = T, 
            cols = viridis::magma(n = 100)
) + look


p + p1 + p2 + p3 + patchwork::plot_layout(nrow = nrow)

}

p <- panel_S4a(size = 1, nrow = 1)

p %>% save_x(data = ., name = paste0("Epithelial_Cell_Calling"), 1.5, 12, 3, svg = F)


panel_4a <- function(){

margin = 0.5

look = theme(panel.background = element_rect(fill  = "white", 
                                             color = "grey60", 
                                             size  = 1), 
             panel.grid = element_blank(), plot.margin = unit(c(margin,margin,margin,margin), "cm"),
             axis.title = element_blank(), 
             axis.text  = element_blank(),
             axis.ticks = element_blank(),
             legend.position = "none")

celltype = "LUM_HR-neg"

Sobj_integrated <- AddMetaData(Sobj_integrated, metadata = as.data.frame(Sobj_integrated@reductions$umap@cell.embeddings))
df_rna = Sobj_integrated@meta.data

set.seed(467)
rows_rna  <- sample(nrow(df_rna)) #, 10000)


P1 <- 
  df_rna[rows_rna, ] %>% #filter(CellType == "Myeloid") %>% 
  
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = Subcluster)) + 
  geom_point(alpha = 0.35, 
             size = 3, 
             stroke = 0) +
  
  #facet_wrap(~ CellType, nrow = 1) +
  
  scale_color_manual(values = pal_d3("category20")(20)) + 
  #scale_color_manual(values = colors[1:4]) + 
  #scale_color_manual(values = colors[7:10]) + 
  
  look + 
  theme(strip.background = element_blank(), panel.background = element_blank(),
        strip.text = element_blank())


P2 <- 
  df_rna[rows_rna, ] %>% #filter(CellType == "Myeloid") %>% 
  
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = factor(Sample, levels = levels_samp))) + 
  
  geom_point(data = filter(select(df_rna[rows_rna, ], -Type)), color = "grey88",
             alpha = 1, 
             size = 3, 
             stroke = 1) +
  
  geom_point(alpha = 1, 
             size = 1.5, 
             stroke = 0) +
  
  facet_wrap(~Type, ncol = 2) +
  
  scale_color_manual(values = cols) + 
  
  look + 
  theme(strip.background = element_blank(), panel.background = element_blank(),
        strip.text = element_blank())


P1 + as.ggplot(P2) + patchwork::plot_layout(ncol = 1, heights = c(2,1)) 

}

p <- panel_4a()

p %>% save_x(data = ., name = paste0("LUM_HR-neg_Subcluster_Umap"), 1.2, 8, 10, svg = F)



panel_4b <- function(){

query <- Bind %>% ungroup() %>% 
  filter(str_detect(Module, "REACTOME|BIOCARTA|WP_|PID|KEGG|HALLMARK"), group1 > 0.2) %>% 
  group_by(CLuster) %>% top_n(2, FC) %>% arrange(CLuster) %>%   
  pull(Module)


mat <- 
  data %>% 
  select(Subcluster, query) %>% 
  pivot_longer(query) %>% 
  group_by(Subcluster, name) %>% 
  summarise(mean = mean(value)) %>% 
  pivot_wider(names_from = Subcluster, values_from = mean) %>% 
  mutate(name = tolower(name)) %>% 
  column_to_rownames("name")


pheatmap(mat[tolower(query),], scale = "row",
         color = viridis::magma(n = 100, begin = 0, end = 0.9), 
         #color = my.colors,
         #breaks = my.breaks,
         clustering_distance_rows = "correlation", 
         cluster_cols = F, 
         cluster_rows = F,
         fontsize_col = 15,
         fontsize_row = 12, 
         treeheight_row = 10, treeheight_col = 0, 
         #annotation_col = Ann,
         main = "Pathway avg_Score")
}

p <- panel_4b()

p %>% save_x(data = ., name = paste0("LUM_HR-neg_Marker_Pathways"), 1, 16, 10, svg = T)


### might be needed for SUP 4 Table
data_diff_CT %>% ungroup() %>% 
  select(-Subcluster) %>% 
  distinct() %>% 
  filter(str_detect(Module, "KEGG"), CF > 0.1, Wilcox <= 0.05) %>% 
  arrange(Wilcox) %>% arrange(FC) %>%  
  print(n = 25) 



BasalChromvar <- function(){

look <-     theme(text = element_text(family = "Lato"),
                  #axis.text.x = element_text(face = "bold", size = 15, angle = -25, hjust = 0, vjust = 1),
                  #axis.title.x = element_blank(),
                  axis.title.y = element_text(face = "bold", size = 18),
                  axis.text.y  = element_text(size = 15), legend.position = "none",
                  
                  axis.title.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  
                  panel.grid = element_blank(),
                  
                  
                  panel.background = element_rect(fill = "white", color = "grey50", size = 3),
                  
                  strip.text = element_text(size = 13, face = "bold", hjust = 0), 
                  legend.key.size = unit(1,"cm"), legend.justification = c(1,1))




motifZscoreMat %>% separate(Var2, into = c("Type", "B"), sep = "-") %>% 
  filter(Var1 %in% c("TEAD1_796", "NFIC_740")) %>%
  
  ggplot(aes(x = Subcluster, y = value, fill = Type)) + 
  
  geom_violin(aes(fill = Type),
              scale = "width", 
              alpha = 1, 
              trim = F,
              lwd = 1.25, 
              color = "grey30") + 
  
  geom_boxplot(aes(color = Type), 
               fill = "grey30", 
               width = 0.1, alpha = 1, 
               lwd = 1.25,
               position=position_dodge(0.9), outlier.alpha = 0) +
  
  stat_summary(fun = "median", 
               geom = "point", size = 1.5,
               position = position_dodge(0.9),
               color = "floralwhite") +
  
  
  scale_fill_manual(values = c("#9d3396", "#eb933e")) +
  scale_color_manual(values = c("grey30", "grey30")) +
  
  facet_wrap(facets = "Var1", scales = "free_y") +
  
  look

}

p2 <- BasalChromvar()

p %>% save_x(data = ., name = paste0("Basal_Chromvar_Scores_by_Subcluster"), 1, 16, 10, svg = T)




p <- p1 + p2 + patchwork::plot_layout(ncol = 1)

p %>% save_x(data = ., name = paste0("Basal_Chromvar_Scores_by_Subcluster"), 1, 16, 10, svg = T)






motifZscoreMat %>% separate(Var2, into = c("Type", "B"), sep = "-", extra = "merge") %>% 
  filter(Var1 %in% c("TEAD1_796", "NFIC_740")) %>% #head() %>% 
  group_by(Subcluster, Var1) %>%  
  do(w = wilcox.test(value~Type, data = ., paired = F)) %>% 
  summarise(Subcluster, Var1, Wilcox = w$p.value)
