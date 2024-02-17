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


save     <- function(name, scale, w, h, svg) {
  
  svg = svg
  
  png.path <- paste0("Submission/Figures/PNG")
  svg.path <- paste0("Submission/Figures/SVG")
  
  
  if (svg == TRUE) {
    ggsave(plot     = image,
           path     = svg.path,
           filename = paste0(name,".svg"),
           width    = w, 
           height   = h,
           scale    = scale,
           dpi      = 300)
  }
  
  
  ggsave(plot     = image,
         path     = png.path,
         filename = paste0(name,".png"),
         width    = w, 
         height   = h,
         scale    = scale,
         dpi      = 300)
  
}
save_x     <- function(data, name, scale, w, h, svg) {
  
  svg = svg
  
  png.path <- paste0("Submission/Figures/PNG")
  svg.path <- paste0("Submission/Figures/SVG")
  
  
  if (svg == TRUE) {
    ggsave(plot     = data,
           path     = svg.path,
           filename = paste0(name,".svg"),
           width    = w, 
           height   = h,
           scale    = scale,
           dpi      = 300)
  }
  
  
  ggsave(plot     = data,
         path     = png.path,
         filename = paste0(name,".png"),
         width    = w, 
         height   = h,
         scale    = scale,
         dpi      = 300)
  
}

cols_type <- c("#9d3396", "#eb933e")
cols_cell   <- c("#44aa99", "#88ccee", "#117733", "#332288", "#ddcc77", "#882255", "#aa4499", "#cc6677", "#999933", "#dddddd")
levels <- c("LUM_HR-pos", "LUM_HR-neg", "Basal", "Fibroblast", "Adipocyte", "Blood_EC", "Lymphatic_EC", "Vasc_Accessory", "Myeloid", "Lymphoid")


data <- Sobj_integrated@meta.data

plot_vln_UMI <- function() {

p1 <- data %>% 
ggplot(aes(x = Type, y = nCount_RNA)) + 
  
  geom_violin(aes(fill = Type),
              scale = "width", 
              alpha = 0.75, 
              trim = F,
              lwd = 1.25, 
              color = "grey30") + 
  
  geom_boxplot(aes(color = Type), 
               fill = "grey30", 
               width = 0.1, alpha = 1, 
               lwd = 1.25,
               position=position_dodge(0.7), outlier.alpha = 0.8) +
  
  stat_summary(fun = "median", 
               geom = "point", size = 1.5,
               position = position_dodge(0.7),
               color = "floralwhite") +
  
  
  scale_fill_viridis_d(option = "C", begin = 0.4, end = 0.8) +
  scale_color_manual(values = c("grey30", "grey30")) +
  
  
  scale_y_log10() +
  
  theme(text = element_text(family = "Lato"),
        #axis.text.x = element_text(face = "bold", size = 15, angle = -25, hjust = 0, vjust = 1),
        #axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 18),
        axis.text.y = element_text(size = 15), legend.position = "none",
        
        axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        
        #panel.background = element_rect(fill = "gray96"),
        panel.background = element_blank(),
        strip.text = element_text(size = 13, face = "bold", hjust = 0), 
        legend.key.size = unit(1,"cm"), legend.justification = c(1,1))



#function() {
#p2 <- 
  data %>% 
  ggplot(aes(x = CellType, y = nCount_RNA, fill = Type)) + 
  
  geom_split_violin(scale = "width", 
                    alpha = 0.75, 
                    trim = F,
                    lwd = 1.25, 
                    color = "grey30") + 
  
  geom_boxplot(aes(color = Type), 
               fill = "grey30", 
               width = 0.1, alpha = 1, 
               lwd = 1.25,
               position=position_dodge(0.7), outlier.alpha = 0.8) +
  
  stat_summary(fun = "median", 
               geom = "point", size = 1.5,
               position = position_dodge(0.7),
               color = "floralwhite") +
  
  
  scale_fill_viridis_d(option = "C", begin = 0.4, end = 0.8) +
  scale_color_manual(values = c("grey30", "grey30")) +
  
  
  scale_y_log10() +
  
  theme(text = element_text(family = "Lato"),
        #axis.text.x = element_text(face = "bold", size = 15, angle = -25, hjust = 0, vjust = 1),
        
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        
        axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        
        panel.background = element_rect(fill = "gray96"),
        strip.text = element_text(size = 13, face = "bold", hjust = 0), 
        legend.key.size = unit(1,"cm"), legend.justification = c(1,1))

  
 p2 <-  
  data %>% 
    ggplot(aes(x = CellType, y = nCount_RNA, fill = Type)) + 
    
    #geom_split_violin(scale = "width", 
    #                  alpha = 0.75, 
    #                  trim = F,
    #                  lwd = 1.25, 
    #                  color = "grey30") + 
    
    geom_boxplot(aes(color = Type), 
                 #fill = "grey30", 
                 width = 0.5, alpha = 0.75, 
                 lwd = 1.25,
                 position=position_dodge(0.7), outlier.alpha = 0.5) +
    
    stat_summary(fun = "median", 
                 geom = "point", size = 1.5,
                 position = position_dodge(0.7),
                 color = "floralwhite") +
    
    
    scale_fill_viridis_d(option = "C", begin = 0.4, end = 0.8) +
    scale_color_manual(values = c("grey30", "grey30")) +
    
    
    scale_y_log10() +
    
    theme(text = element_text(family = "Lato"),
          #axis.text.x = element_text(face = "bold", size = 15, angle = -25, hjust = 0, vjust = 1),
          
          axis.title.y = element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = -25, hjust = 0), 
          #axis.ticks.x = element_blank(),
          
          #panel.background = element_rect(fill = "gray96"),
          panel.background = element_blank(),
          strip.text = element_text(size = 13, face = "bold", hjust = 0), 
          legend.key.size = unit(1,"cm"), legend.justification = c(1,1))  






p1 + p2 + patchwork::plot_layout(ncol = 2, widths = c(1,5))

}


#p <- 
plot_vln_UMI()

p %>% save_x(data = ., name = paste0("UMI_Comparison"), 1, 10, 5, svg = T)








cols <- colorRampPalette( colors = brewer.pal(11,"BrBG") )



mat <- 
df %>% group_by(Subcluster) %>% 
  mutate(Z = (Median - mean(Median)/sd(Median))) %>% 
  select(Motif, Subcluster, Z) %>% 
  pivot_wider(names_from = Subcluster, values_from = Z) %>% 
  column_to_rownames("Motif")

pheatmap(mat, cluster_rows = F)


mat <- df %>% 
  select(Motif, Subcluster, MedianZ) %>% 
  pivot_wider(names_from = Subcluster, values_from = MedianZ) %>% 
  column_to_rownames("Motif")

pheatmap(mat, cluster_rows = F)


mat <- df %>% 
  select(Motif, Subcluster, Median) %>% 
  pivot_wider(names_from = Subcluster, values_from = Median) %>% 
  column_to_rownames("Motif")

pheatmap(mat, cluster_rows = F, scale = "row", color = cols(100))


p %>% save_x(data = ., name = paste0("Basal_Subscluster_Motif_Scores"), 1, 10, 5, svg = T)









Bind <- readRDS("Submission/utilities/Reactome_Splicing_GeneSet_FC.rds")

mat <- 
  Bind %>% 
  select(Gene, avg_logFC, Cluster) %>% 
  pivot_wider(names_from = Cluster, values_from = avg_logFC) %>% 
  column_to_rownames("Gene")


mat[is.na(mat)] <- 0


my.breaks <- c(seq((min(mat)), -0.00001, length.out = 50), seq(0.00001, (max(mat)), length.out = 50))
my.colors <- c(colorRampPalette(colors = c("#002c2d", "#00494b", "#076769", "#208288", "#3ea8a6", "#76c7be", "#bbe4d1", "#ffffe0"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("#ffffe0", "#fed693", "#f5ad52", "#dd8629", "#bd651a", "#9c4511", "#79260b", "#580000"))(length(my.breaks)/2))

#my.colors <- c(colorRampPalette(colors = c("#10194d", "#163771", "#1c5796", "#397aa8", "#579eb9", "#89c0c4", "#bce2cf", "#ffffe0"))(length(my.breaks)/2), 
#               colorRampPalette(colors = c("#ffffe0", "#fad4ac", "#f0a882", "#e47961", "#c65154", "#a52747", "#751232", "#4a001e"))(length(my.breaks)/2))


#image <- 
pheatmap(t(mat), 
         color = my.colors,
         breaks = my.breaks, 
         treeheight_col = 15,
         treeheight_row = 15,
         border_color = NA, 
         fontsize_col = 8,
         fontsize_row = 20, 
         #main = "logFC Reactome_mRNA Splicing")
         main = "logFC Reactome: Eukaryotic Translation Elongation")


save(name = "Reactome: Eukaryotic Translation Elongation_logFC", 1, 16, 6, svg = F)     





# make a random gene set
genes <- rownames(Sobj_integrated) %>% 
  as.data.frame() %>% filter(!str_detect(., "\\.1"), !str_detect(., "\\.2"), !str_detect(., "\\.3"), !str_detect(., "LINC0"), !str_detect(., "-AS")) %>% 
  sample_n(1000) %>% pull(.) %>% as.character()


Bind <- readRDS("Submission/utilities/Random_GeneSet_FC.rds")


mat <- 
  Bind %>% 
  select(Gene, avg_logFC, Cluster) %>% 
  pivot_wider(names_from = Cluster, values_from = avg_logFC) %>% 
  column_to_rownames("Gene")


mat[is.na(mat)] <- 0


my.breaks <- c(seq((min(mat)), -0.00001, length.out = 50), seq(0.00001, (max(mat)), length.out = 50))
my.colors <- c(colorRampPalette(colors = c("#002c2d", "#00494b", "#076769", "#208288", "#3ea8a6", "#76c7be", "#bbe4d1", "#ffffe0"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("#ffffe0", "#fed693", "#f5ad52", "#dd8629", "#bd651a", "#9c4511", "#79260b", "#580000"))(length(my.breaks)/2))


mat_s <- sample_n(mat, 100)

#image <- 
pheatmap(t(mat), 
         color = my.colors,
         breaks = my.breaks, 
         treeheight_col = 15,
         treeheight_row = 15,
         border_color = NA, 
         fontsize_col = 0.000001,
         fontsize_row = 20, 
         main = "logFC Random GeneSet (n = 100)") #color = colorRampPalette(c("Red", "Green", "Blue"), alpha = TRUE)(100), breaks = c(-1, 0, 0.1, 2))


save(name = "Random_GeneSet_logFC", 1, 9, 6, svg = F)     







# HR Receptor average expression
mat <- cluster.averages %>% 
  filter(Gene %in% c("AR", "ESR1", "PGR")) %>% 
  select(Gene, CF, TM, Cluster) %>% 
  pivot_longer(c(CF, TM), names_to = "X") %>% 
  unite(Cluster, Cluster, X) %>% 
  pivot_wider(names_from = Gene) %>% 
  column_to_rownames("Cluster")




pheatmap(mat, 
         scale = "column", 
         color = viridis::magma(n = 100))




query = c("ESR1", "AR", "PGR")

image <- 
Sample_averages %>% filter(Gene %in% query, 
                           #Cluster %in% celltype, 
                           ) %>% 
  ggplot(aes(x = Cluster, y = avg_Expr, fill = Type)) +
  
  #geom_violin(alpha = 0.65, trim = F, scale = "width" , draw_quantiles = T) +
  
  geom_boxplot(color = "grey30", 
               width = 0.75, 
               alpha = 0.8, 
               position=position_dodge(0.9), 
               lwd = 1,
               outlier.alpha = 0) +
  
  geom_point(aes(color = Type), color = "grey30",
             position = position_jitterdodge(jitter.width = 0.05, 
                                             jitter.height = 0,
                                             dodge.width = 0.9), 
             alpha = 0.65,
             size = 2) +
  
  stat_compare_means(paired = F, method = "wilcox", 
                     label.x.npc = 0.25, vjust = -5) +
  
  scale_fill_viridis_d(option = "C", begin = 0.3, end = 0.8) +
  scale_color_viridis_d(begin = 0.3, end = 0.6) +
  #scale_color_manual(values =  c("gray30", "gray30")) +
  
  ggtitle(label = "Sample Average Expression") +
  
  theme(text = element_text(family = "Lato"),
        axis.text = element_text(size=10), 
        #axis.title=element_text(size=14,face="bold"),
        axis.title = element_blank(), 
        axis.text.y = element_text(face = "bold", size = 12),
        legend.position = c(0.95, 0.9), 
        legend.key.size = unit(1, "cm"),
        panel.background = element_rect(fill = "grey95"),
        strip.text  = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(size = 15, angle = -25, hjust = 0, vjust = 1, face = "bold")
  ) +
  #scale_y_sqrt() +
  #scale_y_log10() +
  facet_wrap(vars(Gene), scales = "free_y", ncol = 1)



save(name = "Hormone_Receptor_Expression", 1, 9, 8, svg = F)     







# TF usage in LUM_HR-pos
colnames <- Sobj_integrated@meta.data %>% select(contains("_1")) %>% colnames()
colnames_new <- colnames %>% str_replace(pattern = "_1", replacement = "_grn")
Sobj_integrated@meta.data <- Sobj_integrated@meta.data %>% rename_at(vars(colnames), ~ colnames_new)


features = c("AR_grn", "ESR1_grn", "PGR_grn", "CUX2_grn", "ZNF689_grn", "BATF_grn")
             #, "JUN_grn"  , "NR4A1_grn", ,
             #"CERS4_grn", "ID2_grn", "TBX3_grn", "RPS4X_grn")

plot_list    <- vector("list", length = length(features))


for (i in 1:length(features)) {
  
plot_list[[i]] <- 
  
  FeaturePlot(Sobj_integrated, reduction = "umap",
            #cells = WhichCells(Sobj_integrated, expression = SCORE_1 >= 0.25),
            #features = c("PLXNA4", "HPSE2", "SEMA3C", "EPHA3", "SLC22A3", "ZBTB7C"), 
            features = features[i],
            #split.by = "Type",
            pt.size = 0.25, 
            order = F, 
            #max.cutoff = 0.45,
            cols = viridis::magma(n = 100),
            
            ) + theme(axis.title = element_blank(), 
                      axis.text = element_blank(),
                      axis.ticks = element_blank(), 
                      legend.position = "none"
                      )


}

p <- cowplot::plot_grid(plotlist = plot_list)


p %>% save_x(name = "LUM_Pos_TF-GRN_Modules", 1.5, 10, 6, svg = F)    







# gprofilter of all DE genes total

library(gprofiler2)

Marks  <-  readRDS("Submission/utilities/DE_Test_Everything_TMvCF_20kDownSample.rds")

query        <- filter(Marks, avg_logFC > 0) %>% rownames()
Results      <- gost(query, organism = "hsapiens", evcodes = T)
query        <- filter(Marks, avg_logFC < -0) %>% rownames()
Result_Down <- gost(query, organism = "hsapiens", evcodes = T)


UP   <- head(filter(Results$result, source == "REAC"), 5) %>% mutate(Dir = "UP")
DOWN <- head(filter(Result_Down$result, source == "REAC"), 5) %>% mutate(Dir = "DOWN")

df <- bind_rows(UP, DOWN)
order <- df$term_name
df$term_name <- factor(df$term_name, levels = order)


plot <- 
df %>% ggplot(aes(y = term_name, x = intersection_size, fill = p_value)) + 
  geom_col() + 
  scale_y_discrete(limits = rev(levels(df$term_name))) + 
  
  #scale_fill_manual(values = cols) +
  scale_fill_viridis_c(direction = -1, begin = 0, end = 0.75) +
  
  theme(text = element_text(family = "Lato"), 
        axis.title = element_blank(), 
        axis.text.y = element_text(size = 12)) 

plot %>% save_x(data = ., name = "Global_Pathway_Enrichment", 1, 12, 5, svg = T)    
  





# GTEX TPM of LUM-pos response genes

dflong <- readRDS("Output/Reboot/GTEX_MeadianTPM_long.rds")

function() {

lumpos_u <- Treatment_Response %>% filter(Cluster == "LUM_HR-pos") %>% top_n( 30, avg_logFC) %>% pull(Gene)
lumpos_d <- Treatment_Response %>% filter(Cluster == "LUM_HR-pos") %>% top_n(-15, avg_logFC) %>% pull(Gene)
top40 <- c(lumpos_u, lumpos_d)

lumpos <- Treatment_Response %>% filter(Cluster == "LUM_HR-pos", abs(avg_logFC) > 0.5) %>% pull(Gene)

Ann <- Treatment_Response %>% mutate(AR_Response = ifelse(avg_logFC > 0, "UP", "DOWN")) %>% filter(Cluster == "LUM_HR-pos", Gene %in% lumpos) %>% select(Gene, AR_Response) %>% column_to_rownames("Gene")
Ann2 <- dflong %>% filter(hgnc_symbol %in% c("AR", "ESR1", "PGR")) %>% select(-1, -2) %>% pivot_wider(names_from = hgnc_symbol, values_from = value) %>% column_to_rownames("Tissue")

ann_colors = list(AR_Response = c(UP = "Orange", DOWN = "Purple"))
exclude = c("ENSG00000196873.15", "ENSG00000260807.6")

tissues <- dflong %>% 
  filter(hgnc_symbol %in% c("AR", "ESR1", "PGR")) %>% 
  group_by(hgnc_symbol) %>% 
  top_n(8, value) %>% 
  print(n = 50) %>% pull(Tissue) %>% unique()

tissues <- 
  str_remove(tissues, "Cervix...Ectocervix") %>% str_remove("Artery...Tibial")

#tissues <- c("Testis", "Breast...Mammary.Tissue", "Liver", "Fallopian.Tube", "Ovary", "Prostate", "Uterus", "Vagina", "Artery...Aorta")


mat <- 
dflong %>% 
  filter(hgnc_symbol %in% lumpos, ID %notin% exclude, Tissue %in% tissues) %>% 
  select(-1, -2) %>% 
  pivot_wider(names_from = Tissue, values_from = value) %>% 
  column_to_rownames("hgnc_symbol")


my.breaks <- c(seq((min(mat)), -0.00001, length.out = 50), seq(0.00001, (max(mat)), length.out = 50))
my.colors <- c(colorRampPalette(colors = c("#002c2d", "#00494b", "#076769", "#208288", "#3ea8a6", "#76c7be", "#bbe4d1", "#ffffe0"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("#ffffe0", "#fed693", "#f5ad52", "#dd8629", "#bd651a", "#9c4511", "#79260b", "#580000"))(length(my.breaks)/2))

test <- mat %>% rownames_to_column("Gene") %>% select(1) %>% mutate(X = ifelse(Gene %in% top40, Gene, "")) %>% pull(X)

p <- 
pheatmap(t(mat), scale = "column", border_color = NA, labels_col = test, 
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation", 
         clustering_method = "ward.D2",
         #color = magma(n = 100), 
         color = my.colors, treeheight_col = 50,
         annotation_col = Ann, fontsize_row = 15, 
         annotation_row = Ann2,
         annotation_colors = ann_colors, #angle_col = 45,
         main = "GTEX median TPM") #%>% 
  
  p %>% save_x(data = ., name = "LUM_Pos_GTEX_Comparison", 1, 12, 5, svg = T)    

}








lumpos <- Treatment_Response %>% filter(Cluster == "LUM_HR-pos", abs(avg_logFC) >= 0.35) %>% pull(Gene)

mat <- dflong %>% filter(hgnc_symbol %in% lumpos) %>% select(-1, -2) %>% dplyr::rename(Gene = "hgnc_symbol") %>% filter(Gene %notin% exclude) %>% pivot_wider(names_from = Tissue, values_from = value) %>% column_to_rownames("Gene")

Ann <- Treatment_Response %>% mutate(AR_Response = ifelse(avg_logFC > 0, "UP", "DOWN")) %>% filter(Cluster == "LUM_HR-pos", Gene %in% lumpos) %>% select(Gene, AR_Response) %>% column_to_rownames("Gene")



my.breaks <- c(seq((min(mat)), -0.00001, length.out = 50), seq(0.00001, (max(mat)), length.out = 50))
my.colors <- c(colorRampPalette(colors = c("#002c2d", "#00494b", "#076769", "#208288", "#3ea8a6", "#76c7be", "#bbe4d1", "#ffffe0"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("#ffffe0", "#fed693", "#f5ad52", "#dd8629", "#bd651a", "#9c4511", "#79260b", "#580000"))(length(my.breaks)/2))


image <- 
pheatmap(t(mat[,tissues]), scale = "column", 
         border_color = NA, cutree_cols = 4,
         
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation", clustering_method = "ward.D2", 
         show_colnames = F,   
         color = my.colors, treeheight_col = 5,
         annotation_col = Ann, fontsize_row = 15, 
         #annotation_row = Ann2,
         annotation_colors = ann_colors, 
         angle_col = 45,
         main = "GTEX median TPM")

cutree(image$tree_col, k = 4) %>% as.data.frame() %>% rename(Clust = ".") %>% filter(Clust == 1)


Sobj_integrated <- AddMetaData(Sobj_integrated, metadata = as.data.frame(Sobj_integrated@reductions$umap@cell.embeddings))
Sobj_integrated$CellType <- factor(Sobj_integrated$CellType, levels = levels)

plot_global_umaps <- function(alpha, size, stroke) {
  
  p1 <- 
    df_rna[rows_rna, ] %>% 
    
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = CellType)) + 
    geom_point(alpha = alpha, 
               size = size, 
               stroke = stroke) +
    
    scale_color_manual(values = cols_cell) +
    look
  
  p2 <- 
    df_rna[rows_rna, ] %>% 
    
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = Type)) + 
    geom_point(alpha = alpha, 
               size = size, 
               stroke = stroke) +
    
    scale_color_manual(values = cols_type) +
    look
  
  
  
  p3 <- 
    df_ <- atac[rows_atac, ] %>% 
    
    ggplot(aes(x = UMAP1, y = UMAP2, color = predictedCellTypeGroup)) + 
    geom_point(alpha = alpha, 
               size = size, 
               stroke = stroke) +
    
    scale_color_manual(values = cols_cell) +
    look
  
  
  p4 <- 
    df_atac[rows_atac, ] %>% 
    
    ggplot(aes(x = UMAP1, y = UMAP2, color = SampleType)) + 
    geom_point(alpha = alpha, 
               size = size, 
               stroke = stroke) +
    
    scale_color_manual(values = cols_type) +
    look
  
  
  p1 + p2 + p3 + p4
  
}



df_rna = Sobj_integrated@meta.data
df_atac = AtacSeqEmbeddingBasedOnPeakMatrix_20201007

set.seed(42)
rows_rna  <- sample(nrow(df_rna)) #, 10000)
rows_atac <- sample(nrow(df_atac))#, 10000)






###### Make pie charts with proportions

df_rna = Sobj_integrated@meta.data
df_atac = AtacSeqEmbeddingBasedOnPeakMatrix_20201007


FC_rna <- 
  as.data.frame(prop.table(table(df_rna$CellType, 
                                 df_rna$Type), 
                           margin = 2)) %>%
  mutate(Freq = Freq * 100) %>% 
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  mutate(FC = log(TM/CF)) %>% select(-CF, -TM)


FC_atac <- 
  as.data.frame(prop.table(table(df_atac$predictedCellTypeGroup, 
                                 df_atac$SampleType), 
                           margin = 2)) %>%
  mutate(Freq = Freq * 100) %>% 
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  rename(CF = 2, TM = 3) %>% 
  mutate(FC = log(TM/CF)) %>% select(-CF, -TM)



rna_ready <- 
  as.data.frame(prop.table(table(df_rna$CellType))) %>% 
  left_join(FC_rna) %>% mutate(data = "RNA") 


all_ready <- 
  as.data.frame(prop.table(table(df_atac$predictedCellTypeGroup))) %>% 
  left_join(FC_atac) %>% mutate(data = "ATAC") %>% bind_rows(rna_ready)


all_ready$Var1 <- factor(all_ready$Var1, levels = levels)



plot_props <- function(){
  
  p1 <- 
    all_ready %>% 
    
    ggplot(aes(x = "", y = Freq, fill = Var1)) + 
    geom_bar(stat = "identity", position = position_fill()) +
    #geom_text(aes(label = Freq), position = position_fill(vjust = 0.5)) +
    
    #geom_bar(aes(fill = FC), stat = "identity", position = position_fill(), alpha = 0.5) +
    
    coord_polar(theta = "y") +
    
    scale_fill_manual(values = cols_cell) +
    
    #facet_wrap(~ Sample, nrow = 3)  #+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) + 
    theme(legend.position='bottom') + 
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) + 
    theme(panel.background = element_blank()) + 
    facet_wrap(~ data, nrow = 2)
  
  
  
  
  p2 <- 
    all_ready %>% 
    
    ggplot(aes(x = "", y = Freq, fill = FC)) + 
    geom_bar(stat = "identity", width = 1) +
    
    coord_polar(theta = "y") +
    
    #scale_fill_viridis_b() +
    #scale_fill_brewer(palette = "Spectral") +
    scale_fill_gradient2(low = "#163771", mid = "#ffffe0", high = "#751232") +
    
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) + 
    theme(legend.position='bottom') + 
    scale_y_continuous(trans = "reverse") + 
    theme(panel.background = element_blank()) + 
    facet_wrap(~ data, nrow = 2)
  
  
  p1 + p2
}


p <- plot_props()

p %>% save_x(data = ., name = "CellType_Prop_Piechart_ALL", 1, 10, 10, svg = T)









atac_epi <- df_atac %>% filter(predictedCellTypeGroup %in% c("Basal, LUM_HR-pos", "LUM_HR-neg"))


















celltype = "LUM_HR-pos"

p1 <- 
Treatment_Response %>% 
  filter(Cluster == celltype) %>% 
  top_n(25, avg_logFC) %>% 
  
  ggplot(aes(y = reorder(Gene, avg_logFC), x = avg_logFC)) + 
  geom_col(fill = "#eb933e") + 
  geom_text(aes(label = Gene), hjust = 1.5, size = 8) +
  
  theme(panel.background = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(), 
        panel.grid.major.x = element_line(color = "grey60")
        
        
  )

p2 <- 
Treatment_Response %>% 
  filter(Cluster == celltype) %>% 
  top_n(-25, avg_logFC) %>% 
  
  ggplot(aes(y = reorder(Gene, avg_logFC), x = avg_logFC)) + 
  geom_col(fill = "#9d3396") + 
  geom_text(aes(label = Gene), hjust = -0.5, size = 8) +
  
  theme(panel.background = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey60")
        
        )


p <- p2 + p1 

p %>% save_x(data = ., name = paste0("DE_Genes_Top40_", celltype), 1, 10, 10, svg = T)


CUX2 %>% save_x(data = ., name = paste0("CUX2_Umap"), 1, 10, 8, svg = F)













celltype = "LUM_HR-pos"

genes_u <- Treatment_Response %>% filter(Cluster == celltype) %>% top_n(25, avg_logFC) %>% pull(Gene)
genes_d <- Treatment_Response %>% filter(Cluster == celltype) %>% top_n(-25, avg_logFC) %>% pull(Gene)
genes <- c(genes_u, genes_d)


files <- 
list.files("ATAC_Integration/geneBrowserPlots 2/PDF/") %>% as.data.frame() %>% 
  separate(1, into = c("A", "B"), sep = "Plot-Tracks-Marker-Genes-byCellType_", remove = F) %>% 
  separate(B, into = c("C", "D"), sep = ".pdf") %>% #head() %>% 
  filter(C %in% genes) %>% pull(1) %>% as.character()


dir.create(paste0("ATAC_Integration/geneBrowserPlots 2/PDF/", celltype))

file.copy(from = paste0("ATAC_Integration/geneBrowserPlots 2/PDF/", files), 
          to   = paste0("ATAC_Integration/geneBrowserPlots 2/PDF/", celltype), 
          #overwrite = recursive, 
          recursive = FALSE,
          copy.mode = TRUE, 
          copy.date = FALSE)






Prop <- 
as.data.frame(prop.table(table(Sobj_integrated$CellType, 
                               Sobj_integrated$Type), 
                         margin = 2)) %>% 
  dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq") %>% 
  mutate(Type = substring(.$Sample, 1, 2))# %>% arrange(factor(Cluster, levels = levels), desc(Freq))

p1 <- 
cluster.averages %>% 
  filter(Gene %in% c(ligan, "EGF"), !is.na(avg_logFC)) %>% 
  filter(Gene %notin% c("HSP90AA1")) %>% 
  left_join(Freq) %>% 
  ggplot(aes(y = Cluster, x = Gene, color = avg_logFC)) + 
  
  geom_point(aes(size = Freq), stroke = 0) + 
  ggtitle(label = "RTK ligand output by celltype proportion") +
  
  scale_color_gradient2(low = "#00494b", mid = "#ffffe0", high = "#79260b") +
  
  coord_flip() +
  
  theme(text = element_text(family = "Lato"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12), 
        title = element_text(size = 15, face = "bold"),legend.position = "top",
        panel.grid.major = element_line(color = "grey95"), axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white", color = "grey40", size = 2))
        
  #pivot_longer(cols = c(CF, TM)) %>% unite(Cluster, Cluster, name) %>% ggplot(aes(x = Cluster, y = Gene, color = avg_logFC)) + geom_point(aes(size = value)) + theme(axis.text.x = element_text(angle = -25, hjust = 0))





p2 <- 
cluster.averages %>% 
  filter(Gene %in% c(recpt, "EGFR", "ERBB2")) %>% 
  
  left_join(Freq) %>% pivot_longer(c(CF, TM)) %>% 
  unite(Cluster, Cluster, name) %>% 
  mutate(avg_logFC = replace_na(avg_logFC, 0), expression = Freq * value) %>% 
  
  ggplot(aes(y = Cluster, x = Gene, color = avg_logFC)) + 
  
  geom_point(aes(size = expression), stroke = 2) + 
  ggtitle(label = "RTK receptor expression by celltype") +
  
  scale_color_gradient2(low = "#00494b", mid = "grey90", high = "#79260b") +
  
  coord_flip() +
  
  theme(text = element_text(family = "Lato"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12), legend.position = "top",
        title = element_text(size = 15, face = "bold"), axis.title.y = element_blank(),
        panel.grid.major = element_line(color = "grey95"),
        panel.background = element_rect(fill = "white", color = "grey40", size = 2))



p1 + p2






obs <- read_csv("~/Downloads/obsRaw.csv")

obs$CellType <- factor(obs$CellType, levels = levels)

df <- 
obs %>% 
  select(CellType, Type, ratio = `spliced%`) %>% 
  group_by(CellType, Type) %>% 
  mutate(avg_ratio = mean(ratio)) %>% 
  select(-ratio) %>% distinct()

df2 <- 
obs %>% 
  select(CellType, Type, ratio = `spliced%`) %>% 
  group_by(Type) %>% 
  mutate(avg_ratio = mean(ratio)) %>% 
  select(-ratio, -CellType) %>% distinct()



p2 <- 
df %>% 
  ggplot(aes(y = factor(CellType, levels = rev(levels)), x = avg_ratio, fill = Type)) + 
  
  geom_col(position = position_dodge()) + 
  
  geom_text(aes(label = round(avg_ratio, digits = 1)), 
            position = position_dodge(width = 1)) +
  
  scale_fill_manual(values = cols_type) + 
  
  theme(
    panel.background = element_rect(fill  = "white",
                                    color = "grey60",
                                    size  = 1),
    panel.grid.major.x = element_line(color = "grey60"),
    axis.title = element_blank(),
    legend.position = "none")




p1 <- 
df2 %>% 
  ggplot(aes(y = Type, x = avg_ratio, fill = Type)) + 
  
  geom_col(position = position_dodge()) + 
  
  geom_text(aes(label = round(avg_ratio, digits = 1)), 
            position = position_dodge(width = 1)) +
  
  scale_fill_manual(values = cols_type) + 
  
  theme(
    panel.background = element_rect(fill  = "white",
                                    color = "grey60",
                                    size  = 1),
    panel.grid.major.x = element_line(color = "grey60"),
    axis.title = element_blank(),
    legend.position = "none")


p <- 
p1 + p2 + patchwork::plot_layout(ncol = 1, heights = c(1, 5))


p %>% save_x(data = ., name = "Splicing_Ratio", 1, 5, 10, svg = T)









p <- 
Sobj_integrated@meta.data %>% 
  select(Sample, GPAM_1, Type) %>% 
  group_by(Type, Sample) %>% 
  summarise(avg = mean(GPAM_1)) %>% 
  #filter(Sample %in% keep) %>% 
  
  ggplot(aes(x = Type, y = avg, fill = Type)) +   
  
  geom_boxplot(color = "grey30", 
               width = 0.5, 
               alpha = 0.8, 
               position=position_dodge(0.9), 
               lwd = 1,
               outlier.alpha = 0) +
  
  geom_point(aes(color = Sample),
             position = position_jitterdodge(jitter.width = 0.05, 
                                             jitter.height = 0,
                                             dodge.width = 0.9), 
             alpha = 1,
             size  = 3) +
  
  stat_compare_means(paired = F, method = "wilcox", 
                     label.x.npc = 0.25, vjust = -1) +
  
  scale_fill_viridis_d(option = "C", begin = 0.3, end = 0.8) +
  #scale_color_viridis_d(begin = 0.3, end = 0.6) +
  scale_color_manual(values = cols) +
  
  ggtitle(label = "GPAM Module Score") +
  
  theme(text = element_text(family = "Lato"),
        axis.text = element_text(size=10), 
        #axis.title=element_text(size=14,face="bold"),
        axis.title = element_blank(), 
        #panel.background = element_rect(fill = "grey95"),
        panel.background = element_rect(fill = "white", color = "grey60", size = 2),
        panel.grid = element_line(color = "grey90"),
        strip.text = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(size=14, angle = -25, hjust = 0, vjust = 1, face = "bold")
  ) +
  scale_y_sqrt()


p %>% save_x(data = ., name = "GPMA_Module_Boxplot", 1, 5, 6, svg = T)










p <- 
Prop %>% 
  filter(Cluster == "Lymphatic_EC") %>% 
  ggplot(aes(x = Cluster, y = Freq, fill = Type)) + 
  
  geom_boxplot(color = "grey30", 
               width = 0.5, 
               alpha = 0.8, 
               position=position_dodge(0.9), 
               lwd = 1,
               outlier.alpha = 0) +
  
  geom_point(aes(color = Sample),
             position = position_jitterdodge(jitter.width = 0.05, 
                                             jitter.height = 0,
                                             dodge.width = 0.9), 
             alpha = 1,
             size  = 3) +
  
  stat_compare_means(paired = F, method = "wilcox", 
                     label.x.npc = 0.25, vjust = -1) +
  
  scale_fill_viridis_d(option = "C", begin = 0.3, end = 0.8) +
  #scale_color_viridis_d(begin = 0.3, end = 0.6) +
  scale_color_manual(values = cols) +
  
  ggtitle(label = "Proportion of Vasculature") +
  
  theme(text = element_text(family = "Lato"),
        axis.text = element_text(size=10), 
        #axis.title=element_text(size=14,face="bold"),
        axis.title = element_blank(), 
        #panel.background = element_rect(fill = "grey95"),
        panel.background = element_rect(fill = "white", color = "grey60", size = 2),
        panel.grid = element_line(color = "grey90"),
        strip.text = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(size=14, angle = -25, hjust = 0, vjust = 1, face = "bold")
  )



p %>% save_x(data = ., name = "Lymphatic_EC_Proportion", 1, 5, 6, svg = T)









#plot_prop <- function(query) { 
  
  
  data.subset <- Sobj_integrated@meta.data #%>% filter(CellType %in% query)
  data.subset$Subcluster <- as.character(data.subset$Subcluster)

Prop <- 
  as.data.frame(prop.table(table(data.subset$Subcluster, 
                                 data.subset$Sample), 
                           margin = 2)) %>% 
  dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq") %>% 
  mutate(Type = substring(.$Sample, 1, 2))# %>% arrange(factor(Cluster, levels = levels), desc(Freq))


p3 <- 
Prop %>% filter(Cluster %in% as.character(unique(filter(Sobj_integrated@meta.data, CellType == "Vasc_Accessory")$Subcluster))) %>% 
  ggplot(aes(x = reorder(Cluster, - Freq), y = Freq, fill = Type)) + 
  
  geom_boxplot(color = "grey30", 
               width = 0.5, 
               alpha = 0.8, 
               position=position_dodge(0.9), 
               lwd = 2,
               outlier.alpha = 0) +
  
  geom_point(aes(color = Sample),
             position = position_jitterdodge(jitter.width = 0.05, 
                                             jitter.height = 0,
                                             dodge.width = 0.9), 
             alpha = 1,
             size  = 4) +
  
  stat_compare_means(paired = F, method = "wilcox", 
                     label.x.npc = 0.25, vjust = -1) +
  
  scale_fill_viridis_d(option = "C", begin = 0.3, end = 0.8) +
  #scale_color_viridis_d(begin = 0.3, end = 0.6) +
  scale_color_manual(values = cols) +
  
  ggtitle(label = "Proportion of Vasculature") +
  
  theme(text = element_text(family = "Lato"), 
        title = element_text(size = 30), legend.position = "none",
        axis.text = element_text(size=30), 
        #axis.title=element_text(size=14,face="bold"),
        axis.title = element_blank(), 
        #panel.background = element_rect(fill = "grey95"),
        
        panel.grid.major = element_line(color = "grey90"),
        panel.background = element_rect(fill = "white", color = "grey60", size = 4),
        strip.text = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(size=30, angle = -25, hjust = 0, vjust = 1, face = "bold") 
  )
#}


p <- p1 + p2 + p3



p %>% save_x(data = ., name = "Vasculature_Proportion", 1, 20, 10, svg = T)










gProfiler_KEGG_PI3K_AKT <- read_csv("utilities/gProfiler_KEGG_PI3K-AKT.csv")

PI3_TFs <- gProfiler_KEGG_PI3K_AKT %>% filter(name %in% TFs$Gene) %>% pull(name)


mat <- 
Treatment_Response %>% 
  filter(Gene %in% c(PI3_TFs, "KLF6", "ZBTB7A", "INSR", "MBNL2", "EGR1", "TRPS1", "IRF3")) %>% 
  group_by(Gene) %>% 
  mutate(n = n()) %>% 
  filter(n >= 5) %>% 
  select(Gene, avg_logFC, Cluster) %>% 
  pivot_wider(names_from = Cluster, values_from = avg_logFC)# %>% 
  #column_to_rownames("Gene")


GTEX <- dflong %>% filter(str_detect(Tissue, "Breast|Prostate"), hgnc_symbol %in% mat$Gene) %>% select(hgnc_symbol, Tissue, value) %>% pivot_wider(names_from = Tissue, values_from = value) %>% mutate(FC = log(Prostate/Breast...Mammary.Tissue)) %>% select(Gene = "hgnc_symbol", GTEX = "FC")

mat <- mat %>% left_join(GTEX) %>% column_to_rownames("Gene")

mat[is.na(mat)] <- 0

my.breaks <- c(seq((min(mat)), -0.00001, length.out = 50), seq(0.00001, (max(mat)), length.out = 50))
my.colors <- c(colorRampPalette(colors = c("#002c2d", "#00494b", "#076769", "#208288", "#3ea8a6", "#76c7be", "#bbe4d1", "#ffffe0"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("#ffffe0", "#fed693", "#f5ad52", "#dd8629", "#bd651a", "#9c4511", "#79260b", "#580000"))(length(my.breaks)/2))


p <- pheatmap(mat, breaks = my.breaks, color = my.colors)


p %>% save_x(data = ., name = "TFs_PI3K_Heatmap", 1, 10, 10, svg = T)







p1 <- 
  Sobj_integrated@meta.data %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(color = "grey90", size = 4) + 
  geom_point(data = filter(Sobj_integrated@meta.data, Type == "TM"), 
             aes(color = Sample), size = 2) + 
  scale_color_manual(values = cols[9:18]) + 
  theme(panel.background = element_blank(), axis.title = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank())



p2 + p1 + patchwork::plot_layout(ncol = 1)







df %>% ggplot(aes(x = CellType, y = Insulin, fill = Type)) + 

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
  
  #facet_wrap(vars(Module), scales = "free_y") +
  
  #scale_y_sqrt() +
  
  theme(text = element_text(family = "Lato", size = 20),
        axis.text.x = element_text(angle = -25, hjust = 0, vjust = 1),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray95"),
        panel.border = element_rect(fill = NA, color = "gray30", size = 2),
        strip.text = element_text(size = 13, face = "bold", hjust = 0), 
        legend.key.size = unit(1,"cm"), legend.justification = c(1,1))









p <- 
  Sobj_integrated@meta.data %>% 
  select(Insulin_1, CellType, Type) %>% 
  ggplot(aes(x = factor(CellType, levels = levels), y = Insulin_1, fill = Type)) + 
  
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
  
  theme(text = element_text(family = "Lato"),
        axis.text.x = element_text(face = "bold", size = 15, angle = -25, hjust = 0, vjust = 1),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "gray96"),
        strip.text = element_text(size = 13, face = "bold", hjust = 0), 
        legend.key.size = unit(1,"cm"), legend.justification = c(1,1))


p %>% save_x(data = ., name = "Insulin_Sigaling_Adipocytes", 1, 10, 7, svg = T)


genes %>% length()



