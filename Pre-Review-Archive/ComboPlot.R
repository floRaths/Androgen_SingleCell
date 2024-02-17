library(scales)

levels <- c("CF-3920", "CF-19301", "CF-7780", "CF-2797", "CF-318-813", "CF-0404", "CF-428-112", "CF-249-347", "CF-4014", "TM-9469", "TM-1956", "TM-6544", "TM-6477", "TM-8249", "TM-7567", "TM-9817", "TM-2768", "TM-3937")
cols   <- c("#F26C66", "#FFBB78", "#F57F46", "#7FC210", "#2CA02C", "#BCBD22", "#EDBD1F", "#A8486A", "#C49C94", "#9467BD", "#C5B0D5", "#17BECF", "#C5D4D2", "#438ABB", "#FA7AA3", "#F7B6D2", "#FF9896", "#3FB8AF")


Idents(Sobj_integrated) <- "Sample"
levels(Sobj_integrated) <- levels



combo_plot <- function(Sobj, orient, pt1, pt2) {
  col <- pal_d3("category20")(20)
  

  if(orient == "horz") { 
    
    plot1 <- 
      DimPlot(Sobj, 
              reduction = "umap", 
              group.by  = "Subcluster",
              pt.size = pt1, 
              label = T,
              repel = T, 
              label.size = 6, cols = col
              
      ) + theme(axis.title = element_blank(), 
                axis.ticks = element_blank(), 
                axis.text = element_blank(),
                legend.position = "bottom",
                text =  element_text(family = "Lato")) 
    
    
    plot2 <- 
      DimPlot(Sobj, 
              reduction = "umap", 
              group.by  = "Sample",
              split.by = "Type",
              pt.size = pt2, 
              label = F,
              repel = T, 
              ncol = 1, 
              cols = cols,
              
      ) + theme(axis.title = element_blank(),
                axis.ticks = element_blank(), 
                axis.text = element_blank(),
                legend.position = "bottom",
                text =  element_text(family = "Lato")) 
    
    
    
    plot1 + plot2 + patchwork::plot_layout(widths = c(2,1), ncol = 2)
  } else {

  if(orient == "vert") { 
    plot1 <- 
      DimPlot(Sobj, 
              reduction = "umap", 
              group.by  = "Subcluster",
              pt.size = pt1, 
              label = T,
              repel = T, cols = col,
              label.size = 6,
              
      ) + theme(axis.title = element_blank(),
                axis.ticks = element_blank(), 
                axis.text = element_blank(),
                legend.position = "bottom",
                text =  element_text(family = "Lato")) 
    
    plot2 <- 
      DimPlot(Sobj, 
              reduction = "umap", 
              group.by  = "Sample",
              split.by = "Type",
              pt.size = pt2, 
              label = F,
              repel = T, 
              ncol = 2,
              cols = cols,
              
      ) + theme(axis.title = element_blank(),
                axis.ticks = element_blank(), 
                axis.text = element_blank(),
                legend.position = "bottom",
                text =  element_text(family = "Lato")) 
    
    plot1 + plot2 + patchwork::plot_layout(heights = c(2,1), ncol = 1)
  }}
}

celltype = "Adipocyte"


#image <- 
combo_plot(Sobj = Sobj_integrated, 
           orient = "horz", 
           pt1 = 1, pt2 = 0.5)



image %>% save_x(data = ., name = paste0(celltype, "_Combo_Plot"), 1, 16, 10, svg = F)



features <- c("AHR_Module_1", "ANXA1_Module_1", "FOXO3_Module_1", "JAZF1_Module_1")
plot_list    <- vector("list", length = length(features))

for (i in 1:length(features)) {
  
plot_list[[i]] <- FeaturePlot(Sobj_integrated, reduction = "umap",
            #cells = WhichCells(Sobj_integrated, expression = SCORE_1 >= 0.25),
            #features = c("PLXNA4", "HPSE2", "SEMA3C", "EPHA3", "SLC22A3", "ZBTB7C"), 
            features = features[i],
            #split.by = "Type",
            pt.size = 1, 
            order = T, 
            #max.cutoff = 0.45,
            cols = viridis::magma(n = 100),
            #cols = c("grey", "Green", "Red")
            ) + theme(axis.title = element_blank(),
                      axis.text = element_blank(),
                      axis.ticks = element_blank(), 
                      legend.key.size = unit(0.3, "cm"),
                      legend.position = c(0.85, 1)
                      )
  
}

cowplot::plot_grid(plotlist = plot_list, ncol = 1)







margin = 0.5

look = theme(panel.background = element_rect(fill  = "white", 
                                             color = "grey60", 
                                             size  = 1), 
             panel.grid = element_blank(), plot.margin = unit(c(margin,margin,margin,margin), "cm"),
             axis.title = element_blank(), 
             axis.text  = element_blank(),
             axis.ticks = element_blank(),
             legend.position = "none")





plot_global_umaps(0.5, 2, 0) %>% save_x(data = ., name = "RNA_ATAC_Global_UMAPs", 1.5, 10, 8, svg = T)





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
  
  facet_wrap(~Type, ncol = 1) +
  
  scale_color_manual(values = cols) + 
  
  look + 
  theme(strip.background = element_blank(), panel.background = element_blank(),
               strip.text = element_blank())


P1 + as.ggplot(P2) + patchwork::plot_layout(ncol = 2, widths = c(2,1)) 


P1 + P2 + patchwork::plot_layout(ncol = 2, widths = c(2,1)) 
P1 + P2 + patchwork::plot_layout(ncol = 1, heights = c(2,1)) 




P2 %>% save_x(data = ., name = paste0(celltype, "Subcluster_Umap"), 0.8, 6, 10, svg = F)



