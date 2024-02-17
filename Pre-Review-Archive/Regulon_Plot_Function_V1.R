# Data Loading ------------------------------------------------------------

category <- "CellType_All"
Idents(Sobj_integrated) <- category

input_vargenes <- "4k"

adjacencies <- read_tsv(paste0("~/Documents/SCENIC/pySCENIC_Full/scenicdata/TFs_VarFeatures/Results_", input_vargenes, "/expr_mat.adjacencies.tsv"))

percentile <- 
  adjacencies %>% 
  group_by(TF) %>% 
  mutate(percentile_rank = ntile(importance, 100))


GRN_Scores <- 
  readRDS(paste0("~/Documents/SCENIC/pySCENIC_Full/scenicdata/TFs_VarFeatures/Results_", input_vargenes, "/GRN_Scores_99_percentile.rds"))  %>% 
  dplyr::select(-c(1:22))#%>% column_to_rownames("Cell")

### remove (+) tail from colnames
names(GRN_Scores)[names(GRN_Scores) == colnames(GRN_Scores)] <- 
  paste(substr(colnames(GRN_Scores), 1, nchar(colnames(GRN_Scores))-2))



Sobj_plot   <- AddMetaData(Sobj_integrated, metadata = GRN_Scores)

MetaData    <- readRDS("Output/MetaData_May05.rds")

Wilcox_bind <- readRDS(paste0("Output/Wilcox_bind_", category, "_", input_vargenes, "_May05.rds"))

Marker_Bind <- do.call(rbind.data.frame, readRDS(paste0("Output/TM-CF_Response_", category, "_MAST_list.rds"))) %>% 
  filter(p_val_adj <= 0.05) %>% group_by(Gene) %>% 
  mutate(events = n()) %>% arrange(-events)




#Idents(Sobj_integrated) <- "CellType"
#celltype <- c("LUM_HR-pos")



##################################################################################################################

regulon_plot <- function(celltype, fc.cutoff, n.tfs, coreg.TFs, rank, gap, save) {

  
dir.create(paste0("Figures/5.5_Regulon_Plot_Pipeline/", celltype), showWarnings = F, recursive = T)  
figure.path <- paste0("Figures/5.5_Regulon_Plot_Pipeline/", celltype)
  
  save.curve <- function(name, scale) {
    
    ggsave(plot     = curve,
           path     = figure.path,
           filename = paste0(name,"_Regulon_Curve.png"),
           width    = 16, 
           height   = 9,
           scale    = scale,
           dpi      = 300)
  }
  save.regulon_plot <- function(name, scale) {
    
    ggsave(plot     = out,
           path     = figure.path,
           filename = paste0(name,"_Regulon_Plot.png"),
           width    = 16, 
           height   = 9,
           scale    = scale,
           dpi      = 300)
  }  
  
  
# outlier detection function. adjust k for more stringent outlier detection  
isnt_out_tukey <- function(x, k = 1.5, na.rm = TRUE) {
  quar <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  iqr <- diff(quar)
  
  (quar[1] - k * iqr <= x) & (x <= quar[2] + k * iqr)
  }  
  

celltype <<- celltype

# pull out cell IDs of target population      
cells <<- WhichCells(Sobj_integrated, idents = celltype)

# collect DE-genes and genes from Correlation analyis
targets <<- c(filter(Marker_Bind, CellType %in% celltype, 
                     avg_logFC >= fc.cutoff | avg_logFC <= -fc.cutoff)$Gene, 
              filter(Corr_Bind, CellType %in% celltype)$Gene) %>% unique()

# find the TFs that regulate the DE targets. Targets are in 99th percentile of TF-target pool
# and from that pool we take the top 5 most important TFs

#percentile %>% 
#  filter(target %in% targets, percentile_rank > rank) %>%
#  group_by(target) %>% 
#  top_n(coreg.TFs, importance) %>% 
#  pull(TF) %>% unique() %>% 
#  intersect(pull(filter(Wilcox_bind, CellType %in% celltype, Wilcox <= 0.05), TF)) %>% 
#  unique() %>% 
#  intersect(colnames(GRN_Scores)) %>% 
#  unique() %>% sort()

celltype.tfs <<- 
Average_TF_Expression %>% ungroup() %>% 
  filter(expression == "expressed", CellType_All == celltype, perc >= 30) %>% pull(TF) %>% unique() %>% 
  intersect(colnames(GRN_Scores)) %>% sort()

# filtering out TF regulons where the TM-CF difference is not significant
# now we only keep the significant ones
#celltype.tfs <<- Wilcox_bind %>% pull(filter(Wilcox_bind, CellType %in% celltype, Wilcox <= 0.05), TF) %>% unique() 




#GRN_Summary <- 
#GRN_Scores[cells, ] %>% dplyr::select(celltype.tfs) %>%
#  rownames_to_column("Cell") %>% 
#  left_join(dplyr::select(rownames_to_column(MetaData, "Cell"), Cell, Type), by = "Cell") %>%
#  group_by(Type) %>% summarise_at(vars(celltype.tfs), mean) %>%
#  pivot_longer(cols = celltype.tfs, names_to = "TF") %>% 
#  pivot_wider(names_from = Type, values_from = value) %>%
#  mutate(change = TM - CF) %>% 
#  arrange(-change) 




#############################################################################################################
# building the summary of how strong TF regulon scores are changing between TM CF

GRN_Summary <<- 
  GRN_Scores[cells, ] %>% dplyr::select(celltype.tfs) %>%                        # take the significant tfs and the cells
  rownames_to_column("Cell") %>%                                                      
  left_join(dplyr::select(rownames_to_column(MetaData, "Cell"), # join with metadat (this can be clenaed up)
                          Cell, Type, Sample), by = "Cell") %>%
  group_by(Sample, Type) %>% summarise_at(vars(celltype.tfs), mean) %>%          # calculate mean regulon scores per Sample
  pivot_longer(cols = celltype.tfs, names_to = "TF") %>% 
  group_by(Type, TF) %>%
  mutate(mean.out = mean(value), is_outlier = !isnt_out_tukey(value)) %>%        # detect outliers within their groups and calculate mean of whole Type
  group_by(Type, TF, is_outlier) %>%                                             # separate outliers and non-outliers
  mutate(mean.in = mean(value)) %>%                                              # calculate mean without outliers
  filter(is_outlier == "FALSE") %>%                                              # filter outliers out
  ungroup() %>% 
  dplyr::select(Type, TF, mean.in) %>% distinct() %>%                            # now we have the "uncontaminated" mean values we were looking for
  pivot_wider(names_from = Type, values_from = mean.in) %>%
  mutate(change = TM - CF) %>%                                                   # calculate difference between TM and CF, 
  arrange(-change)                                                               # which allows us to select strong candicted




#############################################################################################################
# this summarizes the mean per sample which will be plotted next to the main plot to see sample variability

GRN_Samp_Summ <<- 
  GRN_Scores[cells, ] %>% dplyr::select(celltype.tfs) %>%
  rownames_to_column("Cell") %>% 
  left_join(dplyr::select(rownames_to_column(MetaData, "Cell"), Cell, Sample), by = "Cell") %>%
  group_by(Sample) %>% summarise_at(vars(celltype.tfs), mean) %>%
  pivot_longer(cols = celltype.tfs, names_to = "TF") %>% 
  pivot_wider(names_from = TF, values_from = value) #%>% column_to_rownames("TF")





# pick out top and bottom most changing TFs .. this needs to be improved
changing.tfs <<- 
rbind(top_n(GRN_Summary, (n.tfs/2), change), top_n(GRN_Summary,-(n.tfs/2), change))

changing.tfs



# plot the change... helps to see if some cutoff need to be adjusted
curve <- 
GRN_Summary %>%
  ggplot(aes(x = reorder(TF, change), y = change)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8), 
                       plot.title = element_text(size = 17, face = "bold")) + 
  ggtitle(paste("GRN Score_changes...", celltype)) + 
  geom_hline(yintercept = 0, color = "grey")

plot(curve)


if(save == T){ save.curve(celltype, 1.5) }

# pick the target genes that belong to the changers and that are DE in the data
variable.targets <- 
percentile %>% 
  filter(TF %in% changing.tfs$TF, percentile_rank > 99) %>% 
  filter(target %in% targets) %>% 
  pull(target) %>% unique() %>% sort()




#############################################################################################################
# Build the Heatmaps

# first heatmap shows the GRNboost TF importance values for the selected pairs
mat <- 
adjacencies %>% 
  filter(TF %in% changing.tfs$TF, target %in% variable.targets) %>% 
  pivot_wider(names_from = TF, values_from = importance) %>%
  column_to_rownames("target")

# second heatmap show the overall score of that regulon in each Sample
mot <- 
GRN_Samp_Summ %>% dplyr::select(Sample, colnames(mat)) %>% column_to_rownames("Sample")


mat[is.na(mat)] <- 0
mot[is.na(mot)] <- 0

mat <<- mat


# annotations ...
#M <- filter(dplyr::select(Marker_Bind, Gene, CellType, avg_logFC), CellType %in% celltype)
#C <- filter(dplyr::select(Corr_Bind,   Gene, CellType, avg_logFC = "Log_FC"), CellType %in% celltype, Gene %notin% M$Gene)

#M <- filter(dplyr::select(Marker_Bind, Gene, CellType, DE.log_FC = "avg_logFC"), CellType %in% celltype)
#C <- filter(dplyr::select(Corr_Bind,   Gene, CellType, Corr.FC = "Log_FC"), CellType %in% celltype)#, Gene %notin% M$Gene)

Ann <- 
filter(dplyr::select(Marker_Bind, Gene, CellType, DE.log_FC = "avg_logFC"), CellType %in% celltype) %>% 
  full_join(filter(dplyr::select(Corr_Bind,   Gene, CellType, Corr.FC = "Log_FC"), CellType %in% celltype), by = c("Gene", "CellType")) %>% 
  dplyr::select(-CellType) %>% column_to_rownames("Gene")

#Ann <- 
#  bind_rows(M, C) %>%
#  dplyr::select(-CellType) %>% column_to_rownames("Gene")

Ann.C <- 
GRN_Summary %>% dplyr::select(-CF, - TM) %>% column_to_rownames("TF")


### main plot
image <- 
  pheatmap(t(mat), scale = "column", fontsize = 10, 
           border_color = NA, 
           cluster_cols = T,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)),
           annotation_col = Ann,
           #annotation_row = Ann.C, 
           silent = T,
           #annotation_colors = ann_colors, 
           treeheight_row = 0, treeheight_col = 0, 
           cutree_rows  = 1, 
           cutree_cols  = 1,
           fontsize_row = 12,
           fontsize_col = 10, main = paste("TF-Target Importance - ", celltype, "             (n = ", length(rownames(mat)),")")
           
  ) 

# now we take the row order from the main plot and reorder the matrix for the second plot (which will not apply clustering)
met <<- 
mot %>% dplyr::select(image$tree_row$order)


# annotation ...
Enn <- 
MetaData %>% dplyr::select(Sample, Type) %>% distinct() %>% column_to_rownames("Sample")

ann_colors = list(
  Type = c(TM = "lightseagreen", CF = "maroon3")
  )


### secondary plot
imoge <- 
  pheatmap(t(met), scale = "row", fontsize = 10, 
           border_color = NA, 
           cluster_rows = F, 
           cluster_cols = F, 
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           color = viridis::magma(n = 100),
           #annotation_row = Ann,
           annotation_col = Enn, silent = T,
           annotation_colors = ann_colors, 
           treeheight_row = 0, treeheight_col = 0, 
           cutree_rows = 1, 
           cutree_cols = 1,
           fontsize_col = 10, 
           legend = F, gaps_col = gap,
           show_rownames = F, main = "TF Regulon Score"
           
  )


# combine the plots in a grid
out <<- plot_grid(NULL, NULL, as.grob(imoge), as.grob(image), 
          ncol = 2, nrow = 2, 
          rel_widths  = c(1,5), 
          rel_heights = c(1,40),
          align = "h"#, #label_x = 0.1,
          #labels = c("TF Module Scores", paste("TF-Target Importance")), hjust = -0.05
          )

# return plot

if(save == T){ save.regulon_plot(celltype, 1.5) }

out

}


regulon_plot("Adipocyte", 
             n.tfs = 30, 
             fc.cutoff = 0.35,
             coreg.TFs = 5, 
             rank      = 99, 
             gap = 9, 
             save = F)






#------------------------------------------

VlnPlot(Sobj_plot, idents = celltype, features = c("CUX2", "PGR", "BATF", "ZNF689"), ncol = 2, pt.size = 0, split.by = "Type", group.by = "Sample")

VlnPlot(Sobj_plot, idents = "LUM_HR-pos", features = c("SF1"), ncol = 1, pt.size = 0, split.by = "Type", group.by = "Sample")


GRN_Summary %>% 
  pivot_longer(cols = c(2,3), names_to = "Type") %>%
  group_by(TF) %>% do(w = wilcox.test(value~Type, data = ., paired=FALSE)) %>% 
  summarise(TF, Wilcox = w$p.value)




Wilcox <<- 
  GRN_Scores[cells, ] %>% dplyr::select(celltype.tfs) %>%
  rownames_to_column("Cell") %>% 
  left_join(dplyr::select(rownames_to_column(MetaData, "Cell"), Cell, Sample, Type), by = "Cell") %>% 
  pivot_longer(cols = celltype.tfs, names_to = "TF") %>% dplyr::select(-1) %>% #arrange(TF)
  group_by(TF) %>% do(w = wilcox.test(value~Type, data = ., paired=FALSE)) %>% 
  summarise(TF, Wilcox = w$p.value)


Wilcox %>% filter(Wilcox <= 0.05) %>% arrange(Wilcox)









Idents(Sobj_integrated) <- "CellType_All"
#celltype <<- "CD8_T"

celltype  <- unique(levels(Sobj_integrated))
Wilcox_list <- vector("list", length = length(celltype))  


for (i in 1:length(celltype)) {
  
cells  <- WhichCells(Sobj_integrated, idents = celltype[i])

Wilcox_list[[i]] <- 
  GRN_Scores[cells, ] %>%
  rownames_to_column("Cell") %>% 
  left_join(dplyr::select(rownames_to_column(MetaData, "Cell"), Cell, Sample, Type), by = "Cell") %>% 
  pivot_longer(cols = colnames(GRN_Scores), names_to = "TF") %>% dplyr::select(-1) %>%
  group_by(TF) %>% do(w = wilcox.test(value~Type, data = ., paired=FALSE)) %>% 
  summarise(TF, Wilcox = w$p.value) %>% 
  add_column(CellType = celltype[i])

}

names(Wilcox_list) <- celltype
Wilcox_bind <- do.call(rbind.data.frame, Wilcox_list)

saveRDS(Wilcox_bind, "Output/Wilcox_bind_CellType_ALL_4k_May05.rds")






