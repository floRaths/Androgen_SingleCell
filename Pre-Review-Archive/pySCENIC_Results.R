library(readr)
library(pheatmap)
library(RColorBrewer)

setwd("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/")

adjacencies <- read_tsv("~/Documents/SCENIC/pySCENIC_Full/scenicdata/CT_All_4kvar/Results/Basal/expr_mat.adjacencies.tsv")
adjacencies <- read_tsv("~/Documents/SCENIC/pySCENIC_Full/scenicdata/TFs_VarFeatures/Results_4k/expr_mat.adjacencies.tsv")
adjacencies <- read_tsv("~/Documents/SCENIC/pySCENIC_Full/scenicdata/TM_CF_Separate/Results/CF/expr_mat.adjacencies.tsv")
auc_mtx     <- read_csv("~/Documents/SCENIC/pySCENIC_Full/scenicdata/4k_Var_Genes/Results/FullDataSet/auc_mtx.csv") %>% column_to_rownames("Cell")



# AUC_Processing ----------------------------------------------------------


### remove (+) tail from colnames
names(auc_mtx)[names(auc_mtx) == colnames(auc_mtx)] <- 
  paste(substr(colnames(auc_mtx),1,nchar(colnames(auc_mtx))-3), "4k", sep = "_")

AUC_4k <- auc_mtx

X2 <- 
rownames_to_column(AUC_05, "Cell") %>% 
  left_join(rownames_to_column(AUC_1k, "Cell"), by = "Cell") %>%
  left_join(rownames_to_column(AUC_2k, "Cell"), by = "Cell") %>%
  left_join(rownames_to_column(AUC_4k, "Cell"), by = "Cell") %>%
  left_join(rownames_to_column(AUC_8k, "Cell"), by = "Cell") %>% 
  left_join(rownames_to_column(AUC_8k_Ank, "Cell"), by = "Cell")
  



regulons    <- read_csv("~/Documents/SCENIC/pySCENIC_Full/scenicdata/8kVargenes_Results/regulons.csv", col_names = c("TF", "MotifID", "AUC", "Annotation",  "Context", "MotifSimilarity",   "NES_OrthologousIden", "RankAtMax", "TargetGenes"), skip = 3)

filter(as.data.frame(colnames(auc_mtx)), colnames(auc_mtx) %in% c("ESR1(+)", "NR4A1(+)", "AR(+)", "PGR(+)"))


cells <- WhichCells(Sobj, idents = "LUM_HR-pos")



mat <- column_to_rownames(filter(auc_mtx, Cell %in% cells), "Cell")
mat <- column_to_rownames(filter(X2, Cell %in% colnames(Sobj_integrated)), "Cell")

regs <- top_n(rownames_to_column(as.data.frame(colMeans(mat)), "Regulon"), 25, colMeans(mat)) 


#image <- 
pheatmap(mat, 
         treeheight_col = 0, 
         cluster_rows = T, 
         cluster_cols = T, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_colnames = T,
         show_rownames = F,
         main = "Immune",
         annotation_row =  dplyr::select(Sobj@meta.data, Type, Sample, CellType),
         annotation_colors = ann_colors
         )



save("Selected_Regulons_LUMpos", 1)


ann_colors = list(
  Type = c(TM = "forestgreen", CF = "goldenrod2")
  #CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
  #GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
)










# Adding AUC umap to Seurat Object ----------------------------------------
Idents(Sobj_integrated) <- "CellType"
celltype <- levels(Sobj_integrated)
celltype <- c("LUM_HR-pos")
cells <- WhichCells(Sobj_integrated, idents = celltype)

#Define cells
regulons <- MainTFs_Bind %>% unite(input, Gene, input) %>% filter(CellType %in% celltype) %>% pull(input) %>% unique()
regulons <- GRN_Scores_99percentile %>% filter(Cell %in% cells) %>% column_to_rownames("Cell") %>% colnames()



cells <- intersect(X2$Cell, WhichCells(Sobj_integrated, idents = celltype))
LUM <- subset(Sobj_integrated, cells = cells)

Sobj_integrated@reductions$SCENIC <- RunUMAP(assay = "RNA", column_to_rownames(filter(X2, Cell %in% cells), "Cell")[regulons], reduction.name = "SCENIC", reduction.key = "SCENIC_")
LUM@reductions$GRN <- RunUMAP(assay = "RNA", column_to_rownames(filter(GRN_Scores_99percentile, Cell %in% cells), "Cell"), reduction.name = "GRN", reduction.key = "GRN_")

#regs <- top_n(rownames_to_column(as.data.frame(colMeans(mat)), "Regulon"), 120, colMeans(mat)) 
#LUM <- AddMetaData(LUM, metadata = select(column_to_rownames(filter(auc_mtx, Cell %in% cells), "Cell"), regs$Regulon))
Sobj_integrated <- AddMetaData(Sobj_integrated, metadata = column_to_rownames(filter(X2, Cell %in% cells), "Cell")[regulons])


#p3 <- 
DimPlot(LUM, reduction = "GRN", pt.size = 1, group.by = "Sample", label = F, split.by = "Type", ncol = 5)
DimPlot(LUM, reduction = "umap", pt.size = 1, group.by = "Sample", split.by = "Type")

plot1 <- DimPlot(LUM, reduction = "SCENIC", pt.size = 1, group.by = "Type")
plot2 <- DimPlot(LUM, reduction = "umap", pt.size = 1, group.by = "Type")

image <- 
CombinePlots(list(plot2, plot1))


#regulons <- MainTFs_Bind %>% unite(input, Gene, input) %>% filter(CellType == "LUM_HR-pos") %>% pull(input) %>% unique()
#regulons <- colnames(LUM@meta.data[22:812])




#for (i in 1:length(regulons)){

quick <-  
  FeaturePlot(LUM, features = regulons[i], 
            reduction = "SCENIC",
            pt.size = 1, 
            order = T, split.by = "Type",
            cols = viridis::magma(n = 100))
save2(paste(regulons[i]), 1)



for (i in 1:length(regulons)){
p1 <- 
FeaturePlot(LUM, features = regulons[i], 
            reduction = "SCENIC",
            pt.size = 1, 
            order = T, split.by = "Type",
            cols = viridis::magma(n = 100), by.col = T, combine = T)
p2 <- 
VlnPlot(LUM, features = regulons[i], 
            #reduction = "SCENIC",
            pt.size = 0, 
            split.by = "Type")

quick <- p1 + p2 + p3 + plot_layout(nrow = 2)
save2(paste(regulons[i]), 1)
}

#image <- 
FeaturePlot(LUM, features = c("nCount_RNA"), 
            reduction = "SCENIC",
            pt.size = 1, 
            sort.cell = T, 
            cols = viridis::magma(n = 100))




save("AUC_umap_vs_RNA_Umap_LUMn", 0.75)







# Adjacencies Processing --------------------------------------------------

celltype <- c("LUM_HR-pos", "LUM_HR-neg", "Basal")

biglist <- c(filter(Corr_Bind, CellType %in% celltype)$Gene, filter(CTPlus_Bind, str_detect(Cluster, "^LUM") | str_detect(Cluster, "^Basa"))$Gene) %>% unique()
biglist <- c(filter(Corr_Bind, CellType %in% celltype)$Gene, filter(Marker_Bind, CellType %in% celltype)$Gene) %>% unique()
biglist <- c(Corr_Bind$Gene, TopTable$Gene, Marker_Bind$Gene, CTPlus_Bind$Gene) %>% unique()


percentile <- 
adjacencies %>% 
  group_by(TF) %>% 
  mutate(percentile_rank = ntile(importance, 100))


Corr_target <- 
  percentile %>% filter(percentile_rank > 90) %>% pull(target) %>% unique() %>% intersect(biglist)

Corr_TFs <- 
  percentile %>% filter(TF %notin% c("ANKRD30A", "ANKRD30B"), target %in% Corr_target) %>% group_by(target) %>% top_n(5, importance) %>% pull(TF) %>% unique()



pivot_wider(adjacencies, names_from = target, values_from = importance, names_repair = "unique")

filter(adjacencies, TF == "NR6A1")


adjacencies %>% group_by(TF) %>% top_n(50, importance) %>% 
  mutate(avg = mean(importance)) %>% arrange(-avg) %>% filter(TF == "ANKRD30A") %>% print(n = 35)


adjacencies %>% group_by(TF) %>% top_n(100, importance) %>% 
  mutate(avg = mean(importance)) %>% filter(target == "ANKRD30A")






topTFs2 <- adjacencies %>% group_by(TF) %>% summarise(mean = mean(importance)) %>% arrange(-mean) %>% top_n(200, mean)





mat <- filter(adjacencies, TF %in% Corr_TFs & target %in% Corr_target) %>% 
  pivot_wider(names_from = target, values_from = importance) %>%
  column_to_rownames("TF")

mat[is.na(mat)] <- 0

image <- 
  pheatmap(mat, scale = "column", fontsize = 10, 
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)),
         annotation_col = Ann,
         #annotation_colors = ann_colors, 
         treeheight_row = 0, treeheight_col = 0, 
         cutree_rows = 1, 
         cutree_cols = 30,
         fontsize_row = 7,
         fontsize_col = 4, 
         labels_col = NULL
         
         )


save2("Heatmap_GRN_Correlation_ConTargets", 3)


Ann <- 
CTPlus_Bind %>% 
  dplyr::select(1, events) %>% 
  arrange(-events) %>% 
  distinct() %>% 
  left_join(TFs) %>% 
  left_join(mean_logFC) %>% 
  as_tibble() %>% column_to_rownames("Gene") %>% dplyr::select(-TF)

Ann2 <- TopTable %>% 
  filter(Gene %notin% CTPlus_Bind$Gene) %>% 
  filter(is.na(DE.evts)) %>% group_by(Gene) %>% 
  summarise(mean_logFC = mean(Corr_FC)) %>% 
  left_join(dplyr::select(Corr_Bind, Gene, events)) %>% 
  distinct()

Ann3 <- rbind(rownames_to_column(Ann, "Gene"), Ann2)

Ann4 <- Marker_Bind %>% 
  filter(Gene %notin% Ann3$Gene) %>% 
  group_by(Gene) %>% 
  summarise(mean_logFC = mean(avg_logFC)) %>% 
  left_join(dplyr::select(Marker_Bind, Gene, events)) %>% 
  distinct()

Ann <- column_to_rownames(rbind(Ann3, Ann4), "Gene")


TF <- 
CTPlus_Bind %>% 
  dplyr::select(1, events) %>% 
  arrange(-events) %>% 
  distinct() %>% 
  left_join(TFs) %>% 
  left_join(mean_logFC) %>% 
  as_tibble() %>% filter(TF == "YES", events >=5) %>% pull(Gene)



Ann <- 
CTPlus_Bind %>% filter(Cluster == "Adipocyte") %>% 
  dplyr::select(1,3) %>% column_to_rownames("Gene")



mean_logFC <- 
  CTPlus_Bind %>% 
  group_by(Gene) %>% 
  summarise(mean_logFC = mean(avg_logFC)) %>% 
  arrange(-mean_logFC) %>% 
  left_join(dplyr::select(CTPlus_Bind, Gene, events)) %>% 
  distinct() %>% arrange(-events)


Ann <- 
  biglist %>% as_tibble() %>% 
  dplyr::select(Gene = "value") %>% 
  left_join(Sign, by = "Gene") %>% 
  dplyr::select(-Percentile, -Disc.Rate, -GMFC) %>%
  pivot_wider(names_from = Family, values_from = Family) %>% 
  unite(col = HR_Rec, 3:5, na.rm = TRUE, remove = T, sep = ", ") %>%
  dplyr::select(-`NA`) %>% left_join(Human_Secretome, by = "Gene") %>%
  left_join(mean_logFC, by = "Gene") %>% 
  left_join(dplyr::select(CTPlus_Bind, 1, 8), by = "Gene") %>% dplyr::select(Gene, mean_logFC) %>%
  distinct() %>%
  column_to_rownames("Gene")







# extract correaltion clusters from heatmap
A <- inner_join(TopTable,  # takes gene avrg_logFC and joins it with dendrogram ID
                rownames_to_column(as.data.frame(cutree(image$tree_col, k = 25)), "Gene"))
names(A)[13] <- "Corr.Cluster" # rename column name
#A <- arrange(A, Corr.Cluster) %>% mutate(CellType = name)

A %>% filter(Gene == "NR4A1")
A %>% filter(Corr.Cluster == "35") %>% pull(Gene) %>% unique()



for (i in 1:length(unique(A$Corr.Cluster))) {
  image <- reactomeBar (name = paste(i), x = filter(A, Corr.Cluster == i)$Gene,  15)
  print(image)
  }



reactomeBar (name = paste(25), x = filter(A, Corr.Cluster == 1)$Gene,  15)




TF_Regulon <- read_excel("utilities/TF_Regulon_Annotation_35.xlsx")

Ann <- 
A %>% dplyr::select(Gene, Corr.Cluster) %>%
  distinct() %>% left_join(TF_Regulon, by = "Corr.Cluster") %>% 
  arrange(Corr.Cluster) %>% dplyr::select(-Corr.Cluster) %>%
  column_to_rownames("Gene")


#ann_colors <- 



colors  <- c(pal_d3("category20")(20), brewer.pal(8, "Paired"), brewer.pal(8, "Set3")) %>% head(35)
  
names(colors) <- Ann$Label %>% unique()

ann_colors = list(
  #Time = c("white", "firebrick"),
  Label = colors
)


ann_colors = list(
  Time = c("white", "firebrick"),
  CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
  GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
)









Ann <- 
  filter(TopTable, Gene %in% Corr_target) %>% 
  select(Gene, CellType, Clust.Name) %>% 
  group_by(CellType) %>%
  pivot_wider(names_from = CellType, values_from = Clust.Name) %>% 
  column_to_rownames("Gene")


Ann <- 
  filter(TopTable, DE.evts > 0, Gene %in% Corr_target) %>% 
  dplyr::select(Gene, Category) %>% distinct() %>%
  #group_by(CellType) %>%
  #pivot_wider(names_from = CellType, values_from = Log_FC) %>% 
  column_to_rownames("Gene")



# Specify colors
ann_colors  = list(
  Adipocyte =    c("HER2_sig_Adp"       = "Blue", "LYPD_Adp" = "Yellow",  "Adipocyte_3" = "Red",  "Adipocyte_4" = "Green"),
  "LUM_HR-neg" = c("HER4+RUNX2_sig_LUn" = "Blue",     "ESR-EGF_LUn" = "Blue",   "IP_Metabolism_LUn" =  NA  ,   "LUM_HR-neg_5" =  NA ,       
  "Muscle_contraction_LUn" = NA,  "Heat_Stress_LUp" =  NA,      "Translation_LUn" =  "Orange" ,     "HER2+4_sig_LUn" =  "Blue",       "HER2_sig_LUn" = "Blue",         
  "LUM_HR-neg_3" =  NA )
  )






filter(Corr_Bind, Gene %in% Corr_target) %>% 
  select(Gene, CellType, Clust.Name)





FindMarkers(Sobj_integrated, test.use = )









