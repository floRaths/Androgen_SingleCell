library(destiny)

meta <- rownames_to_column(Sobj_integrated@meta.data, "X1") %>% mutate(Color = recode(Type, "TM" = "lightseagreen", "CF" = "mediumvioletred"))

celltype <- "LUM_HR-pos"

Sobj_integrated <- subset(Sobj_Full, idents = celltype)

# Reading in Destiny Output --------------------------------------------------------------------
Destiny <- readRDS(paste0("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Seurat_Objects/DataSubsets/Diffusion_Maps/CellType_All/CT_All_4k+TF_norm_",celltype,".csv.rds"))
Destiny <- readRDS(paste0("Output/Reboot/Destiny/Diffusion_LUM-pos_CF0404_CF3920_out.rds"))


# dataset(Destiny) <- dataset(Destiny) %>% left_join(meta, by = "X1")

plot.DiffusionMap(Destiny)
plot.DiffusionMap(Destiny, dims = 1:3, col_by = "Color")

# Adding Diffusion Embeddings to Subsetted Objects -----------------------
rownames(Destiny@eigenvectors) <- rownames(Sobj_integrated@reductions$umap)
Sobj_integrated[["diff"]] <- CreateDimReducObject(embeddings = Destiny@eigenvectors, key = "DC", assay = "RNA")



# Processing Subsets ------------------------------------------------------

Sobj_integrated <- Subsets$Basal_Myo
Sobj_integrated <- Subsets_A$Basal_Myo
Sobj_integrated <- Subsets_B$Basal_Myo
Idents(Sobj_integrated) <- "Treat"

Sobj_integrated <- subset(Sobj_integrated, idents = Idents(Sobj_integrated), downsample = 1600)


Sobj_integrated <- NormalizeData        (Sobj_integrated, normalization.method = "LogNormalize", 
                                         scale.factor = 10000, assay = "RNA")
Sobj_integrated <- FindVariableFeatures (Sobj_integrated, selection.method = "vst", nfeatures = 10000, assay = "RNA")


Sobj_integrated <- ScaleData(Sobj_integrated, vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "integrated")
## pbmc <- ScaleData(pbmc) would make it faster by just scaling variable genes

Sobj_integrated <- RunPCA(Sobj_integrated, features = VariableFeatures(object = Sobj_integrated), 
                          npcs = 50, assay = "integrated")

DefaultAssay(Sobj_integrated) <- "RNA"
Sobj_integrated <- FindNeighbors (Sobj_integrated, dims = 1:20, reduction = "diff", assay = "RNA")
Sobj_integrated <- FindClusters  (Sobj_integrated, resolution = c(0.15), assay = "RNA")
Sobj_integrated <- RunUMAP       (Sobj_integrated, dims = 1:50, assay = "integrated")



#plot1 <- 
DimPlot(Sobj_integrated, reduction = "diff", 
        group.by = "Sample",
        split.by = "Type", 
        pt.size = 1, dims = 2:3)

plot2 <- 
  DimPlot(Sobj_integrated, reduction = "diff", 
          group.by = "Type",
          split.by = "Type", 
          pt.size = 1, dims = 2:3)

image <- CombinePlots(list(plot1, plot2), ncol = 1)

save("Destiny_Basal", 1)  



image = CombinePlots(plots = list(A, B), ncol = 1)
ggsave(file = "./Figrues/10.24_Nuclei-Alignment_MAST_Responses/test.eps", plot = image, width=16, height=9, )








DimPlot(Sobj_integrated, reduction =  "diff", dims = 1:2, group.by = "Type", split.by = "Sample", pt.size = 0.5, ncol = 5)

DimPlot(Sobj_integrated, reduction =  "umap", 
        #dims = 2:3, 
        group.by = "RNA_snn_res.0.15", split.by = "Type",
        pt.size = 1, cols = col)



c("CUX2_1", "ZNF689_1", "PGR_1", "BATF_1", "ESR1_1", "CERS4_1")


plot6 <- 
  FeaturePlot(Sobj_integrated,
              reduction = "diff", 
              features = c("CERS4_1"), pt.size = 1, 
              order = T, 
              dims = 2:3,
              #split.by = "Type",
              cols = viridis::magma(n = 100)) + 
  xlim (c(-0.025, 0.055)) + 
  ylim (c(-0.035, 0.06))



plot1 + plot2 + plot3 + plot4 + plot5  + plot6 + patchwork::plot_layout(ncol = 2)



plot2 <- 
  FeaturePlot(Sobj_integrated, 
              reduction = "diff", 
              features = c("NEBL"), pt.size = 1, 
              sort.cell = F, 
              dims = 2:3,
              #split.by = "Type",
              cols = viridis::magma(n = 100)) + 
  xlim (c(-0.05, 0.05)) + 
  ylim (c(-0.05, 0.05))



plot3 <- 
FeaturePlot(Sobj_integrated, 
            reduction = "diff", 
            features = c("ERBB4"), pt.size = 1, 
            sort.cell = F, 
            dims = 1:2,
            #split.by = "Type",
            cols = viridis::magma(n = 100)) + 
  xlim (c(-0.05, 0.05)) + 
  ylim (c(-0.05, 0.05))

plot4 <- 
  FeaturePlot(Sobj_integrated, 
              reduction = "diff", 
              features = c("ERBB4"), pt.size = 1, 
              sort.cell = F, 
              dims = 2:3,
              #split.by = "Type",
              cols = viridis::magma(n = 100)) + 
  xlim (c(-0.05, 0.05)) + 
  ylim (c(-0.05, 0.05))





image <- CombinePlots(list(plot1, plot2, plot3, plot4), ncol = 2)
save("Destiny_Basal_ERBB4_NEBL_sort=F", 1)  















Diff <- Destiny@eigenvectors[,1:3]  %>% 
  as.data.frame() %>% rownames_to_column("Cell") %>% 
  left_join(rownames_to_column(Sobj_Full@meta.data, "Cell"), by = "Cell")




x <- Diff$DC1
y <- Diff$DC2
z <- Diff$DC3




# x, y and z coordinates
x <- sep.l <- iris$Sepal.Length
y <- pet.l <- iris$Petal.Length
z <- sep.w <- iris$Sepal.Width

scatter3D(x, y, z, colvar = as.numeric(Diff$Type), 
          col = c("#D95F02", "#7570B3"), pch = 19, alpha  = 0.5, theta = 20, phi = 25, cex = 0.5,
          colkey = list(at = c(2, 3), labels = c("setosa", "versicolor")))




scatter3D(x, y, z, bty = "g", pch = 18, 
          col.var = as.integer(iris$Species), 
          col = c("#1B9E77", "#D95F02", "#7570B3"),
          pch = 18, ticktype = "detailed",
          colkey = list(at = c(2, 3, 4), side = 1, 
                        addlines = TRUE, length = 0.5, width = 0.5,
                        labels = c("setosa", "versicolor", "virginica")) )








library(plyr)
library(grid)



levels <- levels(Sobj_Full)

for (i in 1:length(levels)) {

  celltype <- "DC"
  
  Sobj_integrated <- subset(Sobj_Full, idents = levels[i])

# Reading in Destiny Output --------------------------------------------------------------------
Destiny <- readRDS(paste0("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Seurat_Objects/DataSubsets/Diffusion_Maps/CellType_All/CT_All_4k+TF_norm_",celltype,".csv.rds"))

dev.off()
plot(Destiny)

rm(image)

grid.echo()

image <- grid.grab()


save(paste0("DiffMap_", celltype), 1)

genes <- colnames(Destiny@data_env$data)[-1]


x <- ldply(1:length(genes), function(s)c(genes[s], 
                         cor(Destiny$DC1, pull(dplyr::select(Destiny@data_env$data, s+1),1), method = "pearson"), 
                         cor(Destiny$DC2, pull(dplyr::select(Destiny@data_env$data, s+1),1), method = "pearson"),
                         cor(Destiny$DC3, pull(dplyr::select(Destiny@data_env$data, s+1),1), method = "pearson")))

x <- x %>% dplyr::select(Gene = 1, DC1 = 2, DC2 = 3, DC3 = 4) #%>% top_n(25, DC3) %>% arrange(DC3)

saveRDS(x, paste0("Output/Destiny/Correlations/Correlation_", celltype, ".rds"))

}




celltype <- "LUM_HR-pos"



query <- readRDS(paste0("Output/Destiny/Correlations/Correlation_", celltype, ".rds")) %>% top_n(100, DC2) %>% pull(Gene) %>% sort()
query2 <- `Correlation_LUM_HR-pos` %>% top_n(-150, DC2) %>% pull(Gene) %>% sort()
#query <- Correlation_Macrophage %>% top_n(-100, DC1) %>% pull(Gene)

reactomeBar(name = paste0(celltype, " - DC3 correlating"), query, 15)


p1 <- reactomeBar(name = paste0(celltype, " - DC2 correlating"), query, 15)
p2 <- reactomeBar(name = paste0(celltype, " - DC2 anti-correlating"), query2, 15)



query <- Marker_Bind %>% ungroup() %>% filter(CellType == "Macrophage") %>% top_n(-50, avg_logFC) %>% pull(Gene)












