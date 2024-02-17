library(DoubletFinder)

Idents(Sobj_integrated) <- "Sample"
#CF3920 <- subset(Sobj_integrated, idents = c("CF_3920"))
seu_kidney <- subset(Sobj_integrated, idents = c("TM-8249"))
seu_kidney <- DietSeurat(seu_kidney, counts = T, data = T, scale.data = F, assays = "RNA")
#CF3920 <- readRDS("Seurat_Objects/Doublet_Alignment.rds")



DefaultAssay(seu_kidney) <- "RNA"

seu_kidney <- NormalizeData(seu_kidney)
seu_kidney <- ScaleData(seu_kidney)
seu_kidney <- FindVariableFeatures(seu_kidney, selection.method = "vst", nfeatures = 2000)
seu_kidney <- RunPCA(seu_kidney)
seu_kidney <- RunUMAP(seu_kidney, dims = 1:50)
seu_kidney <- FindNeighbors (seu_kidney, reduction  = "pca", dims = 1:50)
seu_kidney <- FindClusters  (seu_kidney, resolution = c(0.2))

DimPlot(seu_kidney, pt.size = 0.5, label = T)



## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:20, sct = FALSE)
sweep.stats_kidney    <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney          <- find.pK(sweep.stats_kidney)

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(CF3920, PCs = 1:10, sct = FALSE)
gt.calls <- CF3920@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"]
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
bcmvn_kidney <- find.pK(sweep.stats_kidney)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seu_kidney@meta.data$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi       <- round(0.15*length(colnames(seu_kidney)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj   <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:50, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)


DimPlot(seu_kidney, pt.size = 1, group.by = "DF.classifications_0.25_0.005_1545", label = T)





TM9469meta.adj <- 
  separate(rownames_to_column(subset(Sobj_integrated, idents = "TM_9469")@meta.data), rowname, into = c("A", "B")) %>% 
  mutate(B = 11) %>% 
  unite(col = "ROW", c("A", "B")) %>% 
  column_to_rownames("ROW")



TM <- select(TM9469meta.adj, Class = "DF.classifications_0.25_0.18_556")
CF <- select(CF3920meta.adj, Class = "DF.classifications_0.25_0.01_944")

MetaAdj <- rbind(TM9469meta.adj, CF3920meta.adj)

Test <- AddMetaData(Sobj_integrated, metadata = rbind(TM, CF))
Test <- AddMetaData(Test, metadata = TM)


plot <- DimPlot(Test, pt.size = 0.5, 
        #group.by = "Sample", 
        #split.by = "Class", 
        ncol = 2)



Sobj_integrated <- AddMetaData(Sobj_integrated, metadata = MetaAdj)

Sobj_integrated@meta.data$Group_New   <- replace_na(Sobj_integrated@meta.data$Group_New,   replace = "Discard")
Sobj_integrated@meta.data$DoubletSelect <- replace_na(Sobj_integrated@meta.data$DoubletSelect, replace = "Singlet")
Sobj_integrated@meta.data$DoubletComb   <- replace_na(Sobj_integrated@meta.data$DoubletComb,   replace = "Singlet")

Idents(Sobj_integrated) <- "DoubletComb"

Sobj_integrated <- subset(Sobj_integrated, idents = "Singlet")

FeaturePlot(Test, features = "nCount_RNA", pt.size = 1)


DimPlot(Sobj_integrated, split.by = "DoubletComb", pt.size = .5)


select.cells <- CellSelector(plot = plot)

Idents(Sobj_integrated, cells = select.cells) <- "NewCells"

DimPlot(Sobj_integrated, pt.size = .5)

DimPlot(Sobj_integrated, group.by = "DoubletComb", pt.size = .5)

Sobj_integrated@meta.data$DoubletExt <- recode(Sobj_integrated@meta.data$DoubletExt, NewCells = "Doublet",
                                             CF_3920 = "Singlet", 
                                             TM_9469 = "Singlet", 
                                             TM_1956 = "Singlet",
                                             CF_19301 = "Singlet")



Subset <- subset(Sobj_integrated, idents = c("CF_3920", "TM_9469"))

WhichCells(Sobj_integrated, idents = "Doublet")


Sobj_integrated@meta.data$DoubletComb <- Idents(Sobj_integrated)


Assign <- select(MetaAdj, DoubletFind = "Doublet", DoubletSelect = "DoubletExt", DoubletComb = "DoubletComb")

saveRDS(Assign, "Output/Doublet_Assign.rds")










Sobj_integrated <- DietSeurat(Sobj_integrated, counts = T, data = T, scale.data = F, assays = "RNA")

Idents(Sobj_integrated) <- "Sample"
DefaultAssay(Sobj_integrated) <- "RNA"


Samples           <- levels(Sobj_integrated)
meta.data.list    <- vector("list", length = length(levels(Sobj_integrated)))

for (i in seq_along(levels(Sobj_integrated))) {
  
seu_kidney <- subset(Sobj_integrated, idents = Samples[i])
seu_kidney <- DietSeurat(seu_kidney, counts = T, data = T, scale.data = F, assays = "RNA")
#CF3920 <- readRDS("Seurat_Objects/Doublet_Alignment.rds")


seu_kidney <- NormalizeData(seu_kidney)
seu_kidney <- ScaleData(seu_kidney)
seu_kidney <- FindVariableFeatures(seu_kidney, selection.method = "vst", nfeatures = 2000)
seu_kidney <- RunPCA(seu_kidney)
seu_kidney <- RunUMAP(seu_kidney, dims = 1:50)
seu_kidney <- FindNeighbors (seu_kidney, reduction  = "pca", dims = 1:50)
seu_kidney <- FindClusters  (seu_kidney, resolution = c(0.2))

DimPlot(seu_kidney, pt.size = 0.5, group.by = "Sample", label = T)



## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:25, sct = FALSE)
sweep.stats_kidney    <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney          <- find.pK(sweep.stats_kidney)


## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
#homotypic.prop <- modelHomotypic(seu_kidney@meta.data$integrated_snn_res.0.1)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi       <- round(0.15*length(colnames(seu_kidney)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
#nExp_poi.adj   <- round(nExp_poi*(1-homotypic.prop))
pK <- bcmvn_kidney %>% arrange(-BCmetric) %>% head(1) %>% pull(pK) %>% as.numeric()
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:50, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

meta.data.list[[i]] <- seu_kidney@meta.data

}

DimPlot(seu_kidney, pt.size = 1, group.by = "DF.classifications_0.25_0.04_1626", label = T)







for (i in 1:length(meta.data.list)) {
  meta.data.list[[i]] <- meta.data.list[[i]] %>% select(14, Doublet = 19)
}

Bind <- do.call(rbind.data.frame, meta.data.list)



Scrublet <- 
Scrublet %>% 
  separate("ID", into = c("A", "B")) %>% 
  add_column(Suffix = "1") %>% 
  unite("X", c("Suffix", "B"), sep = "_") %>% 
  unite("ID", c("A", "X"), sep = "-") %>% 
  add_column("Scrublet" = "Singlet")

Meta2 <- 
Meta %>% 
  left_join(Scrublet, by = "ID") %>% 
  mutate(Scrublet3x = replace_na(Scrublet, "Doublet")) %>% 
  column_to_rownames("ID")
