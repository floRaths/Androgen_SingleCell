
# Extracting Normalized Count Matrices ------------------------------------

Pat_A <- subset(Hormone.Integrated, idents = "Patient_A")
Pat_B <- subset(Hormone.Integrated, idents = "Patient_B")

Idents(Sobj_integrated) <- "Cluster"
Idents <- levels(Sobj_integrated)

Matrices <- vector("list", length = length(Idents))

for (i in 1:length(Idents)) {
  Matrices[[i]] <- as.matrix(subset(Sobj_integrated, idents = Idents[i])@assays$RNA@data) %>% t
}

names(Matrices) <- Idents




# Running Diffusion -------------------------------------------------------

library(destiny)

Pat_A_Matrices <- readRDS("/home/rathsf/Organoid_Project/Analysis/Destiny_Output/Per_Patient/Pat_A_matrices.rds")

Idents <- c("Basal", "Basal_Myo", "Luminal_HR-neg", "Stroma", "Basal/Luminal", "Low.UMI", "Basal.CC", "Luminal_HR-pos")
DiffMaps <- vector("list", length = length(Idents))

for (i in 1:length(Idents)) {
  DiffMaps[[i]] <- DiffusionMap(Matrices[[i]])
}

saveRDS(Pat_A_DiffMaps, "/home/rathsf/Organoid_Project/Analysis/Destiny_Output/Per_Patient/Pat_A_DiffMaps.rds")

names(Pat_A_DiffMaps) <- Idents

# Creating a List of subsetted Seurat Objects -----------------------------


Sobj <- Pat_A
Idents(Sobj) <- "Cluster"

### Preparation of our file lists that will contain our seurat objects
Names   <- levels(Sobj)
Subsets <- vector("list", length = length(Names))

for (i in 1:length(Names)) {
  Subsets[[i]] <- subset(Sobj, idents = Names[i])
  
}

names(Subsets) <- Names


# Reading in Destiny Output --------------------------------------------------------------------

Names   <- levels(Sobj)
Destiny <- vector("list", length = length(Names))

for (i in 1:8) {
  Destiny[[i]] <- readRDS(paste0("~/Box Sync/Knott_Lab/Flo/Projects/Organoids/Analysis/2019.02.17_Oganoid_Hormone_TimeCourse/Outputs/V3/Destiny_Output/Destiny_", Names[i],".rds"))
}

names(Destiny) <- Names



# Adding Diffusion Embeddings to Subsetted Objects -----------------------

Destiny <- AdipoDiff
Subsets <- Adipo

for (i in 1:8) {
  rownames(Destiny[[i]]@eigenvectors) <- rownames(Subsets[[i]]@reductions$pca)
  Subsets[[i]][["diff"]] <- CreateDimReducObject(embeddings = Destiny[[i]]@eigenvectors, key = "DC", assay = "RNA")
}


rownames(Destiny@eigenvectors) <- rownames(Subsets@reductions$umap)
Subsets[["diff"]] <- CreateDimReducObject(embeddings = Destiny@eigenvectors, key = "DC", assay = "RNA")


# Processing Subsets ------------------------------------------------------

Sobj <- Subsets$Basal_Myo
Sobj <- Subsets_A$Basal_Myo
Sobj <- Subsets_B$Basal_Myo
Idents(Sobj) <- "Treat"

Sobj <- subset(Sobj, idents = Idents(Sobj), downsample = 1600)


Sobj <- NormalizeData        (Sobj, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
Sobj <- FindVariableFeatures (Sobj, selection.method = "vst", nfeatures = 10000, assay = "RNA")


Sobj <- ScaleData(Sobj, vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "integrated")
## pbmc <- ScaleData(pbmc) would make it faster by just scaling variable genes

Sobj <- RunPCA(Sobj, features = VariableFeatures(object = Sobj), npcs = 50, assay = "integrated")

DefaultAssay(Sobj) <- "RNA"
Sobj <- FindNeighbors (Sobj, dims = 1:20, reduction = "diff", assay = "RNA")
Sobj <- FindClusters  (Sobj, resolution = c(0.15), assay = "RNA")
Sobj <- RunUMAP       (Sobj, dims = 1:50, assay = "integrated")


DimPlot(Sobj, reduction = "diff", split.by = "Treat", group.by = "Sample", pt.size = 1)
DimPlot(Sobj, reduction = "diff", split.by = "Sample", pt.size = 1, dims = 1:2)
DimPlot(Sobj, reduction = "diff", split.by = "Treat", pt.size = 1, dims = 1:2)



DimPlot(Sobj, reduction = "diff", dims = 1:2, group.by = "RNA_snn_res.0.15", pt.size = 1)
DimPlot(Sobj, reduction = "diff", dims = 2:3, group.by = "RNA_snn_res.0.15", pt.size = 1)
DimPlot(Sobj, reduction = "diff", dims = 3:4, group.by = "RNA_snn_res.0.15", pt.size = 1)


DimPlot(Sobj, reduction = "diff", dims = 1:2, split.by = "Treat", group.by = "RNA_snn_res.0.15", pt.size = 1)
DimPlot(Sobj, reduction = "diff", dims = 2:3, split.by = "Treat", group.by = "RNA_snn_res.0.15", pt.size = 1)
DimPlot(Sobj, reduction = "diff", dims = 1:2, split.by = "Treat", group.by = "Sample", pt.size = 1)
DimPlot(Sobj, reduction = "diff", dims = 2:3, split.by = "Treat", group.by = "RNA_snn_res.0.15", pt.size = 1)
DimPlot(Sobj, reduction = "diff", dims = 3:4, split.by = "Treat", group.by = "RNA_snn_res.0.15", pt.size = 1)
DimPlot(Sobj, reduction = "diff", dims = 4:5, split.by = "Treat", group.by = "RNA_snn_res.0.15", pt.size = 1)





image = CombinePlots(plots = list(A, B), ncol = 1)
ggsave(file = "./Figrues/10.24_Nuclei-Alignment_MAST_Responses/test.eps", plot = image, width=16, height=9, )








DimPlot(Sobj, reduction =  "diff", dims = 2:3, group.by = "RNA_snn_res.0.15", pt.size = 1, cols = col)

DimPlot(Sobj, reduction =  "umap", 
        #dims = 2:3, 
        group.by = "RNA_snn_res.0.15", split.by = "Type",
        pt.size = 1, cols = col)




FeaturePlot(Sobj_integrated, reduction = "diff", 
            features = "PIP", pt.size = 1, 
            sort.cell = F, 
            dims = 1:2,
            #split.by = "Type",
            cols = viridis::magma(n = 100)) + 
  xlim (c(-0.02, 0.06)) + 
  ylim (c(-0.05, 0.06))









Idents(Sobj_integrated) <- "CellType"

Sobj_integrated <- subset(Sobj_integrated, idents = "Lum_HR-pos")



