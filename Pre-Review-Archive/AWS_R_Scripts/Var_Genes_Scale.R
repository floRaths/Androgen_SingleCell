library(Seurat)
library(tidyverse)
#library(future)

#options(future.globals.maxSize= 3670016000)
#plan("multiprocess", workers = 32)

Sobj_integrated <- readRDS("~/mnt/Sobj_Integrated_Scaled_Proc.rds")
DefaultAssay(Sobj_integrated) <- "RNA"

Sobj_integrated <- DietSeurat(Sobj_integrated, counts = T, data = T, scale.data = F, assays = "RNA")

Sobj_integrated <- NormalizeData        (Sobj_integrated, verbose = T)
Sobj_integrated <- FindVariableFeatures (Sobj_integrated, selection.method = "vst", nfeatures = 8000, verbose = T)

all.genes <- rownames(Sobj_integrated)
Sobj_integrated <- ScaleData(Sobj_integrated, 
                             #features = all.genes, 
                             vars.to.regress = c("nCount_RNA", "percent.mt"), 
                             verbose = T)

Sobj_integrated <- RunPCA        (Sobj_integrated, npcs = 50, verbose = FALSE)
Sobj_integrated <- RunUMAP       (Sobj_integrated, reduction  = "pca", dims = 1:50)
Sobj_integrated <- FindNeighbors (Sobj_integrated, reduction  = "pca", dims = 1:50)
Sobj_integrated <- FindClusters  (Sobj_integrated, resolution = c(0.1, 0.2, 0.5, 1))


saveRDS(Sobj_integrated, "~/mnt/VarGenes_Scaled.rds")
