library(Seurat)
library(tidyverse)

Sobj_integrated <- readRDS("~/mnt/Sobj_Integrated_Scaled_Proc.rds")
MetaAdj         <- readRDS("~/mnt/Doublet_Assign.rds")


Sobj_integrated <- AddMetaData(Sobj_integrated, metadata = MetaAdj)

Sobj_integrated@meta.data$DoubletFind   <- replace_na(Sobj_integrated@meta.data$DoubletFind,   replace = "Singlet")
Sobj_integrated@meta.data$DoubletSelect <- replace_na(Sobj_integrated@meta.data$DoubletSelect, replace = "Singlet")
Sobj_integrated@meta.data$DoubletComb   <- replace_na(Sobj_integrated@meta.data$DoubletComb,   replace = "Singlet")

Idents(Sobj_integrated) <- "DoubletComb"
DefaultAssay(Sobj_integrated) <- "RNA"



Sobj_integrated <- subset(Sobj_integrated, idents = "Singlet")

Sobj_integrated <- DietSeurat(Sobj_integrated, counts = T, data = T, scale.data = F, assays = "RNA")

Batch_list      <- SplitObject(Sobj_integrated, split.by = "Batch")



for (i in 1:length(Batch_list)) {

    DefaultAssay                            (Batch_list[[i]]) <- "RNA"
    Batch_list[[i]] <- NormalizeData        (Batch_list[[i]], verbose = T)
    Batch_list[[i]] <- FindVariableFeatures (Batch_list[[i]], selection.method = "vst", nfeatures = 3000, verbose = T)

  }
  
Sobj_anchors     <- FindIntegrationAnchors (object.list = Batch_list,    dims = 1:50)
Subset_list      <- IntegrateData          (anchorset   = Sobj_anchors, dims = 1:50)
rm(Sobj_anchors)

  
DefaultAssay                 (Subset_list) <- "integrated"
Subset_list <- ScaleData     (Subset_list, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = T)
Subset_list <- RunPCA        (Subset_list, npcs = 50, verbose = FALSE)
Subset_list <- RunUMAP       (Subset_list, reduction  = "pca", dims = 1:50)
Subset_list <- FindNeighbors (Subset_list, reduction  = "pca", dims = 1:50)
Subset_list <- FindClusters  (Subset_list, resolution = c(0.2, 0.5, 1))


saveRDS(Subset_list, "~/mnt/Nuclei_All_Samples_Batch_Integrated_DoubRM.rds")
