library(Seurat)
library(tidyverse)
library(future)

plan("multiprocess", workers = 32)


Sobj_integrated <- readRDS("~/mnt/Sobj_Integrated_Scaled_Proc.rds")
DefaultAssay(Sobj_integrated) <- "RNA"

Sobj_integrated <- DietSeurat(Sobj_integrated, counts = T, data = T, scale.data = F, assays = "RNA")

Subset_list     <- SplitObject(Sobj_integrated, split.by = "Type")

for (i in 1:length(Subset_list)) {
  
  Subset_list[[i]] <- SplitObject(Subset_list[[i]], split.by = "Batch")
  
  for (s in 1:length(Subset_list[[i]])) {
    DefaultAssay                                  (Subset_list[[i]][[s]]) <- "RNA"
    Subset_list[[i]][[s]] <- NormalizeData        (Subset_list[[i]][[s]], verbose = T)
    Subset_list[[i]][[s]] <- FindVariableFeatures (Subset_list[[i]][[s]], selection.method = "vst", nfeatures = 3000, verbose = T)
  }
  
  Sobj_anchors    <- FindIntegrationAnchors (object.list = Subset_list[[i]],    dims = 1:50)
  Subset_list[[i]] <- IntegrateData          (anchorset   = Sobj_anchors, dims = 1:50)
  rm(Sobj_anchors)
  
  DefaultAssay                      (Subset_list[[i]]) <- "integrated"
  Subset_list[[i]] <- ScaleData     (Subset_list[[i]], vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = T)
  Subset_list[[i]] <- RunPCA        (Subset_list[[i]], npcs = 50, verbose = FALSE)
  Subset_list[[i]] <- RunUMAP       (Subset_list[[i]], reduction  = "pca", dims = 1:50)
  Subset_list[[i]] <- FindNeighbors (Subset_list[[i]], reduction  = "pca", dims = 1:50)
  Subset_list[[i]] <- FindClusters  (Subset_list[[i]], resolution = c(0.2, 0.5, 1))
}

saveRDS(Subset_list, "~/mnt/Subset_List_Batch_Integrated.rds")


Subset_Integration.R
Subset_Integration.R
