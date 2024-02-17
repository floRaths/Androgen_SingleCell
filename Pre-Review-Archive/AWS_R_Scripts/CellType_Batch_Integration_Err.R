library(Seurat)
library(tidyverse)

Sobj_integrated               <- readRDS("~/mnt/Sobj_Integrated_Scaled_Proc.rds")
DefaultAssay(Sobj_integrated) <- "RNA"
Idents(Sobj_integrated)       <- "CellType"

Sobj_integrated <- DietSeurat(Sobj_integrated, counts = T, data = T, scale.data = F, assays = "RNA")
Sobj_integrated <- subset(Sobj_integrated, idents = c("Lum_HR-pos", "Lum_HR-neg", "Basal", "Endothelial", "Fibroblast", "Adipocyte"))
#Sobj_integrated <- subset(Sobj_integrated, idents = c("Fibroblast", "Adipocyte", "Myeloid", "Lymphoid"))

Subset_list     <- SplitObject(Sobj_integrated, split.by = "CellType")

for (i in 1:length(Subset_list)) {
  
  skip_to_next <- FALSE
  
  Subset_list[[i]] <- SplitObject(Subset_list[[i]], split.by = "Batch")
  
  
  for (s in 1:length(Subset_list[[i]])) {
    
    skip_to_next <- FALSE
    
    DefaultAssay                                           (Subset_list[[i]][[s]]) <- "RNA"
    
    tryCatch(Subset_list[[i]][[s]] <- NormalizeData        (Subset_list[[i]][[s]], verbose = T), error = function(e) { skip_to_next <<- TRUE})
    if(skip_to_next) { next }
    
    tryCatch(Subset_list[[i]][[s]] <- FindVariableFeatures (Subset_list[[i]][[s]], selection.method = "vst", nfeatures = 3000, verbose = T), error = function(e) { skip_to_next <<- TRUE})
    if(skip_to_next) { next }
    
  }
  
  
  tryCatch(Sobj_anchors     <- FindIntegrationAnchors (object.list = Subset_list[[i]],    dims = 1:50) , error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }
  tryCatch(Subset_list[[i]] <- IntegrateData          (anchorset   = Sobj_anchors, dims = 1:50)        , error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }
  rm(Sobj_anchors)
  
  DefaultAssay                      (Subset_list[[i]]) <- "integrated"
  Subset_list[[i]] <- ScaleData     (Subset_list[[i]], vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = T)
  Subset_list[[i]] <- RunPCA        (Subset_list[[i]], npcs = 50, verbose = FALSE)
  Subset_list[[i]] <- RunUMAP       (Subset_list[[i]], reduction  = "pca", dims = 1:50)
  Subset_list[[i]] <- FindNeighbors (Subset_list[[i]], reduction  = "pca", dims = 1:50)
  Subset_list[[i]] <- FindClusters  (Subset_list[[i]], resolution = c(0.2, 0.5, 1))


  saveRDS(Subset_list[[i]], paste0("~/mnt/CellType_Batch_Int_", levels(Sobj_integrated)[i], ".rds"))

}

saveRDS(Subset_list, "~/mnt/List_CellType_Batch_Integrated.rds")
