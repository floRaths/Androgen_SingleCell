library(Seurat)
#library(tidyverse)

Sobj_integrated <- readRDS("/home/rathsf/scNuclei-Libraries/Analysis/Top_Samples/Nuclei_Top_Samples.rds")
DefaultAssay(Sobj_integrated) <- "RNA"

#Idents(Sobj_integrated) <- "Batch"
#Sobj_integrated <- subset(Sobj_integrated, idents = "Epithelial")

#Sobj_integrated <- DietSeurat(Sobj_integrated, counts = T, data = T, scale.data = F, assays = "RNA")

Batch_list     <- SplitObject(Sobj_integrated, split.by = "Batch")

for (i in 1:length(Batch_list)) {

    DefaultAssay                            (Batch_list[[i]]) <- "RNA"
    Batch_list[[i]] <- NormalizeData        (Batch_list[[i]], verbose = T)
    Batch_list[[i]] <- FindVariableFeatures (Batch_list[[i]], selection.method = "vst", nfeatures = 3000, verbose = T)
    }
  

Sobj_anchors          <- FindIntegrationAnchors (object.list = Batch_list,    dims = 1:20)
Sobj_integrated       <- IntegrateData          (anchorset   = Sobj_anchors,  new.assay.name = "intgr", dims = 1:50)

saveRDS(Sobj_anchors, "~/mnt/Seurat_Objects/Sobj_BatchAnchors_2000.rds")

rm(Sobj_anchors)

  
DefaultAssay                     (Sobj_integrated) <- "intgr"
Sobj_integrated <- ScaleData     (Sobj_integrated, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = T)
Sobj_integrated <- RunPCA        (Sobj_integrated, npcs = 50, verbose = FALSE)
Sobj_integrated <- RunUMAP       (Sobj_integrated, reduction  = "pca", dims = 1:50)
Sobj_integrated <- FindNeighbors (Sobj_integrated, reduction  = "pca", dims = 1:50)
Sobj_integrated <- FindClusters  (Sobj_integrated, resolution = c(0.2, 0.5, 0.75, 1, 1.5))


saveRDS(Sobj_integrated, "~/mnt/Seurat_Objects/Sobj_BatchIntegrated_2000.rds")
