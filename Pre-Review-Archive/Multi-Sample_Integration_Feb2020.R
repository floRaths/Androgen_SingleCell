library(Seurat)
library(tidyverse)

Sobj_list <- readRDS("/home/rathsf/scNuclei-Libraries/Analysis/Feb15_All_Nuclei_Sobj-List.rds")


Sobj_list <- SplitObject(Sobj_integrated, split.by = "Batch")
LUMA_Sample_list <- SplitObject(Cluster_Subset_list, split.by = "Sample")



### Initiating integration procedure

for (i in 1:length(Sobj_list)) {
  Sobj_list[[i]] <- NormalizeData        (Sobj_list[[i]], verbose = T)
  Sobj_list[[i]] <- FindVariableFeatures (Sobj_list[[i]], selection.method = "vst", nfeatures = 4000, verbose = T)
}


Sobj_anchors    <- FindIntegrationAnchors (object.list = Sobj_list,    dims = 1:50, anchor.features = 2000)
saveRDS(Sobj_anchors, "~/scNuclei-Libraries/Analysis/Reboot/Seurat_Objects/Sobj_anchors_2000_new.rds")

Sobj_integrated <- IntegrateData          (anchorset   = Sobj_anchors, dims = 1:50)
saveRDS(Sobj_integrated, "~/scNuclei-Libraries/Analysis/Reboot/Seurat_Objects/Sobj_integrated_2000_new.rds")

rm(Sobj_anchors, Sobj_list)

DefaultAssay (Sobj_integrated) <- "integrated"
Sobj_integrated <- ScaleData  (Sobj_integrated, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = T) ### with regression
#saveRDS(Sobj_integrated, "~/home/ubuntu/mnt/Seurat_Objects/Sobj_Integrated_Scaled.rds")

### continuing with regressed dataset
Sobj_integrated <- RunPCA        (Sobj_integrated, npcs = 50, verbose = FALSE)
Sobj_integrated <- RunUMAP       (Sobj_integrated, reduction  = "harmony", dims = 1:50)
Sobj_integrated <- FindNeighbors (Sobj_integrated, reduction  = "harmony", dims = 1:50)
Sobj_integrated <- FindClusters  (Sobj_integrated, resolution = c(0.1, 0.2, 0.5, 0.75, 1, 1.5))

saveRDS(Sobj_integrated, "~/scNuclei-Libraries/Analysis/Reboot/Seurat_Objects/Sobj_integrated_scaled_2000_new.rds")


#plot <- 
DimPlot(Sobj_integrated, 
        #cells = WhichCells(Sobj_integrated, idents = c("TM-3937", "CF-2797")),
        reduction = "umap", 
        group.by  = "integrated_snn_res.0.2",
        #group.by  = "CellType_All",
        split.by  = "Type",
        pt.size = 1, 
        label = T, 
        repel = T, 
        label.size = 5,
        #ncol = 5, 
        #cols = col
        ) 

#plot <- 
FeaturePlot(Sobj_integrated, reduction = "umap",
            #cells = WhichCells(Sobj_integrated, expression = SCORE_1 >= 0.25),
            features = c("KIT", "VWF", "PTPRC", "ESR1", "RELN", "TP63", "ANKRD30A", "COL6A3", "nCount_RNA"), 
            #features = c("MZB1", "MS4A1", "IGHM", "ITGAM", "CD44", "TCF7", "KRT15", "nCount_RNA"), 
            #features = c("ANKRD30A"), 
            #split.by = "Type",
            pt.size = 0.5, 
            order = T, 
            #min.cutoff = 0.5,
            cols = viridis::magma(n = 100), 
            #ncol = 4
            )



