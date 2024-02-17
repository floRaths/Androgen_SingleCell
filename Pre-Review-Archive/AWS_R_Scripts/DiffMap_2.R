library(Seurat)
library(tidyverse)
library(destiny)


Sobj_integrated               <- readRDS("~/mnt/Sobj_Integrated_Scaled_Proc.rds")
DefaultAssay(Sobj_integrated) <- "RNA"
Idents(Sobj_integrated)       <- "CellType"

Sobj_integrated <- DietSeurat(Sobj_integrated, counts = T, data = T, scale.data = F, assays = "RNA")

#Sobj_integrated <- subset(Sobj_integrated, idents = c("Lum_HR-pos", "Lum_HR-neg", "Basal", "Endothelial", "Fibroblast", "Adipocyte", "Myeloid", "Lymphoid"))

CellType = "LUM_HR-neg"


Sobj_integrated <- subset(Sobj_integrated, idents = CellType)

Sobj_integrated <- NormalizeData        (Sobj_integrated, verbose = T, assay = "RNA")

Matrices <- as.matrix(Sobj_integrated@assays$RNA@data) %>% t

DiffMaps <- DiffusionMap(Matrices)

saveRDS(DiffMaps, paste0("~/mnt/DiffMap_", CellType, ".rds"))
