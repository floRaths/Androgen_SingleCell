library(Seurat)
library(tidyverse)
library(destiny)


Sobj_integrated               <- readRDS("home/rathsf/scNuclei-Libraries/Analysis/Top_Samples/Nuclei_Top_Samples.rds")
DefaultAssay(Sobj_integrated) <- "RNA"
Idents(Sobj_integrated)       <- "CellType_All"

Sobj_integrated <- DietSeurat(Sobj_integrated, counts = T, data = T, scale.data = F, assays = "RNA")

#Sobj_integrated <- subset(Sobj_integrated, idents = c("Lum_HR-pos", "Lum_HR-neg", "Basal", "Endothelial", "Fibroblast", "Adipocyte", "Myeloid", "Lymphoid"))

CellType = "LUM_HR-pos"


Sobj_integrated <- subset(Sobj_integrated, idents = CellType)

Sobj_integrated <- NormalizeData        (Sobj_integrated, verbose = T, assay = "RNA")

Matrices <- t(as.matrix(Sobj_integrated@assays$integrated@data))

DiffMaps <- DiffusionMap(Matrices)

saveRDS(DiffMaps, paste0("home/rathsf/scNuclei-Libraries/Analysis/Top_Samples/DiffMap_", CellType, ".rds"))
