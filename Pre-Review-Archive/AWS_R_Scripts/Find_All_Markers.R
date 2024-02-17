library(Seurat)
library(tidyverse)

Sobj_integrated               <- readRDS("~/mnt/Sobj_Integrated_Scaled_Proc.rds")
Idents(Sobj_integrated)       <- "CellType"
DefaultAssay(Sobj_integrated) <- "RNA"


AllMarks <- FindAllMarkers(Sobj_integrated, 
                           test.use = "MAST", 
                           assay    = "RNA", 
                           slot     = "data", 
                           min.pct  = 0.35, 
                           verbose  = T)


saveRDS(AllMarks, "~/mnt/AllMarks_Celltype_MAST.rds")
