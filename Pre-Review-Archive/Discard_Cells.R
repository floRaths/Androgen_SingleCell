setwd("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/")

Discard <- readRDS("Output/MetaData_DiscardCells.rds")


Sobj_integrated <- AddMetaData(Sobj_integrated, select(Discard, Discard = "Discard"))
Sobj_integrated@meta.data$Discard       <- replace_na(Sobj_integrated@meta.data$Discard,       replace = "Keep")


Idents(Sobj_integrated) <- "Discard"
Sobj_integrated <- subset(Sobj_integrated, idents = "Keep")


Sobj_integrated <- DietSeurat(Sobj_integrated, counts = T, data = T, scale.data = F, assays = "RNA")
