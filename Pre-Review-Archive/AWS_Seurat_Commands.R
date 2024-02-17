library(Seurat)
library(tidyverse)
library(future)

plan("multiprocess", workers = 6)
options(future.globals.maxSize= 2621440000) #850(MB)*1024^2 = globals VALUE


Sobj_list <- readRDS("~/mnt/Sobj_List.rds")
Sobj_integrated <- readRDS("~/mnt/Sobj_Integrated_Scaled_Proc.rds")

MetaData <- readRDS("~/mnt/MetaData_After_Trimming.rds")
Sobj_integrated <- AddMetaData(Sobj_integrated, select(MetaData, Cluster = "Cluster", Group = "Group", CellType = "CellType", CellType_Plus = "CellType_Plus"))


### assigning Patient_ID/replicate_ID to the Data
Idents(Sobj_integrated) <- "orig.ident"
new.cluster.ids <- substring(levels(Sobj_integrated), 5, 100) #remove the "Nuc-" prefix and only keep saple name
names(new.cluster.ids) <- levels(Sobj_integrated)
Sobj_integrated <- RenameIdents(Sobj_integrated, new.cluster.ids)
Sobj_integrated@meta.data$Sample <- Idents(Sobj_integrated)


Idents(Sobj_integrated) <- "orig.ident"
new.cluster.ids <- substr(levels(Sobj_integrated), 5, nchar(levels(Sobj_integrated))-8)
names(new.cluster.ids) <- levels(Sobj_integrated)
Sobj_integrated <- RenameIdents(Sobj_integrated, new.cluster.ids)
Sobj_integrated@meta.data$Sample <- Idents(Sobj_integrated)



### assigning Treatment groups to the Data
Idents(Sobj_integrated) <- "orig.ident"
new.cluster.ids <- substring(levels(Sobj_integrated), 5, 6) #only keep CF or TM
names(new.cluster.ids) <- levels(Sobj_integrated)
Sobj_integrated <- RenameIdents(Sobj_integrated, new.cluster.ids)
Sobj_integrated@meta.data$Type <- Idents(Sobj_integrated)


### assigning Batches to the Data
Idents(Sobj_integrated) <- "integrated_snn_res.0.2"
new.cluster.ids <- c("CD4_T-Cells",        #0
                     "Macrophages",        #1
                     "Monocyte",      #2
                     "CD8_T-Cells",     #3
                     "NK_Cells",     #4
                     "Dendritic Cell",     #5
                     "Maybe DC",     #6
                     "Unknown_Vasc",     #7
                     "B-Cell (MZ)",     #8
                     "Unknown_Neuro1",     #9
                     "Unknown_Neuro2",     #10
                     "B-Cell",      #11
                     "HSC",     #12
                     "Undown_Epith_Doubl?",     #13
                     "MAST_Cell",     #14
                     "Unknown"     #15
                     
)

names(new.cluster.ids) <- levels(Sobj_integrated)
Sobj_integrated <- RenameIdents(Sobj_integrated, new.cluster.ids)
Sobj_integrated@meta.data$Cluster <- Idents(Sobj_integrated)


### assigning prep method groups to the Data
Idents(Sobj_integrated) <- "orig.ident"
new.cluster.ids <- c("Sort",     #1
                     "Sort",     #2
                     "Sort",     #3
                     "Sort",     #4
                     "Sort",     #5
                     "Sort",     #6
                     "Count",     #7
                     "Sort",     #8
                     "Sort",     #9
                     "Sort",     #10
                     "Sort",     #11
                     "Sort",     #12
                     "Sort",     #13
                     "Sort",     #14
                     "Sort",     #15
                     "Sort",     #16
                     "Sort",     #17
                     "Sort",     #18
                     "Sort",     #19
                     "Sort",     #20
                     "Count",     #21
                     "Sort"      #22
)


names(new.cluster.ids) <- levels(Sobj_integrated)
Sobj_integrated <- RenameIdents(Sobj_integrated, new.cluster.ids)
Sobj_integrated@meta.data$Prep <- Idents(Sobj_integrated)



### assigning prep method groups to the Data
Idents(Sobj_integrated) <- "orig.ident"
new.cluster.ids <- c("Batch_5",     #1
                     "Batch_2",     #2
                     "Batch_7",     #3
                     "Batch_6",     #4
                     "Batch_4",     #5
                     "Batch_4",     #6
                     "Batch_1",     #7
                     "Batch_7",     #8
                     "Batch_6",     #9
                     "Batch_3",     #10
                     "Batch_5",     #11
                     "Batch_7",     #12
                     "Batch_2",     #13
                     "Batch_4",     #14
                     "Batch_4",     #15
                     "Batch_3",     #16
                     "Batch_3",     #17
                     "Batch_6",     #18
                     "Batch_5",     #19
                     "Batch_5",     #20
                     "Batch_1",     #21
                     "Batch_6"      #22
)


names(new.cluster.ids) <- levels(Sobj_integrated)
Sobj_integrated <- RenameIdents(Sobj_integrated, new.cluster.ids)
Sobj_integrated@meta.data$Batch <- Idents(Sobj_integrated)



## assigning cluster IDs on 0.1 resolution
Idents(Sobj_integrated) <- "integrated_snn_res.0.1"
new.cluster.ids <- levels(Sobj_integrated)
new.cluster.ids <- recode(new.cluster.ids, 
                          "0" =  "Basal",
                          "1" =  "Luminal_HR-pos",
                          "2" =  "Luminal_HR-neg",
                          "3" =  "Endothelial",
                          "4" =  "Fibroblast",
                          "5" =  "Endothelial_2",
                          "6" =  "Dying?",
                          "7" =  "Adipocyte",
                          "8" =  "Macrophage",
                          "9" =  "T-Cell",
                          "10" = "Vascular_Accessory",
                          "11" = "Neuronal_1",
                          "12" = "Neuronal_2",
                          "13" = "B-Cell"
)

names(new.cluster.ids) <- levels(Sobj_integrated)
Sobj_integrated <- RenameIdents(Sobj_integrated, new.cluster.ids)
Sobj_integrated@meta.data$Cluster <- Idents(Sobj_integrated)


## assigning broader Groups on 0.1 resolution
Idents(Sobj_integrated) <- "integrated_snn_res.0.1"
new.cluster.ids <- levels(Sobj_integrated)
new.cluster.ids <- recode(new.cluster.ids, 
                          "0" =  "Epithelial",
                          "1" =  "Epithelial",
                          "2" =  "Epithelial",
                          "3" =  "Vascular",
                          "4" =  "Stroma",
                          "5" =  "Vascular",
                          "6" =  "Epithelial",
                          "7" =  "Stroma",
                          "8" =  "Immune",
                          "9" =  "Immune",
                          "10" = "Vascular",
                          "11" = "Neuronal",
                          "12" = "Epithelial",
                          "13" = "Immune"
)

names(new.cluster.ids) <- levels(Sobj_integrated)
Sobj_integrated <- RenameIdents(Sobj_integrated, new.cluster.ids)
Sobj_integrated@meta.data$Group <- Idents(Sobj_integrated)







for (i in 1:length(Sobj_list)) {
  Sobj_integrated <- NormalizeData        (Sobj_integrated, verbose = T)
  Sobj_integrated <- FindVariableFeatures (Sobj_integrated, selection.method = "vst", nfeatures = 3000, verbose = T)
  Sobj_integrated <- ScaleData     (Sobj_integrated, 
                                    #vars.to.regress = c("nCount_RNA", "percent.mt"), 
                                    verbose = T) ### with regression
  Sobj_integrated <- RunPCA        (Sobj_integrated, npcs = 50, verbose = T)
  Sobj_integrated <- RunUMAP       (Sobj_integrated, reduction  = "pca", dims = 1:50)
  Sobj_integrated <- FindNeighbors (Sobj_integrated, reduction  = "pca", dims = 1:50)
  Sobj_integrated <- FindClusters  (Sobj_integrated, resolution = c(0.2))
}










library(Seurat)
library(tidyverse)

Meta <- readRDS("Output/MetaData_After_Trimming.rds")

MetaAdd <- select(Meta, Cluster = "Cluster", Group = "Group", CellType = "CellType", Detail = "CellType_Plus")

Sobj_integrated$integrated_snn_res.0.1 <- NULL
Sobj_integrated$integrated_snn_res.0.2 <- NULL
Sobj_integrated$integrated_snn_res.0.5 <- NULL
Sobj_integrated$integrated_snn_res.0.75 <- NULL
Sobj_integrated$integrated_snn_res.1 <- NULL
Sobj_integrated$integrated_snn_res.1.5 <- NULL
Sobj_integrated$integrated_snn_res.2 <- NULL
Sobj_integrated$seurat_clusters <- NULL
Sobj_integrated$LUMA1 <- NULL
Sobj_integrated$LUPR1 <- NULL
Sobj_integrated$MASC1 <- NULL
Sobj_integrated$STRM1 <- NULL
Sobj_integrated$S.Score <- NULL
Sobj_integrated$G2M.Score <- NULL
Sobj_integrated$Phase <- NULL
Sobj_integrated$TEST1 <- NULL
Sobj_integrated$Cluster <- NULL
Sobj_integrated$Group <- NULL

Sobj_integrated <- AddMetaData(Sobj_integrated, metadata = MetaAdd)

Sobj_integrated@meta.data$Group   <- replace_na(Sobj_integrated@meta.data$Group,   replace = "Discard")

Idents(Sobj_integrated) <- "Group"

Sobj_integrated <- subset(Sobj_integrated, idents = "Discard", invert = T)


Sobj_integrated <- subset(Sobj_integrated, idents = c("Lum_HR-pos", "Lum_HR-neg", "Basal", "Endothelial", "Fibroblast", "Adipocyte", "Myeloid", "Lymphoid"))




CellType = "Lum_HR-pos"


















