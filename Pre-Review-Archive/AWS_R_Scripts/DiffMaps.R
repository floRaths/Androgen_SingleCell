library(Seurat)
library(tidyverse)
library(destiny)
#library(future)


#options(future.globals.maxSize= 3670016000)
#plan("multiprocess", workers = 32)

Sobj_integrated <- readRDS("~/mnt/Sobj_Integrated_Scaled_Proc.rds")
DefaultAssay(Sobj_integrated) <- "RNA"

#Sobj_integrated <- DietSeurat(Sobj_integrated, counts = T, data = T, scale.data = F, assays = "RNA")

Idents(Sobj_integrated) <- "Cluster"
Idents <- levels(Sobj_integrated)

Matrices <- vector("list", length = length(Idents))

for (i in 1:length(Idents)) {
  Matrices[[i]] <- as.matrix(subset(Sobj_integrated, idents = Idents[i])@assays$RNA@data) %>% t
}

names(Matrices) <- Idents

DiffMaps <- vector("list", length = length(Idents))

for (i in 1:length(Idents)) {
  DiffMaps[[i]] <- DiffusionMap(Matrices[[i]])
}


saveRDS(DiffMaps, "~/mnt/DiffMap.rds")
