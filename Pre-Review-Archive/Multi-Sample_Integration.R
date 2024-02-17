library(Seurat)
library(tidyverse)

### This workflow automates the processing of multiple 10x runs into one integrated dataset
setwd("~/Box Sync/Knott_Lab/Flo/Projects/Organoids/Analysis/2019.02.17_Oganoid_Hormone_TimeCourse/")


### Preparation of our file lists that will contain our seurat objects
samples <- c("A1", "B1") # this is just a list of the sample names
rawdata_list <- vector("list", length = length(samples))                 # this will contain our count matrices
Sobj_list    <- vector("list", length = length(samples))                 # our Seurat objects will be generated in here

### following loop reads in rawdata from cellranger output folders according to the sample name and then creates an Sobj in the list
for (i in seq_along(samples)) {
  rawdata_list[[i]] <- Read10X            (data.dir = file.path("CellRanger_Out/V3/", samples[i] , "/filtered_feature_bc_matrix/"))
  Sobj_list[[i]]    <- CreateSeuratObject (rawdata_list[[i]], project = paste0("Organoid_Hormone_", samples[i]), min.cells = 1, min.features = 1)
}

names(Sobj_list) <- samples # assign the correct names to the Sobj_List levels
rm(rawdata_list)


### Standart preprocessing and filtering 

# mito genes
for (i in seq_along(samples)) {
  Sobj_list[[i]][["percent.mt"]] <- PercentageFeatureSet(Sobj_list[[i]], pattern = "^MT-")
}

# VlnPlot of UMI and mito distribution
for (i in seq_along(samples)) {
  print(VlnPlot     (Sobj_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
}

# filtering out high mito and high UMI cells
for (i in seq_along(samples)) {
  Sobj_list[[i]] <- subset(Sobj_list[[i]], subset = percent.mt < 7.5 & nCount_RNA < 50000)
}


