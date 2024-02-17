library(Seurat)
library(tidyverse)

### This workflow automates the processing of multiple 10x runs into one integrated dataset
setwd("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/")


### Preparation of our file lists that will contain our seurat objects
list.files("./CellRanger_Out/RNA/Final_Alignment/")
samples <- list.files("./CellRanger_Out/RNA/Final_Alignment/")
#samples <- samples[-20] #remove unused folder
#samples <- c("CF-318-813", "CF-2797", "CF-0404", "TM-9817", "TM-8249", "CF_3920", "TM-3937", "CF-428-112", "TM-7567", "TM_9469", "TM-2768", "CF-7780", "TM-6477", "CF_19301")
#samples <- c("CF_19301", "CF_3920", "CF-0404", "CF-249-347", "CF-2797", "CF-318-813", "CF-428-112", "CF-7780", "TM_1956", "TM_9469") # this is just a list of the sample names as found in the file structure

rawdata_list <- vector("list", length = length(samples))                 # this will contain our count matrices
Sobj_list    <- vector("list", length = length(samples))                 # our Seurat objects will be generated in here

### following loop reads in rawdata from cellranger output folders according to the sample name and then creates an Sobj in the list
for (i in seq_along(samples)) {
  rawdata_list[[i]] <- Read10X            (data.dir = file.path("CellRanger_Out/RNA/Final_Alignment", samples[i] , "filtered_feature_bc_matrix/"))
  Sobj_list   [[i]] <- CreateSeuratObject (rawdata_list[[i]], project = paste0("Nuc_", samples[i]), min.cells = 3, min.features = 200)
}

names(Sobj_list) <- samples # assign the correct names to the Sobj_List levels
rm(rawdata_list)




### Standard preprocessing and filtering 

# mito genes
for (i in seq_along(samples)) {
  Sobj_list[[i]][["percent.mt"]] <- PercentageFeatureSet(Sobj_list[[i]], pattern = "^MT-")
}

# VlnPlot of UMI and mito distribution
plot_list    <- vector("list", length = length(samples))

for (i in seq_along(samples)) {
  #plot_list[[i]] <- VlnPlot(Sobj_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot_list[[i]] <- VlnPlot(Sobj_list[[i]], features = c("nCount_RNA"), ncol = 3)
}

CombinePlots(plots = plot_list)

ggsave(file = "./Figures/12.05_Transgender_Normal_Nuclei_Alignment/UMI_Counts.png", plot = image, width = 16, height = 9, scale = 2)




for (i in seq_along(samples)) {
  print(VlnPlot     (Sobj_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  #print(VlnPlot     (Sobj_list[[i]], features = c("nFeature_RNA"), ncol = 3))
}

# filtering out high mito and high UMI cells
for (i in seq_along(samples)) {
  Sobj_list[[i]] <- subset(Sobj_list[[i]], subset = percent.mt < 2.5 & nCount_RNA < 30000)
}

saveRDS(Sobj_list, "Seurat_Objects/Top-Sample_Processing/Nuclei_Top-Samples_List.rds")

