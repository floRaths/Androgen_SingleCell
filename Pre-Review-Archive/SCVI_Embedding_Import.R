Sobj_integrated <- readRDS("Seurat_Objects/Global_DataSets/Nuclei_All_Samples_Batch_Reduced.rds")
Sobj_integrated <- readRDS("Seurat_Objects/DataSubsets/CellType_Batch_Integration/Added_DimReducts/Neuronal_DimReds.rds")
saveRDS(Sobj_integrated,   paste0("Seurat_Objects/DataSubsets/CellType_Batch_Integration/Added_DimReducts/", path2, "_DimReds.rds"))


Sobj_integrated <- AddMetaData(Sobj_integrated, metadata = Sobj_Full@meta.data)


Sobj_integrated <- subset(Sobj_Full, idents = "Vasc_Acc")

path1 <- "CellType"
path2 <- "Vasc_Acc"
path3 <- "_umap"

SCVI_umap <- read_csv(paste0("~/Documents/SCVI/misc_single_cell/", path1, "/SCVI_", path2 , path3, ".tsv"), col_names = "SCVI") %>% 
  separate(SCVI, into = c("SCVI_1", "SCVI_2"), sep = "\t", convert = T)

SCVI_meta <- read_csv(paste0("~/Documents/SCVI/misc_single_cell/", path1, "/adata.obs_", path2, ".csv"))

SCVI_umap <- SCVI_umap %>% bind_cols(dplyr::select(SCVI_meta, index)) %>% filter(index %in% colnames(Sobj_integrated)) %>% column_to_rownames("index") %>% as.matrix()

Sobj_integrated[["SCVI"]] <- CreateDimReducObject(embeddings = as.matrix(SCVI_umap), key = "SCVI", assay = "RNA")





group <- "Sample"
label <- T

#DimPlot(Sobj_integrated, reduction = "umap",   group.by = group, split.by = "Type", ncol = 1, label = F, repel = T)   + theme(legend.position = "none") + ggtitle("Seurat Integration")
#DimPlot(Sobj_integrated, reduction = "SCVI",   group.by = group, split.by = "Type", ncol = 1, label = F, repel = T)   + theme(legend.position = "none") + ggtitle("SCVI")
#DimPlot(Sobj_integrated, reduction = "SCENIC", group.by = group, split.by = "Type", ncol = 1, label = F, repel = T)   + ggtitle("Scenic AUC_Scores")

image <- 
DimPlot(Sobj_integrated, reduction = "umap",   group.by = group, split.by = "Type", ncol = 2, label = T, repel = T, pt.size = 1)  #     + theme(legend.position = "none") + ggtitle("Seurat Integration")
image <- 
DimPlot(Sobj_integrated, reduction = "SCVI",   group.by = group, split.by = "Type", ncol = 2, label = label, repel = T, pt.size = 1)#   + theme(legend.position = "none") + ggtitle("SCVI")
DimPlot(Sobj_integrated, reduction = "SCENIC", group.by = group, split.by = "Type", ncol = 1, label = label, repel = T)   + ggtitle("Scenic AUC_Scores")

quick <- 
  p1 + p2 + p3 + plot_layout(ncol = 3)

save2(paste0("Batch_Correction_", path2, "_", group), 1)

celltype <- path2
save(name = paste0("UMAP_Seurat_", celltype), 1)




levels(Sobj_Full)
Idents(Sobj_Full) <- "CellType_All"

lymph <- WhichCells(Sobj_Full, idents = c("B_Cell", "CD8_T", "CD4_T", "NK", "HSC"))
myelo <- WhichCells(Sobj_Full, idents = c("Monocyte", "Macrophage", "moDC", "DC"))








group <- c("Sample", "Group", "Batch", "Prep", "Type", "CellType", "CellType_Plus")

for (i in 1:length(group)) {
  
  dir.create   (paste0("Figures/4.29_SCVI_and_Workflow/PNG/",path2), recursive = T, showWarnings = F)
  figure.path = paste0("Figures/4.29_SCVI_and_Workflow/PNG/",path2)
  
  p1 <- DimPlot(Sobj_integrated, reduction = "umap",   group.by = group[i], split.by = "Type", ncol = 1, label = F, repel = T)   + theme(legend.position = "none") + ggtitle("Seurat Integration")
  p2 <- DimPlot(Sobj_integrated, reduction = "SCVI",   group.by = group[i], split.by = "Type", ncol = 1, label = F, repel = T)   + theme(legend.position = "none") + ggtitle("SCVI")
  p3 <- DimPlot(Sobj_integrated, reduction = "SCENIC", group.by = group[i], split.by = "Type", ncol = 1, label = F, repel = T)   + ggtitle("Scenic AUC_Scores")
  
  quick <- 
    p1 + p2 + p3 + plot_layout(ncol = 3)
  
  save2(paste0("Batch_Correction_", path2, "_", group[i]), 1)
  
  }








