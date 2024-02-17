Sobj_list_Unsampled <- Sobj_list
Sobj_list[[1]] <- subset(Sobj_list[[1]], ident = "Organoid_Hormone_CF_3920_premRNA", downsample = 5583)

### Initiating integration procedure

for (i in 1:length(Sobj_list)) {
  Sobj_list[[i]] <- NormalizeData        (Sobj_list[[i]], verbose = T)
  Sobj_list[[i]] <- FindVariableFeatures (Sobj_list[[i]], selection.method = "vst", nfeatures = 5000, verbose = T)
}


Sobj_anchors    <- FindIntegrationAnchors (object.list = Sobj_list,    dims = 1:50)
Sobj_integrated <- IntegrateData          (anchorset   = Sobj_anchors, dims = 1:50)
rm(Sobj_anchors, Sobj_list)

DefaultAssay (Sobj_integrated) <- "integrated"
Sobj_integrated <- ScaleData  (Sobj_integrated, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = T) ### with regression
#saveRDS(Sobj_integrated, file = "Outputs/DataSets/Sobj_integrated_from_V3-Scaled_and_regressed.rds")

### continuing with regressed dataset
Sobj_integrated <- RunPCA        (Sobj_integrated, npcs = 50, verbose = FALSE)
Sobj_integrated <- RunUMAP       (Sobj_integrated, reduction  = "pca", dims = 1:50)
Sobj_integrated <- FindNeighbors (Sobj_integrated, reduction  = "pca", dims = 1:50)
Sobj_integrated <- FindClusters  (Sobj_integrated, resolution = c(0.2, 0.5, 1))






### assigning Patient_ID/replicate_ID to the Data
Idents(Sobj_integrated) <- "orig.ident"
new.cluster.ids <- c("CF_3920", "TM_9469", "TM_1956", "CF_19301")
names(new.cluster.ids) <- levels(Sobj_integrated)
Sobj_integrated <- RenameIdents(Sobj_integrated, new.cluster.ids)
Sobj_integrated@meta.data$Sample <- Idents(Sobj_integrated)

### assigning Treatment groups to the Data
Idents(Sobj_integrated) <- "orig.ident"
new.cluster.ids <- c("Manual","Manual","Sorted","Sorted")
names(new.cluster.ids) <- levels(Sobj_integrated)
Sobj_integrated <- RenameIdents(Sobj_integrated, new.cluster.ids)
Sobj_integrated@meta.data$Counting <- Idents(Sobj_integrated)

Idents(Sobj_integrated) <- "integrated_snn_res.0.2"
new.cluster.ids <- levels(Sobj_integrated)
new.cluster.ids <- recode(new.cluster.ids, 
                          "0" =  "Luminal_HR-pos",
                          "1" =  "Basal_Myo",
                          "2" =  "Luminal_HR-neg",
                          "3" =  "Endothelial",
                          "4" =  "Fibroblast",
                          "5" =  "Adipocyte",
                          "6" =  "Macrophage",
                          "7" =  "T-Cell",
                          "8" =  "Neuronal",
                          "9" =  "Luminal_HR-pos")

names(new.cluster.ids) <- levels(Sobj_integrated)
Sobj_integrated <- RenameIdents(Sobj_integrated, new.cluster.ids)
Sobj_integrated@meta.data$Cluster <- Idents(Sobj_integrated)


# Secondary processing
Modules_List <- readRDS("~/Box/Knott_Lab/Flo/Projects/Organoids/Analysis/2019.02.17_Oganoid_Hormone_TimeCourse//utilities/Modules_List.rds")
Sobj_integrated <- AddModuleScore  (Sobj_integrated, assay = "RNA", list((top_n(Modules_List$LUMA, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "LUMA")
Sobj_integrated <- AddModuleScore  (Sobj_integrated, assay = "RNA", list((top_n(Modules_List$LUPR, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "LUPR")
Sobj_integrated <- AddModuleScore  (Sobj_integrated, assay = "RNA", list((top_n(Modules_List$MASC, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "MASC")
Sobj_integrated <- AddModuleScore  (Sobj_integrated, assay = "RNA", list((top_n(Modules_List$STRM, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "STRM")
rm(Modules_List)

s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Sobj_integrated <- CellCycleScoring(Sobj_integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = F, assay = "RNA")
rm(s.genes, g2m.genes)


Sobj_integrated <- AddModuleScore  (Sobj_integrated, assay = "RNA", list((top_n(filter(Final_MarkerSet, Name == "Macrophage"), 50, avg_logFC)$gene)), nbin = 25, name = "TEST")


DimPlot(Sobj_integrated, pt.size = 1, cells = WhichCells(Sobj_integrated, idents = 4), split.by = NULL, group.by = "integrated_snn_res.0.2", label = T, ncol = 2)

DimPlot(Sobj_integrated, pt.size = 1, group.by = "Cluster", label = T, split.by = "Sample", ncol = 2)
DimPlot(Sobj_integrated, pt.size = 1, group.by = "Sample", label = T)

select.cells <- CellSelector(plot = plot)
Idents(Sobj_integrated, cells = select.cells) <- "pot.doublets"


VlnPlot(Sobj_integrated, idents = c("pot.doublets", 0, 2), features = "nCount_RNA")


FeaturePlot(Sobj_integrated, features = c("LUMA1", "LUPR1", "MASC1", "STRM1"), pt.size = 1, cols = c("lightblue", "mediumvioletred"), reduction = "umap", min.cutoff = "q6", order = T)
FeaturePlot(Sobj_integrated, features = c("ESAM", "LUMA1"), pt.size = 1, cols = c("lightblue", "mediumvioletred"), reduction = "umap", min.cutoff = "q6", order = T)


Idents(Sobj_integrated) <- "Cluster"
Unknown <- FindMarkers(Sobj_integrated, ident.1 =  "Unknown", test.use = "MAST", assay = "RNA", slot = "data")

saveRDS(Sobj_integrated, "./Outputs/DataSets/Sobj_integratedd_Processed.rds")
Sobj_integrated <- readRDS("./Outputs/DataSets/Sobj_integratedd_Processed.rds")































