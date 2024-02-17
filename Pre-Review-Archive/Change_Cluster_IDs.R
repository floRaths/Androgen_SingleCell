Mobj <- subset(pancreas.query, idents = c("CF-3920_PreMrna", "Estrogen_Growth-TC_Pat-A"), downsample = 7000)

DefaultAssay(Mobj) <- "integrated"

#Mobj <- ScaleData  (Mobj, features = all.genes, verbose = T)
Mobj <- ScaleData  (Mobj, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = T)

Mobj <- RunPCA     (Mobj, npcs = 50, verbose = FALSE)
Mobj <- RunUMAP    (Mobj, reduction = "pca", dims = 1:50)

Mobj <- FindNeighbors (Mobj, reduction = "pca", dims = 1:50)
Mobj <- FindClusters  (Mobj, resolution = c(0.15, 0.5, 1))

DimPlot(Mobj, pt.size = 1, group.by = "Cluster", split.by = "Sample", label = T)
DimPlot(CF3920_Subset, pt.size = 1, group.by = "Cluster", split.by = "Hashtag", label = F)
DimPlot(s19301_Subset, pt.size = 1, group.by = "Cluster", split.by = "Hashtag", label = F)
DimPlot(Mobj, pt.size = 1, group.by = "integrated_snn_res.0.15", label = T)



# Secondary processing
Modules_List <- readRDS("~/Box Sync/Knott_Lab/Flo/Projects/Organoids/Analysis/2019.05.28_Estrogen_24-48h/utilities/Modules_List.rds")
Mobj <- AddModuleScore  (Mobj, assay = "RNA", list((top_n(Modules_List$LUMA, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "LUMA")
Mobj <- AddModuleScore  (Mobj, assay = "RNA", list((top_n(Modules_List$LUPR, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "LUPR")
Mobj <- AddModuleScore  (Mobj, assay = "RNA", list((top_n(Modules_List$MASC, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "MASC")
Mobj <- AddModuleScore  (Mobj, assay = "RNA", list((top_n(Modules_List$STRM, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "STRM")
rm(Modules_List)

Mobj <- AddModuleScore  (Mobj, assay = "RNA", list((top_n(filter(Final_MarkerSet, Name == "Basal"), 50, avg_logFC)$gene)), nbin = 25, name = "TEST")


FeaturePlot(Mobj, features = c("LUMA1", "LUPR1", "MASC1", "STRM1"), pt.size = 1, cols = c("lightblue", "mediumvioletred"), reduction = "umap", min.cutoff = "q6")
FeaturePlot(Mobj, features = c("ESR1"), pt.size = 1, cols = c("lightblue", "mediumvioletred"), reduction = "umap", min.cutoff = "q6", split.by = "orig.ident", order = T)

s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Mobj <- CellCycleScoring(Mobj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
rm(s.genes, g2m.genes)



Idents(Mobj) <- "integrated_snn_res.0.15"

new.cluster.ids <- levels(Mobj)

new.cluster.ids <- recode(new.cluster.ids, 
                          "0" =  "Fibroblast",
                          "1" =  "Basal",
                          "2" =  "Luminal_HR-pos",
                          "3" =  "Endothelial",
                          "4" =  "Luminal_HR-neg",
                          "5" =  "Basal_Myo",
                          "6" =  "Fibroblast",
                          "7" =  "Pre-Adipocyte",
                          "8" =  "Basal_Myo",
                          "9" =  "Fibroblast",
                          "10" = "Pot.Doublets",
                          "11" = "T-Cell",
                          "12" = "Endothelial2",
                          "13" = "Macrophage",
                          "14" = "Luminal_HR-pos",
                          "15" = "Basal",
                          "16" = "Adipocyte",
                          "17" = "B-Cell")


names(new.cluster.ids) <- levels(Mobj)
Mobj <- RenameIdents(Mobj, new.cluster.ids)
Mobj@meta.data$Cluster <- Idents(Mobj)





### assigning Patient_ID/replicate_ID to the Data
Idents(Mobj) <- "orig.ident"
new.cluster.ids <- c("CF-3920_Organoid", "CF-3920_Nuclei")
names(new.cluster.ids) <- levels(Mobj)
Mobj <- RenameIdents(Mobj, new.cluster.ids)
Mobj@meta.data$Sample <- Idents(Mobj)

### assigning Treatment groups to the Data
Idents(Mobj) <- "orig.ident"
new.cluster.ids <- c("Sample_CF-3920","Sample_CF-3920")
names(new.cluster.ids) <- levels(Mobj)
Mobj <- RenameIdents(Mobj, new.cluster.ids)
Mobj@meta.data$Patient <- Idents(Mobj)

### assigning Treatment groups to the Data
Idents(Mobj) <- "orig.ident"
new.cluster.ids <- c("Organoid","Nuclei")
names(new.cluster.ids) <- levels(Mobj)
Mobj <- RenameIdents(Mobj, new.cluster.ids)
Mobj@meta.data$SampleType <- Idents(Mobj)





A <- DimPlot (CF3920_Nuclei, pt.size = 1,  group.by = "Cluster", label = T)
B <- DimPlot(s19_301_Nuclei, pt.size = 1,  group.by = "Cluster", label = T)

image = CombinePlots(plots = list(A, B), ncol = 1)

DimPlot(Sobj_Pre, pt.size = 1,  group.by = "RNA_snn_res.0.25", label = T)


#save the plot in a variable image to be able to export to svg
image = DimPlot(Nuclei_Integrated, split.by = "Sample", ncol = 2, group.by = "Cluster", pt.size = 1)
#This actually save the plot in a image
ggsave(file = "./Figrues/10.24_Nuclei-Alignment_MAST_Responses/test.eps", plot = image, width=16, height=9)




Idents(Sobj) <- "RNA_snn_res.0.25"

new.cluster.ids <- levels(Sobj)

new.cluster.ids <- recode(new.cluster.ids, 
                          "0" = "Basal_Myo",
                          "1" = "Basal_Myo",
                          "2" = "Fibroblast",
                          "3" = "Basal",
                          "4" = "Luminal_HR-pos",
                          "5" = "Basal.CC",
                          "6" = "Luminal_HR-neg",
                          "7" = "Low.UMI",
                          "8" = "Endothelial")#,
                          #"9" = "Pot.Doublets",
                          "10" = "Endothelial2",
                          "11" = "Adipocyte",
                          "12" = "Pot.Doublet")

names(new.cluster.ids) <- levels(Sobj)
Sobj <- RenameIdents(Sobj, new.cluster.ids)
Sobj@meta.data$Cluster <- Idents(Sobj)







