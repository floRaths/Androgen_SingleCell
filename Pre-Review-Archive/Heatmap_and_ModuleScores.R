

SCT <- do.call(rbind.data.frame, SCT_24_48_Resp[[1]])
RNA <- do.call(rbind.data.frame, RNA_24_48_Resp[[1]])




rownames(a) <- c()
im <- as.matrix(dplyr::select(a, "gene", "avg_diff", "Name") %>% column_to_rownames("gene"))
as.data.frame(im) %>% column_to_rownames("gene")
column_to_rownames(im, "gene")

pheatmap(im)


DoHeatmap(Downsample, features = a$gene, group.by = "Treat", assay = "SCT", slot = "scale.data",
          cells = WhichCells(Downsample, idents = 1))


DoHeatmap(Downsample, features = Garner_Markers$gene, group.by = "seurat_clusters", assay = "SCT", slot = "scale.data",
          cells = WhichCells(Downsample, idents = levels(Idents(Downsample)), downsample = 500))


FeaturePlot (Downsample, features = c("KRT23"), pt.size = 1, cols = c("lightblue", "mediumvioletred"), reduction = "umap", order = F, min.cutoff = "q9")









pancreas.integrated <- AddModuleScore  (pancreas.integrated, assay = "SCT", list((top_n(filter(All_Mrks_Cntrl, Cluster.Name == "MaSC.ECM"),  30, avg_logFC)$gene)), nbin = 25, name = "MaSC.ECM")
pancreas.integrated <- AddModuleScore  (pancreas.integrated, assay = "SCT", list((top_n(filter(All_Mrks_Cntrl, Cluster.Name == "Basal"),     30, avg_logFC)$gene)), nbin = 25, name = "BASAL")
pancreas.integrated <- AddModuleScore  (pancreas.integrated, assay = "SCT", list((top_n(filter(All_Mrks_Cntrl, Cluster.Name == "LUPR"),      50, avg_logFC)$gene)), nbin = 25, name = "LUPR")
pancreas.integrated <- AddModuleScore  (pancreas.integrated, assay = "SCT", list((top_n(filter(All_Mrks_Cntrl, Cluster.Name == "CC.High"),   30, avg_logFC)$gene)), nbin = 25, name = "CC.High")
pancreas.integrated <- AddModuleScore  (pancreas.integrated, assay = "SCT", list((top_n(filter(All_Mrks_Cntrl, Cluster.Name == "LUMA"),      50, avg_logFC)$gene)), nbin = 25, name = "LUMA")
pancreas.integrated <- AddModuleScore  (pancreas.integrated, assay = "SCT", list((top_n(filter(All_Mrks_Cntrl, Cluster.Name == "STRM.(Col)"),50, avg_logFC)$gene)), nbin = 25, name = "STRM.Col")
pancreas.integrated <- AddModuleScore  (pancreas.integrated, assay = "SCT", list((top_n(filter(All_Mrks_Cntrl, Cluster.Name == "STRM.(ECM)"),50, avg_logFC)$gene)), nbin = 25, name = "STRM.ECM")
pancreas.integrated <- AddModuleScore  (pancreas.integrated, assay = "SCT", list((top_n(filter(All_Mrks_Cntrl, Cluster.Name == "Low.UMI"),   50, avg_logFC)$gene)), nbin = 25, name = "Low.UMI")
pancreas.integrated <- AddModuleScore  (pancreas.integrated, assay = "SCT", list((top_n(filter(All_Mrks_Cntrl, Cluster.Name == "Krt.Env"),   50, avg_logFC)$gene)), nbin = 25, name = "Krt.Env")



Modules_List <- readRDS("~/Box Sync/Knott_Lab/Flo/Projects/Organoids/Analysis/2019.05.28_Estrogen_24-48h/utilities/Modules_List.rds")
pancreas.integrated <- AddModuleScore  (pancreas.integrated, list((top_n(Modules_List$LUMA, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "LUM")
pancreas.integrated <- AddModuleScore  (pancreas.integrated, list((top_n(Modules_List$LUPR, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "LUP")
pancreas.integrated <- AddModuleScore  (pancreas.integrated, list((top_n(Modules_List$MASC, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "MAS")
pancreas.integrated <- AddModuleScore  (pancreas.integrated, list((top_n(Modules_List$STRM, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "STR")




top_n(filter(RNA, Cluster == 4, p_val_adj < 0.05), 25, avg_logFC)
top_n(filter(SCT, Cluster == 4, p_val_adj < 0.05), 25, avg_diff)




X_NOR <- Combined_RNA[[3]]$`Luminal_HR-pos`

X_NOR <- LUMpos
X_NOR <- BasMyo
X_NOR <- LUMneg

EnhancedVolcano(X_NOR,
                lab = X_NOR$gene,
                x = 'avg_logFC',
                y = 'p_val_adj',
                xlim = c(-1.1, 1.1),
                FCcutoff = 0.25,
                transcriptPointSize = 1.5, pCutoff = 0.05, title = "Luminal-HR-neg", subtitle = "Response to Estrogen and Progesterone")


EnhancedVolcano(Y_NOR,
                lab = Y_NOR$gene,
                x = 'avg_logFC',
                y = 'p_val_adj',
                
                FCcutoff = 0.1,
                transcriptPointSize = 1.5, pCutoff = 0.05)



DimPlot(Hormone.integrated, reduction = "umap", group.by = "Sample", label = T)

FeaturePlot(Control_DataSets_Downsampled, features = "ESAM", pt.size = 1, order = T)

Hormone.integrated@meta.data$ <- TimeCourse_Final@meta.data$















