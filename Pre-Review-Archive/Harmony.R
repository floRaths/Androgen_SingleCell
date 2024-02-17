### Run Harmony on Seurat Object
library(harmony)
library(Seurat)



DefaultAssay (Sobj_integrated) <- "RNA"

Sobj_integrated <- NormalizeData(Sobj_integrated, verbose = T)
Sobj_integrated <- FindVariableFeatures (Sobj_integrated, selection.method = "vst", nfeatures = 2000, verbose = T)

Sobj_integrated <- ScaleData  (Sobj_integrated, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = T) ### with regression

### continuing with regressed dataset
Sobj_integrated <- RunPCA        (Sobj_integrated, npcs = 50, verbose = T)
Sobj_integrated <- RunHarmony    (Sobj_integrated, "Batch")
Sobj_integrated <- RunUMAP       (Sobj_integrated, reduction  = "harmony", dims = 1:50)
#Sobj_integrated <- FindNeighbors (Sobj_integrated, reduction  = "harmony", dims = 1:50)
#Sobj_integrated <- FindClusters  (Sobj_integrated, resolution = c(0.1, 0.2, 0.5))


#plot1 <- 
DimPlot(Sobj_integrated, 
        #cells = WhichCells(Sobj_integrated, idents = c("TM-3937", "CF-2797")),
        reduction = "umap", 
        group.by  = "Type",
        #group.by  = "CellType",
        #split.by = "Type",
        pt.size = 1, 
        label = F,
        repel = T, 
        label.size = 4,
        #split.by  = "Sample", ncol = 6, 
        
        ) 


plot1 + plot2 + patchwork::plot_layout(widths = c(1,2))


plot2 <- 
DimPlot(Sobj_integrated, 
        #cells = WhichCells(Sobj_integrated, idents = c("TM-3937", "CF-2797")),
        reduction = "umap", 
        #group.by  = "RNA_snn_res.0.2",
        group.by  = "Sample",
        split.by = "Type",
        pt.size = 0.5, 
        label = T,
        repel = T, 
        label.size = 4,
        #split.by  = "Sample", ncol = 6, 
        #cols = viridis::magma(n = 7)
) 

plot3 <- 
DimPlot(Sobj_integrated, 
        #cells = WhichCells(Sobj_integrated, idents = c("TM-3937", "CF-2797")),
        reduction = "umap", 
        #group.by  = "RNA_snn_res.0.2",
        group.by  = "Inclusion",
        split.by = "Type",
        pt.size = 0.5, 
        label = T,
        repel = T, 
        label.size = 10,
        #split.by  = "Sample", ncol = 6, 
        #cols = viridis::magma(n = 7)
) 

plot4 <- 
DimPlot(Sobj_integrated, 
        #cells = WhichCells(Sobj_integrated, idents = c("TM-3937", "CF-2797")),
        reduction = "umap", 
        #group.by  = "RNA_snn_res.0.2",
        group.by  = "Contr.Samples",
        split.by = "Type",
        pt.size = 0.5, 
        label = T,
        repel = T, 
        label.size = 10,
        #split.by  = "Sample", ncol = 6, 
        cols = viridis::magma(n = 8)
) 



plot1 + plot2 + plot3 + plot4 + plot_layout(ncol = 2)





tfs <- table(Secretory_System_Components$Subsystem) %>% names()



for (i in seq_along(tfs)) {
        
        skip_to_next <- FALSE
        
        module <- percentile %>% filter(TF == tfs[i], percentile_rank > 99) %>% pull(target)
        
        tryCatch(Sobj_integrated <- AddModuleScore  (Sobj_integrated, 
                                    assay = "RNA", 
                                    #list(filter(Corr_Bind, str_detect(Clust.Name, "^KLK") | str_detect(Clust.Name, "^AR"))$Gene), 
                                    list(module), 
                                    nbin = 25, name = paste0(tfs[i],"_")), 
                 
                 error = function(e) { skip_to_next <<- TRUE})
        
        if(skip_to_next) { next }
        }




