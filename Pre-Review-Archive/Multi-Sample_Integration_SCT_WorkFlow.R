library(Seurat)
library(tidyverse)


Sobj_list <- "MY OBJECT LIST"

for (i in 1:length(Sobj_list)) {
  Sobj_list[[i]] <- SCTransform(Sobj_list[[i]], vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = T)
}




Sobj.features      <- SelectIntegrationFeatures (object.list = Sobj_list, nfeatures = 3000)
Sobj_list          <- PrepSCTIntegration        (object.list = Sobj_list, anchor.features = Sobj.features, verbose = T)
Sobj.SCTanchors    <- FindIntegrationAnchors    (object.list = Sobj_list,     normalization.method = "SCT", anchor.features = Sobj.features, verbose = T)
Sobj.SCTintegrated <- IntegrateData             (anchorset = Sobj.SCTanchors, normalization.method = "SCT", verbose = T)


Sobj.SCTintegrated <- RunPCA        (Sobj.SCTintegrated, verbose = T)
Sobj.SCTintegrated <- FindNeighbors (Sobj.SCTintegrated, reduction = "pca", dims = 1:50)
Sobj.SCTintegrated <- FindClusters  (Sobj.SCTintegrated, resolution = c(0.2, 0.5, 1))
Sobj.SCTintegrated <- RunUMAP       (Sobj.SCTintegrated, dims = 1:50)

DimPlot(Sobj.SCTintegrated, reduction = "umap")
