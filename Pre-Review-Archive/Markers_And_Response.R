### Finding Conserved Markers


pancreas.integrated <- Est24_48_Integrated

Clusters    <- c( 0,     1,    2,     3,         4,       5,        6,        7,       8, 9, 10)
Names <-       c("Myo","Stroma","LUMA","Stroma2","Low.UMI","Myo.CC","NEAT1","Stroma3","Endothelial", "LUPR", "T-Cells")
Conserved <- vector("list", length = 11)    
Idents(pancreas.integrated) <- "integrated_snn_res.0.2"

for (i in 1:11) {
  Conserved[[i]]  <- FindConservedMarkers(pancreas.integrated, ident.1 = Clusters[i],  
                                          grouping.var = "orig.ident", verbose = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i])
}

Conserved_24_48 <- Conserved
pancreas.integrated <- pan.int.full





Clusters    <- c( 0,     1,    2,     3,         4,       5,        6,        7)
Names <-       c("Basal",     #0
                 "Basal_Myo", #1
                 "Luminal_HR-neg",          #2
                 "Stroma",        #3
                 "Basal/Luminal",      #4
                 "Low.UMI",      #5
                 "Basal.CC",      #6
                 "Luminal_HR-pos")         #7

Names <-       c("Basal/Luminal", "Low.UMI", "Basal.CC", "Luminal_HR-pos")         #7


Markers_SCT <- vector("list", length = length(Clusters))    
Markers_RNA <- vector("list", length = length(Clusters))    
Idents(A1) <- "integrated_snn_res.0.15"

for (i in 1:length(Clusters)) {
  Markers_SCT[[i]]  <- FindMarkers(A1, ident.1 = Clusters[i], verbose = T,
                                   assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
  Markers_RNA[[i]]  <- FindMarkers(A1, ident.1 = Clusters[i], verbose = T,
                                   assay = "RNA")                      %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
}




Conserved_SCT <- vector("list", length = length(Clusters))    
Conserved_RNA <- vector("list", length = length(Clusters))    
Idents(Sobj) <- "integrated_snn_res.0.15"


for (i in 1:length(Clusters)) {
  Conserved_SCT[[i]]  <- FindConservedMarkers(Sobj, ident.1 = Clusters[i], grouping.var = "orig.ident", verbose = T,
                                              assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
  
  Conserved_RNA[[i]]  <- FindConservedMarkers(Sobj, ident.1 = Clusters[i], grouping.var = "orig.ident", verbose = T,
                                              assay = "RNA") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
}


for (i in 1:length(Clusters)) {
  Conserved_SCT[[i]]  <- FindConservedMarkers(Sobj, ident.1 = Clusters[i], grouping.var = "orig.ident", verbose = T,
                                              assay = "SCT", slot = "scale.data") #%>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
  
  Conserved_RNA[[i]]  <- FindConservedMarkers(Sobj, ident.1 = Clusters[i], grouping.var = "orig.ident", verbose = T,
                                              assay = "RNA") #%>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
}

#########################Integr###############################################

pancreas.integrated <- Full_DataSet_SCT_Integr
Idents(pancreas.integrated) <- "SCT_snn_res.0.2"
DefaultAssay(pancreas.integrated) <- "RNA"

#### assigning clusters to treatments which allows detection of treatment responses per cluster #####
pancreas.integrated@meta.data$cluster.treat <- paste0(pancreas.integrated@active.ident, "_", pancreas.integrated@meta.data$Hashtag)
Idents(pancreas.integrated) <- "cluster.treat"


### Preparation of our file lists that will contain our seurat objects
Clusters    <- c( 0,     1,    2,     3,         4,       5,        6,        7,       8,             9,      10, 11)
Names <-       c("Myo","Stroma","LUMA","Stroma2","Low.UMI","Myo.CC","NEAT1","Stroma3","Endothelial", "LUPR", "T-Cells")

Resp_Est24 <- vector("list", length = 12)                 # our Seurat objects will be generated in here
Resp_Est48 <- vector("list", length = 12)


### following loop reads in rawdata from cellranger output folders according to the sample name and then creates an Sobj in the list
for (i in 1:11) {
  Resp_Est24[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Est24"), ident.2 = paste0(Clusters[i],"_Cntrl"),
                               print.bar = T, assay = "RNA") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i])
  
  Resp_Est48[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Est48"),  ident.2 = paste0(Clusters[i],"_Cntrl"),
                               print.bar = T, assay = "RNA") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i])
  
  }

names(Resp_Est24) <- Names
names(Resp_Est48) <- Names

RNA_24_48_Resp <- list(Resp_Est24, Resp_Est48)


#########################Integr###############################################

pancreas.integrated <- Est24_A
Idents(pancreas.integrated) <- "integrated_snn_res.0.2"
DefaultAssay(pancreas.integrated) <- "SCT"

#### assigning clusters to treatments which allows detection of treatment responses per cluster #####
pancreas.integrated@meta.data$cluster.treat <- paste0(pancreas.integrated@active.ident, "_", pancreas.integrated@meta.data$Hashtag)
Idents(pancreas.integrated) <- "cluster.treat"


### Preparation of our file lists that will contain our seurat objects
Clusters    <- c( 0,     1,    2,     3,         4,       5,        6,        7,       8,             9,      10)
Names <-       c("Myo","Stroma","LUMA","Stroma2","Low.UMI","Myo.CC","NEAT1","Stroma3","Endothelial", "LUPR", "T-Cells")

Resp_Est24 <- vector("list", length = 11)                 # our Seurat objects will be generated in here
Resp_Est48 <- vector("list", length = 11)


### following loop reads in rawdata from cellranger output folders according to the sample name and then creates an Sobj in the list
for (i in 1:9) {
  Resp_Est24[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Est24"), ident.2 = paste0(Clusters[i],"_Cntrl"),
                                 print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
  
  Resp_Est48[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Est48"),  ident.2 = paste0(Clusters[i],"_Cntrl"),
                                 print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
  
}

names(Resp_Est24) <- Names
names(Resp_Est48) <- Names

Est_24_48_A_Resp <- list(Resp_Est24, Resp_Est48)

#########################Integr###############################################

pancreas.integrated <- Est24_B
Idents(pancreas.integrated) <- "integrated_snn_res.0.2"
DefaultAssay(pancreas.integrated) <- "SCT"

#### assigning clusters to treatments which allows detection of treatment responses per cluster #####
pancreas.integrated@meta.data$cluster.treat <- paste0(pancreas.integrated@active.ident, "_", pancreas.integrated@meta.data$Hashtag)
Idents(pancreas.integrated) <- "cluster.treat"


### Preparation of our file lists that will contain our seurat objects
Clusters    <- c( 0,     1,    2,     3,         4,       5,        6,        7,       8,             9,      10)
Names <-       c("Myo","Stroma","LUMA","Stroma2","Low.UMI","Myo.CC","NEAT1","Stroma3","Endothelial", "LUPR", "T-Cells")

Resp_Est24 <- vector("list", length = 11)                 # our Seurat objects will be generated in here
Resp_Est48 <- vector("list", length = 11)


### following loop reads in rawdata from cellranger output folders according to the sample name and then creates an Sobj in the list
for (i in 1:9) {
  Resp_Est24[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Est24"), ident.2 = paste0(Clusters[i],"_Cntrl"),
                                 print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
  
  Resp_Est48[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Est48"),  ident.2 = paste0(Clusters[i],"_Cntrl"),
                                 print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
  
}

names(Resp_Est24) <- Names
names(Resp_Est48) <- Names

Est_24_48_B_Resp <- list(Resp_Est24, Resp_Est48)






















#########################Integr###############################################

pancreas.integrated <- pan.int.full
Idents(pancreas.integrated) <- "integrated_snn_res.0.2"
DefaultAssay(pancreas.integrated) <- "SCT"

#### assigning clusters to treatments which allows detection of treatment responses per cluster #####
pancreas.integrated@meta.data$cluster.treat <- paste0(pancreas.integrated@active.ident, "_", pancreas.integrated@meta.data$Treat)
Idents(pancreas.integrated) <- "cluster.treat"


### Preparation of our file lists that will contain our seurat objects
Clusters    <- c(0, 1, 2, 3, 4, 5, 6, 7, 8)
Names <-       c("MaSC","Myo","LUPR","Myo.MaSC","Stroma","Myo.LUM","Low.UMI","Myo.CC","LUMA")# this is just a list of the cluster names

Resp_Est <- vector("list", length = 9)                 # our Seurat objects will be generated in here
Resp_Pro <- vector("list", length = 9)
Resp_E.P <- vector("list", length = 9)
Resp_E.S <- vector("list", length = 9)

### following loop reads in rawdata from cellranger output folders according to the sample name and then creates an Sobj in the list
for (i in 1:9) {
  Resp_Est[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Estrogen"), ident.2 = paste0(Clusters[i],"_Control"),
                               print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
  
  Resp_Pro[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Progest"),  ident.2 = paste0(Clusters[i],"_Control"),
                               print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
  
  Resp_E.P[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Est+Pro"),  ident.2 = paste0(Clusters[i],"_Control"),
                               print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
  
  Resp_E.S[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Est_Shrt"), ident.2 = paste0(Clusters[i],"_Control"),
                               print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
}


names(Resp_Est) <- Names
names(Resp_Pro) <- Names
names(Resp_E.P) <- Names
names(Resp_E.S) <- Names

Integr_Resp <- list(Resp_Est, Resp_Pro, Resp_E.P, Resp_E.S)






#########################A_Int###############################################

pancreas.integrated <- Pat_A
Idents(pancreas.integrated) <- "integrated_snn_res.0.2"
DefaultAssay(pancreas.integrated) <- "SCT"

#### assigning clusters to treatments which allows detection of treatment responses per cluster #####
pancreas.integrated@meta.data$cluster.treat <- paste0(pancreas.integrated@active.ident, "_", pancreas.integrated@meta.data$Treat)
Idents(pancreas.integrated) <- "cluster.treat"


### Preparation of our file lists that will contain our seurat objects
Clusters    <- c(0, 1, 2, 3, 4, 5, 6, 7, 8)
Names <-       c("MaSC","Myo","LUPR","Myo.MaSC","Stroma","Myo.LUM","Low.UMI","Myo.CC","LUMA")# this is just a list of the cluster names

Resp_Est <- vector("list", length = 9)                 # our Seurat objects will be generated in here
Resp_Pro <- vector("list", length = 9)
Resp_E.P <- vector("list", length = 9)
Resp_E.S <- vector("list", length = 9)

### following loop reads in rawdata from cellranger output folders according to the sample name and then creates an Sobj in the list
for (i in 1:9) {
  Resp_Est[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Estrogen"), ident.2 = paste0(Clusters[i],"_Control"),
                               print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
  
  Resp_Pro[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Progest"),  ident.2 = paste0(Clusters[i],"_Control"),
                               print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
  
  Resp_E.P[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Est+Pro"),  ident.2 = paste0(Clusters[i],"_Control"),
                               print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
  
  Resp_E.S[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Est_Shrt"), ident.2 = paste0(Clusters[i],"_Control"),
                               print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
}


names(Resp_Est) <- Names
names(Resp_Pro) <- Names
names(Resp_E.P) <- Names
names(Resp_E.S) <- Names

Pat_A_Resp <- list(Resp_Est, Resp_Pro, Resp_E.P, Resp_E.S)






#########################B_Int###############################################

pancreas.integrated <- Pat_B
Idents(pancreas.integrated) <- "integrated_snn_res.0.2"
DefaultAssay(pancreas.integrated) <- "SCT"

#### assigning clusters to treatments which allows detection of treatment responses per cluster #####
pancreas.integrated@meta.data$cluster.treat <- paste0(pancreas.integrated@active.ident, "_", pancreas.integrated@meta.data$Treat)
Idents(pancreas.integrated) <- "cluster.treat"


### Preparation of our file lists that will contain our seurat objects
Clusters    <- c(0, 1, 2, 3, 4, 5, 6, 7, 8)
Names <-       c("MaSC","Myo","LUPR","Myo.MaSC","Stroma","Myo.LUM","Low.UMI","Myo.CC","LUMA")# this is just a list of the cluster names
Resp_Est <- vector("list", length = 9)                 # our Seurat objects will be generated in here
Resp_Pro <- vector("list", length = 9)
Resp_E.P <- vector("list", length = 9)
Resp_E.S <- vector("list", length = 9)

### following loop reads in rawdata from cellranger output folders according to the sample name and then creates an Sobj in the list
for (i in 1:9) {
  Resp_Est[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Estrogen"), ident.2 = paste0(Clusters[i],"_Control"),
                               print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
  
  Resp_Pro[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Progest"),  ident.2 = paste0(Clusters[i],"_Control"),
                               print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
  
  Resp_E.P[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Est+Pro"),  ident.2 = paste0(Clusters[i],"_Control"),
                               print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
  
  Resp_E.S[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Clusters[i],"_Est_Shrt"), ident.2 = paste0(Clusters[i],"_Control"),
                               print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Clusters[i], Name = Names[i])
}


names(Resp_Est) <- Names
names(Resp_Pro) <- Names
names(Resp_E.P) <- Names
names(Resp_E.S) <- Names

Pat_B_Resp <- list(Resp_Est, Resp_Pro, Resp_E.P, Resp_E.S)






for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- SCTransform(pancreas.list[[i]], vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
}


pancreas.features <-   SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
pancreas.list <-       PrepSCTIntegration(object.list = pancreas.list, anchor.features = pancreas.features, verbose = T)
pancreas.anchors <-    FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT", anchor.features = pancreas.features, verbose = T)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", verbose = T)


pancreas.integrated <- RunPCA(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:50)
pancreas.integrated <- FindNeighbors (pancreas.integrated, reduction = "pca", dims = 1:50)
pancreas.integrated <- FindClusters  (pancreas.integrated, resolution = c(0.2))

DimPlot(pancreas.integrated, reduction = "umap", group.by = "integrated_snn_res.0.2", 
        label = T, split.by = "orig.ident", ncol = 5, pt.size = 1)


















































for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- SCTransform(pancreas.list[[i]], vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
}


pancreas.features <-   SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
pancreas.list <-       PrepSCTIntegration(object.list = pancreas.list, anchor.features = pancreas.features, verbose = T)
pancreas.anchors <-    FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT", anchor.features = pancreas.features, verbose = T)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", verbose = T)


pancreas.integrated <- RunPCA(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:50)
pancreas.integrated <- FindNeighbors (pancreas.integrated, reduction = "pca", dims = 1:50)
pancreas.integrated <- FindClusters  (pancreas.integrated, resolution = c(0.2))

DimPlot(pancreas.integrated, reduction = "umap", group.by = "integrated_snn_res.0.2", 
        label = T, split.by = "orig.ident", ncol = 5, pt.size = 1)



















