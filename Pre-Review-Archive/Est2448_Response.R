

pancreas.integrated <- SCT_Integr
Idents(pancreas.integrated) <- "Cluster"
Names <-      levels(SCT_Integr)

#### assigning clusters to treatments which allows detection of treatment responses per cluster #####
pancreas.integrated@meta.data$cluster.treat <- paste0(pancreas.integrated@active.ident, "_", pancreas.integrated@meta.data$Hashtag)
Idents(pancreas.integrated) <- "cluster.treat"

### Preparation of our file lists that will contain our seurat objects
Resp_Est24_RNA <- vector("list", length = length(levels(SCT_Integr)))                 
Resp_Est48_RNA <- vector("list", length = length(levels(SCT_Integr)))
Resp_Est24_SCT <- vector("list", length = length(levels(SCT_Integr)))                 
Resp_Est48_SCT <- vector("list", length = length(levels(SCT_Integr)))

### following loop reads in rawdata from cellranger output folders according to the sample name and then creates an Sobj in the list
for (i in 1:14) {
  #Resp_Est24_SCT[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Names[i],"_Est24"), ident.2 = paste0(Names[i],"_Cntrl"),
                                 #print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Names[i])
  
  #Resp_Est48_SCT[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Names[i],"_Est48"), ident.2 = paste0(Names[i],"_Cntrl"),
                                 #print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Names[i])
  
  Resp_Est24_RNA[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Names[i],"_Est24"), ident.2 = paste0(Names[i],"_Cntrl"),
                                 print.bar = T, assay = "RNA") %>% rownames_to_column("gene") %>% add_column(Cluster = Names[i])
  
  Resp_Est48_RNA[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Names[i],"_Est48"), ident.2 = paste0(Names[i],"_Cntrl"),
                                 print.bar = T, assay = "RNA") %>% rownames_to_column("gene") %>% add_column(Cluster = Names[i])
}


names(Resp_Est24_RNA) <- Names
names(Resp_Est48_RNA) <- Names
names(Resp_Est24_SCT) <- Names
names(Resp_Est48_SCT) <- Names

Normalized_Est24_48_RNA_Resp <- list(Resp_Est24_RNA, Resp_Est48_RNA)
Combined_Est24_48_SCT_Resp <- list(Resp_Est24_SCT, Resp_Est48_SCT)


#################################
#################################
#################################


pancreas.integrated <- Test
Idents(pancreas.integrated) <- "Cluster"
Names <-      levels(SCT_Integr)

#### assigning clusters to treatments which allows detection of treatment responses per cluster #####
pancreas.integrated@meta.data$cluster.treat <- paste0(pancreas.integrated@active.ident, "_", pancreas.integrated@meta.data$Treat)
Idents(pancreas.integrated) <- "cluster.treat"

### Preparation of our file lists that will contain our seurat objects
Resp_Est24_RNA <- vector("list", length = length(levels(SCT_Integr)))                 
Resp_Est48_RNA <- vector("list", length = length(levels(SCT_Integr)))
Resp_Est24_SCT <- vector("list", length = length(levels(SCT_Integr)))                 
Resp_Est48_SCT <- vector("list", length = length(levels(SCT_Integr)))

### following loop reads in rawdata from cellranger output folders according to the sample name and then creates an Sobj in the list
for (i in 1:14) {
  Resp_Est24_SCT[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Names[i],"_Est24"), ident.2 = paste0(Names[i],"_Cntrl"),
                                     print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Names[i])
  
  Resp_Est48_SCT[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Names[i],"_Est48"), ident.2 = paste0(Names[i],"_Cntrl"),
                                     print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Names[i])
  
  Resp_Est24_RNA[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Names[i],"_Est24"), ident.2 = paste0(Names[i],"_Cntrl"),
                                     print.bar = T, assay = "RNA") %>% rownames_to_column("gene") %>% add_column(Cluster = Names[i])
  
  Resp_Est48_RNA[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Names[i],"_Est48"), ident.2 = paste0(Names[i],"_Cntrl"),
                                     print.bar = T, assay = "RNA") %>% rownames_to_column("gene") %>% add_column(Cluster = Names[i])
}


names(Resp_Est24_RNA) <- Names
names(Resp_Est48_RNA) <- Names
names(Resp_Est24_SCT) <- Names
names(Resp_Est48_SCT) <- Names

Pat_A_Est24_48_RNA_Resp <- list(Resp_Est24_RNA, Resp_Est48_RNA)
Pat_A_Est24_48_SCT_Resp <- list(Resp_Est24_SCT, Resp_Est48_SCT)


#################################
#################################
#################################


pancreas.integrated <- Pat_B
Idents(pancreas.integrated) <- "Cluster"
Names <-      levels(SCT_Integr)

#### assigning clusters to treatments which allows detection of treatment responses per cluster #####
pancreas.integrated@meta.data$cluster.treat <- paste0(pancreas.integrated@active.ident, "_", pancreas.integrated@meta.data$Hashtag)
Idents(pancreas.integrated) <- "cluster.treat"

### Preparation of our file lists that will contain our seurat objects
Resp_Est24_RNA <- vector("list", length = length(levels(SCT_Integr)))                 
Resp_Est48_RNA <- vector("list", length = length(levels(SCT_Integr)))
Resp_Est24_SCT <- vector("list", length = length(levels(SCT_Integr)))                 
Resp_Est48_SCT <- vector("list", length = length(levels(SCT_Integr)))

### following loop reads in rawdata from cellranger output folders according to the sample name and then creates an Sobj in the list
for (i in 1:14) {
  Resp_Est24_SCT[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Names[i],"_Est24"), ident.2 = paste0(Names[i],"_Cntrl"),
                                     print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Names[i])
  
  Resp_Est48_SCT[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Names[i],"_Est48"), ident.2 = paste0(Names[i],"_Cntrl"),
                                     print.bar = T, assay = "SCT", slot = "scale.data") %>% rownames_to_column("gene") %>% add_column(Cluster = Names[i])
  
  Resp_Est24_RNA[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Names[i],"_Est24"), ident.2 = paste0(Names[i],"_Cntrl"),
                                     print.bar = T, assay = "RNA") %>% rownames_to_column("gene") %>% add_column(Cluster = Names[i])
  
  Resp_Est48_RNA[[i]] <- FindMarkers(pancreas.integrated, ident.1 = paste0(Names[i],"_Est48"), ident.2 = paste0(Names[i],"_Cntrl"),
                                     print.bar = T, assay = "RNA") %>% rownames_to_column("gene") %>% add_column(Cluster = Names[i])
}


names(Resp_Est24_RNA) <- Names
names(Resp_Est48_RNA) <- Names
names(Resp_Est24_SCT) <- Names
names(Resp_Est48_SCT) <- Names

Pat_B_Est24_48_RNA_Resp <- list(Resp_Est24_RNA, Resp_Est48_RNA)
Pat_B_Est24_48_SCT_Resp <- list(Resp_Est24_SCT, Resp_Est48_SCT)


