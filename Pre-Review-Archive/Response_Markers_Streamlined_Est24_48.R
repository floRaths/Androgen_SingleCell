#########################B_Int###############################################
Sobj <- Full_DataSet_SCT_Annotated
DefaultAssay(Sobj) <- "RNA"
Idents(Sobj) <- "Cluster"
Names <-       levels(Sobj)

#### assigning Names to treatments which allows detection of treatment responses per cluster #####
Sobj@meta.data$cluster.treat <- paste0(Sobj@active.ident, "_", Sobj@meta.data$Hashtag)
Idents(Sobj) <- "cluster.treat"


### Preparation of our file lists that will contain our seurat objects
Resp_Est24 <- vector("list", length = length(Names))    
Resp_Est48 <- vector("list", length = length(Names))


### following loop reads in rawdata from cellranger output folders according to the sample name and then creates an Sobj in the list
for (i in 1:length(Names)) {
  Resp_Est24[[i]] <- FindMarkers(Sobj, ident.1 = paste0(Names[i],"_Est24"), 
                                       ident.2 = paste0(Names[i],"_Cntrl"),
                                 assay = "RNA", slot = "data", test.use = "MAST", logfc.threshold = 0.01) %>% 
    rownames_to_column("gene") %>% 
    add_column(Cluster  = Names[i])
  
  Resp_Est48[[i]] <- FindMarkers(Sobj, ident.1 = paste0(Names[i],"_Est48"),  
                                       ident.2 = paste0(Names[i],"_Cntrl"),
                                 assay = "RNA", slot = "data", test.use = "MAST", logfc.threshold = 0.01) %>% 
    rownames_to_column("gene") %>% 
    add_column(Cluster =  Names[i])
  
}


names(Resp_Est24) <- Names
names(Resp_Est48) <- Names

Treats <- c("Est_24", "Est_48")

      Combined_MAST_Resp <- list(Resp_Est24, Resp_Est48)
names(Combined_MAST_Resp) <- Treats




# Patient A ---------------------------------------------------------------

Idents(Full_DataSet_SCT_Annotated) <- "orig.ident"

Pat_A <- subset(Full_DataSet_SCT_Annotated, idents = "Estrogen_24-48h_Pat-A")
Pat_B <- subset(Full_DataSet_SCT_Annotated, idents = "Estrogen_24-48h_Pat-B")


Sobj <- Pat_A
DefaultAssay(Sobj) <- "RNA"
Idents(Sobj) <- "Cluster"
Names <-       levels(Sobj)

#### assigning Names to treatments which allows detection of treatment responses per cluster #####
Sobj@meta.data$cluster.treat <- paste0(Sobj@active.ident, "_", Sobj@meta.data$Hashtag)
Idents(Sobj) <- "cluster.treat"


### Preparation of our file lists that will contain our seurat objects
Resp_Est24 <- vector("list", length = length(Names))    
Resp_Est48 <- vector("list", length = length(Names))


### following loop reads in rawdata from cellranger output folders according to the sample name and then creates an Sobj in the list
for (i in 1:length(Names)) {
  Resp_Est24[[i]] <- FindMarkers(Sobj, ident.1 = paste0(Names[i],"_Est24"), 
                                 ident.2 = paste0(Names[i],"_Cntrl"),
                                 assay = "RNA", slot = "data", test.use = "MAST", logfc.threshold = 0.01) %>% 
    rownames_to_column("gene") %>% 
    add_column(Cluster  = Names[i])
  
  Resp_Est48[[i]] <- FindMarkers(Sobj, ident.1 = paste0(Names[i],"_Est48"),  
                                 ident.2 = paste0(Names[i],"_Cntrl"),
                                 assay = "RNA", slot = "data", test.use = "MAST", logfc.threshold = 0.01) %>% 
    rownames_to_column("gene") %>% 
    add_column(Cluster =  Names[i])
  
}


names(Resp_Est24) <- Names
names(Resp_Est48) <- Names

Treats <- c("Est_24", "Est_48")

Pat_A_MAST_Resp <- list(Resp_Est24, Resp_Est48)
names(Pat_A_MAST_Resp) <- Treats






# patient B ---------------------------------------------------------------

Sobj <- Pat_B
DefaultAssay(Sobj) <- "RNA"
Idents(Sobj) <- "Cluster"
Names <-       levels(Sobj)

#### assigning Names to treatments which allows detection of treatment responses per cluster #####
Sobj@meta.data$cluster.treat <- paste0(Sobj@active.ident, "_", Sobj@meta.data$Hashtag)
Idents(Sobj) <- "cluster.treat"


### Preparation of our file lists that will contain our seurat objects
Resp_Est24 <- vector("list", length = length(Names))    
Resp_Est48 <- vector("list", length = length(Names))


### following loop reads in rawdata from cellranger output folders according to the sample name and then creates an Sobj in the list
for (i in 1:length(Names)) {
  Resp_Est24[[i]] <- FindMarkers(Sobj, ident.1 = paste0(Names[i],"_Est24"), 
                                 ident.2 = paste0(Names[i],"_Cntrl"),
                                 assay = "RNA", slot = "data", test.use = "MAST", logfc.threshold = 0.01) %>% 
    rownames_to_column("gene") %>% 
    add_column(Cluster  = Names[i])
  
  Resp_Est48[[i]] <- FindMarkers(Sobj, ident.1 = paste0(Names[i],"_Est48"),  
                                 ident.2 = paste0(Names[i],"_Cntrl"),
                                 assay = "RNA", slot = "data", test.use = "MAST", logfc.threshold = 0.01) %>% 
    rownames_to_column("gene") %>% 
    add_column(Cluster =  Names[i])
  
}


names(Resp_Est24) <- Names
names(Resp_Est48) <- Names

Treats <- c("Est_24", "Est_48")

Pat_B_MAST_Resp <- list(Resp_Est24, Resp_Est48)
names(Pat_B_MAST_Resp) <- Treats
