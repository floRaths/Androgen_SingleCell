#########################B_Int###############################################
Sobj <- Hormone.Integrated
DefaultAssay(Sobj) <- "RNA"
Idents(Sobj) <- "Cluster"
Names <-       levels(Sobj)

#### assigning Names to treatments which allows detection of treatment responses per cluster #####
Sobj@meta.data$cluster.treat <- paste0(Sobj@active.ident, "_", Sobj@meta.data$Treat)
Idents(Sobj) <- "cluster.treat"


### Preparation of our file lists that will contain our seurat objects
Resp_Est <- vector("list", length = length(Names))    
Resp_Pro <- vector("list", length = length(Names))
Resp_E.P <- vector("list", length = length(Names))
Resp_E.S <- vector("list", length = length(Names))

### following loop reads in rawdata from cellranger output folders according to the sample name and then creates an Sobj in the list
for (i in 1:length(Names)) {
  Resp_Est[[i]] <- FindMarkers(Sobj, ident.1 = paste0(Names[i],"_Estrogen"), 
                                     ident.2 = paste0(Names[i],"_Control"),
                                assay = "RNA", slot = "data", test.use = "MAST", logfc.threshold = 0.05) %>% rownames_to_column("gene") %>% add_column(Cluster  = Names[i])
  
  Resp_Pro[[i]] <- FindMarkers(Sobj, ident.1 = paste0(Names[i],"_Progest"),  
                                     ident.2 = paste0(Names[i],"_Control"),
                                assay = "RNA", slot = "data", test.use = "MAST", logfc.threshold = 0.05) %>% rownames_to_column("gene") %>% add_column(Cluster =  Names[i])
  
  Resp_E.P[[i]] <- FindMarkers(Sobj, ident.1 = paste0(Names[i],"_Est+Pro"),  
                                     ident.2 = paste0(Names[i],"_Control"),
                                assay = "RNA", slot = "data", test.use = "MAST", logfc.threshold = 0.05) %>% rownames_to_column("gene") %>% add_column(Cluster =  Names[i])
  
  Resp_E.S[[i]] <- FindMarkers(Sobj, ident.1 = paste0(Names[i],"_Est_Shrt"),
                                     ident.2 = paste0(Names[i],"_Control"),
                                assay = "RNA", slot = "data", test.use = "MAST", logfc.threshold = 0.05) %>% rownames_to_column("gene") %>% add_column(Cluster =  Names[i])
}


names(Resp_Est) <- Names
names(Resp_Pro) <- Names
names(Resp_E.P) <- Names
names(Resp_E.S) <- Names
Treats <- c("Estrogen", "Progesterone", "Est+Pro", "Est_Short")

      Combined_MAST_Resp <- list(Resp_Est, Resp_Pro, Resp_E.P, Resp_E.S)
names(Combined_MAST_Resp) <- Treats




LUMneg <- FindMarkers(Sobj, ident.1 = paste0(Names[i],"_Est+Pro"),  
                             ident.2 = paste0(Names[i],"_Control"),
                             assay = "RNA", slot = "data", test.use = "MAST", logfc.threshold = 0.05) %>% rownames_to_column("gene") %>% add_column(Cluster =  Names[i])



rm(Resp_Est, Resp_Pro, Resp_E.P, Resp_E.S)