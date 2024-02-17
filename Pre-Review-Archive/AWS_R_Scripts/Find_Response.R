library(Seurat)
library(tidyverse)

Sobj               <- readRDS("~/mnt/Sobj_Scaled_Proc.rds")


######### Response by CellType
{

Idents(Sobj)       <- "CellType"
DefaultAssay(Sobj) <- "RNA"

Marker_list    <- vector("list", length = length(levels(Sobj)))

for (i in seq_along(Marker_list)) {
  Marker_list[[i]]    <- FindMarkers(subset(Sobj, idents = levels(Sobj)[i]),     
            group.by = "hash.ID",  
            ident.1  = "AR-high", 
            ident.2  = "Control", 
            assay    = "RNA", slot = "data",
            logfc.threshold = 0.25,
            test.use = "MAST", ) %>% 
    rownames_to_column("Gene") %>% 
    add_column(Cluster  = levels(Sobj)[i])
}


Bind <- do.call(rbind.data.frame, Marker_list)

}


######### Subcluster Treatment Response
{
  
  Idents(Sobj)       <- "CellType"
  DefaultAssay(Sobj) <- "RNA"
  
  Marker_list    <- vector("list", length = length(levels(Sobj)))
  
  for (i in seq_along(levels(Sobj))) {
    
    subset <- subset(Sobj, idents = levels(Sobj)[i])
    Idents(subset)       <- "Subcluster"
    
    SubMarker_list    <- vector("list", length = length(levels(subset)))
    
    for (j in seq_along(Marker_list)) {
      SubMarker_list[[j]]    <- FindMarkers(subset(subset, idents = levels(subset)[j]),
                                         group.by = "Type",
                                         ident.1  = "TM",
                                         ident.2  = "CF", 
                                         assay    = "RNA", 
                                         slot     = "data",
                                         logfc.threshold = 0, 
                                         test.use = "MAST") %>% 
        rownames_to_column("Gene") %>% 
        add_column(Subcluster  = levels(subset)[j])#, CellType = levels(Sobj)[i])
    }
    Marker_list[[i]] <- do.call(rbind.data.frame, SubMarker_list)
  }
  Bind <- do.call(rbind.data.frame, Marker_list)
}



######### Subcluster Markers
{

Idents(Sobj)       <- "CellType"
DefaultAssay(Sobj) <- "RNA"

Marker_list    <- vector("list", length = length(levels(Sobj)))

for (i in seq_along(levels(Sobj))) {
  
  subset <- subset(Sobj, idents = levels(Sobj)[i])
  Idents(subset)       <- "Subcluster"
  
  Marker_list[[i]]    <- FindAllMarkers(subset,     
                                     assay = "RNA", 
                                     slot  = "data",
                                     logfc.threshold = 0.25,
                                     test.use = "MAST", 
                                     only.pos = T) %>% 
    rownames_to_column("Gene") %>% 
    add_column(CellType  = levels(Sobj)[i])
}

Bind <- do.call(rbind.data.frame, Marker_list)

}
