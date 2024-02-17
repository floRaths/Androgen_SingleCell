
### proof of error catch (b and x dont exist...so error should stop loop)
for (i in 1:10) {
  
  skip_to_next <- FALSE
  
  # Note that print(b) fails since b doesn't exist
  
  tryCatch(print(b), error = function(e) { skip_to_next <<- TRUE})
  tryCatch(print(x), error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }     
}





for (i in 1:length(Subset_list)) {
  
  skip_to_next <- FALSE
  
  Subset_list[[i]] <- SplitObject(Subset_list[[i]], split.by = "Batch")
  
  
      for (s in 1:length(Subset_list[[i]])) {
        DefaultAssay                                  (Subset_list[[i]][[s]]) <- "RNA"
        Subset_list[[i]][[s]] <- NormalizeData        (Subset_list[[i]][[s]], verbose = T)
        Subset_list[[i]][[s]] <- FindVariableFeatures (Subset_list[[i]][[s]], selection.method = "vst", nfeatures = 3000, verbose = T)
        }

    
  tryCatch(Sobj_anchors     <- FindIntegrationAnchors (object.list = Subset_list[[i]],    dims = 1:50), 
           error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }
  
  tryCatch(Subset_list[[i]] <- IntegrateData          (anchorset   = Sobj_anchors, dims = 1:50), 
           error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }
  
  rm(Sobj_anchors)
  
  DefaultAssay                      (Subset_list[[i]]) <- "integrated"
  Subset_list[[i]] <- ScaleData     (Subset_list[[i]], vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = T, features = )
  Subset_list[[i]] <- RunPCA        (Subset_list[[i]], npcs = 50, verbose = FALSE)
  Subset_list[[i]] <- RunUMAP       (Subset_list[[i]], reduction  = "pca", dims = 1:50)
  Subset_list[[i]] <- FindNeighbors (Subset_list[[i]], reduction  = "pca", dims = 1:50)
  Subset_list[[i]] <- FindClusters  (Subset_list[[i]], resolution = c(0.2, 0.5, 1))

}





