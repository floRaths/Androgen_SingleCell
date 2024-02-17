


filter(dplyr::select(as.data.frame(VariableFeatures(Sobj)), Gene = 1), Gene %in% variableTFs$Gene) %>% arrange(Gene)


Sobj <- NormalizeData(Sobj_integrated, assay = "RNA", verbose = T)

Sobj <- FindVariableFeatures(Sobj, assay = "RNA", selection.method = "vst", nfeatures = 12000)

features <- c(VariableFeatures(Sobj), variableTFs$Gene) %>% unique()

saveRDS(subset(Sobj, features = VariableFeatures(Sobj))@assays$RNA@data, "~/Documents/SCENIC/pySCENIC_Full/scenicdata/TFs_VarFeatures/TF_12k_Varft.rds")


#mat, "~/Documents/SCENIC/pySCENIC_Full/scenicdata/mat_8k_norm_varft.rds")


exmat <- readRDS("~/int/exprMat_1k.rds")

write.csv(rownames_to_column(as.data.frame(t(Matrix::as.matrix(exmat))), "Gene"), "int/exmat_1k_t.csv")
y <- rownames_to_column(as.data.frame(  Matrix::as.matrix(exmat)),  "Gene")




write.csv(t(Matrix::as.matrix(exmat)), "/home/ubuntu/mnt/pySCENIC_Full/scenicdata/mat_TF_4k_cleaned.csv")




celltype = levels(Sobj_integrated)

for (i in 1:length(celltype)) {
  
  Subset <- subset(Sobj_integrated, idents = celltype[i])
  
  genes <- LUMexpr %>% filter(CellType == celltype[i], perc >= 5) %>% pull(TF) %>% unique() 
  
  Subset@assays$RNA@data[genes, ] %>% as.matrix() %>% t() %>% write.csv(paste0("/home/rathsf/scNuclei-Libraries/Analysis/Reboot/Output/CellType_SCENIC/Matrices/Matrix_", celltype[i], ".csv"))
  
}




hgnc_tfs <- read_csv("~/Documents/SCENIC/pySCENIC_Full/scenicdata/Resources/hs_hgnc_tfs.txt", col_names = "TF")



