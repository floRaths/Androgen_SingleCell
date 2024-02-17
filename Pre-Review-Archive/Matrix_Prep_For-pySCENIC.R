
category <- "CellType_ALL"

Idents(Sobj_integrated) <- category

levels <- levels(Sobj_integrated)
S      <- SplitObject(Sobj_integrated, split.by = category)


for (i in 1:length(levels)) {
  DefaultAssay(S[[levels[i]]]) <- "RNA"
  S[[levels[i]]]               <- FindVariableFeatures(S[[levels[i]]], assay = "RNA", nfeatures = 4000, selection.method = "vst")
  S[[levels[i]]]               <- NormalizeData       (S[[levels[i]]], assay = "RNA", verbose = T)
  saveRDS(subset(S[[levels[i]]], features = VariableFeatures(S[[levels[i]]]))@assays$RNA@data, paste0("~/Documents/SCENIC/pySCENIC_Full/scenicdata/4k_Var_Genes/", levels[i], "_mat_4k_norm_varft.rds"))
  }


files <- list.files("/home/ubuntu/mnt/pySCENIC_Full/scenicdata/4k_Var_Genes/")

for (i in 1:length(files)) {
  mat <-               readRDS(paste0("/home/ubuntu/mnt/pySCENIC_Full/scenicdata/4k_Var_Genes/Groups/", files[i]))
  write.csv(t(as.matrix(mat)), paste0("/home/ubuntu/mnt/pySCENIC_Full/scenicdata/4k_Var_Genes/Groups/", files[i] ,".csv"))
}