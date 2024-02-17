#celltype <- "Adipocyte"

library(Seurat)
library(tidyverse)

S               <- readRDS("~/mnt/Nuclei_All_Samples_Batch_Reduced_May06.rds")
#Marker_Bind                   <- readRDS("~/mnt/Marker_Bind_CT_All_significant_only_May06.rds")


category <- "CellType_All"

Idents(S) <- category


S = SplitObject(Sobj_integrated, split.by = "CellType")

#celltype <- levels(S)
celltype <- levels(Sobj_integrated)

#for (i in 1:1) {
for (i in 1:length(celltype)) {
  #genes <- TFs %>% filter(CellType == celltype[i]) %>% pull(Gene) %>% unique()
  DefaultAssay(S[[celltype[i]]]) <- "RNA"
  genes <- List %>% filter(CellType == celltype[i], expression == "expressed", perc >= 5) %>% pull(TF) %>% unique()
  #S[[celltype[i]]]  <- NormalizeData       (S[[celltype[i]]], assay = "RNA", verbose = T)
  #S[[celltype[i]]]  <- FindVariableFeatures(S[[celltype[i]]], assay = "RNA", nfeatures = 10000, selection.method = "vst")
  
  #features <- c(S[[celltype[i]]]@assays$RNA@var.features, genes) %>% unique() %>% sort()
  
  S[[celltype[i]]]  <- ScaleData(S[[celltype[i]]], 
                                 assay = "RNA", 
                                 features = genes,
                                 vars.to.regress = c("nCount_RNA", "percent.mt"), 
                                 verbose = T)
  
  #saveRDS            (S[[celltype[i]]],                        paste0("~/mnt/data/CT_All_Indiv.Scaled_", celltype[i], ".rds"))
  write.csv(as.matrix(S[[celltype[i]]]@assays$RNA@scale.data), paste0("scNuclei-Libraries/Analysis/Reboot/Output/Indiv.Scaled_5perc_", celltype[i], "_scaled_data.csv", sep = ''))
  write.csv          (S[[celltype[i]]]@meta.data,              paste0("scNuclei-Libraries/Analysis/Reboot/Output/Indiv.Scaled_5perc_", celltype[i], "_metadata.csv",    sep = ''))
  
  }



Sobj_plot <- readRDS("~/Documents/SCENIC/pySCENIC_Full/scenicdata/CT_All_4kvar/Seurat_Subsets/CT_All_Subset_LUM_HR-pos.rds")

Sobj_plot <- ScaleData(Sobj_plot, 
                       features = features,
                       assay = "RNA",
                       vars.to.regress = c("nCount_RNA", "percent.mt"),
                       verbose = T)


features <- 
  Marker_Bind %>% filter(CellType == "LUM_HR-pos") %>% pull(Gene) %>% c(Sobj_plot@assays$RNA@var.features) %>% unique()





#Average_TF_CF <- 
levels <- levels(Sobj_integrated)
Expression_list    <- vector("list", length = length(levels(Sobj_integrated)))

for (i in 1:length(levels)) {

cells_CF <- WhichCells(Sobj_integrated, idents = levels[i])

Expression_list[[i]] <- 
  Sobj_integrated@assays$RNA@counts[rownames(Sobj_integrated), cells_CF] %>% 
  as.data.frame() %>% 
  rownames_to_column("TF") %>% 
  as_tibble() %>% 
  pivot_longer(cols = all_of(cells_CF)) %>% 
  left_join(dplyr::select(rownames_to_column(Sobj_integrated@meta.data, "name"), name, Type, CellType_All), by = "name") %>%
  mutate(expression = ifelse(value >= 1, "expressed", "absent")) %>% 
  group_by(Type, TF, expression) %>%
  summarise(n = n()) %>% 
  group_by(Type, TF) %>% 
  mutate(total = sum(n)) %>% mutate(perc = (n/total)*100) %>% add_column("CellType" = levels[i]) 

#saveRDS(x, paste0("~/mnt/Gene_Expression_", levels[i], ".rds"))

}

saveRDS(x, paste0("~/mnt/Gene_Expression_", levels[i], ".rds"))







Average_TF_Expression <- rbind(Average_TF_CF, Average_TF_TM) %>% ungroup()


Average_TF_Expression %>% filter(expression == "expressed", perc >= 50) %>% group_by(Type, TF) %>% summarise(n = n()) %>% pivot_wider(names_from = Type, values_from = n) %>% arrange(-CF)








celltype <- "LUM_HR-pos"

adjacencies <- read_tsv(paste0("~/Documents/SCENIC/pySCENIC_Full/scenicdata/CT_All_4kvar/Results/", celltype, "/expr_mat.adjacencies.tsv"))

percentile <- 
  adjacencies %>% 
  group_by(TF) %>% 
  mutate(percentile_rank = ntile(importance, 100))


TF_candidates <- 
Average_TF_Expression %>% ungroup() %>% 
  filter(CellType_All == celltype, expression == "expressed", perc >= 30) %>% 
  arrange(Type, -perc) %>% pull(TF) %>% unique() %>% sort()



Marker_Bind %>% filter(CellType == celltype) %>% 
  pull(Gene) %>% 
  unique() %>% 
  intersect(TF_candidates)




x <- Marker_Bind %>% filter(CellType == celltype) %>% dplyr::select(target = 1,3, 6)


query <- "TRPS1"

#p4 <- 
percentile %>% ungroup %>%
  left_join(x, by = "target") %>% 
  filter(TF == query, !is.na(avg_logFC), percentile_rank > 95) %>%
  ggplot(aes(x = avg_logFC, y = importance, label = target)) + 
  geom_point() + 
  geom_text(nudge_y = 1, size = 3, check_overlap = T) + ggtitle(paste(query, " in ", celltype))

#p1 + p2 + p3 + p4





percentile %>% ungroup %>%
  left_join(LUMp_DE, by = "target") %>% group_by()








files <- list.files("Output/Gene_Expression_Ratios/")

file_list <- vector("list", length = length(files))  

for (i in 1:length(files)) {
  file_list[[i]] <- readRDS(paste0("Output/Gene_Expression_Ratios/", files[i])) %>% dplyr::rename("Gene" = "TF")# 
  
}

names(file_list) <- files %>% substr(start = 17, stop = 50)
 

file_list$Adipocyte.rds %>% filter(expression ==  "expressed", perc >= 2.5) %>% pull(Gene) %>% unique() %>% length()



GeneEx_Bind <- do.call(rbind.data.frame, file_list)


all_genes <- GeneEx_Bind %>% pull(Gene) %>% unique()

x <- GeneEx_Bind %>% ungroup() %>% filter(expression ==  "expressed") %>% group_by(Type, CellType) %>% filter(perc >= 2.5) #%>% pull(Gene) %>% unique() %>% length()


out <- genes %>% outersect(all_genes)




outersect <- function(x, y) {
  sort(c(x[!x%in%y],
         y[!y%in%x]))
}





subset(Sobj_Full, idents = "LUM_HR-pos", features = genes)




c(VariableFeatures(Sobj_Sub), TFs$Gene)




Sobj_integrated@assays$RNA@var.features




