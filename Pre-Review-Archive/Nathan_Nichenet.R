
### upregulated AR_responsive genes
filter(rownames_to_column(AR_Marks, "Gene"), p_val_adj <= 0.05, avg_logFC >= 0.5)


receptor_modules <- c("~/Documents/NicheNet/receptor_genes_long/")

list.files(receptor_modules)
receptors <- substr(list.files(receptor_modules), 1,nchar(list.files(receptor_modules))-4)
receptor_list <- vector("list", length = length(receptors))  

#read_tsv(paste0("utilities/receptor_genes/AR.txt")

for (i in 1:length(list.files(receptor_modules))) {
  receptor_list[[i]] <- read_tsv(paste0(receptor_modules, 
                                        receptors[i], ".txt"), 
                                 col_names = "Gene") %>% mutate(Receptor = receptors[i])
}

names(receptor_list) <- receptors


Recp_Bind <- do.call(rbind.data.frame, receptor_list) %>% dplyr::select(Receptor, Gene)



filter(Recp_Bind, Gene == "SRRM2")





### list of receptor genes without .txt
receptors <- substr(list.files("utilities/receptor_genes/"), 1,nchar(list.files("utilities/receptor_genes/"))-4)





read_tsv("utilities/receptor_genes/PIP.txt", col_names = "Gene")





sender <- 
  filter(rownames_to_column(AR_Marks, "Gene"), p_val_adj <= 0.05, avg_logFC >= 0.5) %>% 
    left_join(Human_Secretome, by = "Gene") %>% arrange(Category) %>% head(10)


filter(weighted_networks$lr_sig, from == "KLK3") %>% arrange(-weight)

filter(weighted_networks$lr_sig, to == "RYR2") %>% arrange(-weight)



filter(TopTable, 
       Gene %in% filter(weighted_networks$lr_sig, from %in% sender$Gene, weight >= 0.25)$to) %>%
  arrange(Gene)





detach("package:sceasy", unload = TRUE)

library(sceasy)
library(reticulate)



path_to_python <- "/Users/rathsf/miniconda3/bin/python3.7"
reticulate::use_python(path_to_python)
reticulate::py_config()

reticulate::use_condaenv("/Users/rathsf/miniconda3/envs/sceasy2")

use_condaenv('sceasy2')
loompy <- reticulate::import('loompy')

sceasy::convertFormat(Sobj, from="seurat", to="anndata", outFile='~/Documents/SCVI/Sobj_Anndata.h5ad', main_layer = "counts")

sceasy::convertFormat(obj = Sobj_integrated, 




                      
                      
                      
                      
                      
                      
                      
                      
Nicheput <- read_csv("~/Documents/NicheNet/output_short.csv")
head(Nicheput)


colnames <- colnames(column_to_rownames(Nicheput, "barcodes"))
Meta     <- Sobj_integrated@meta.data %>% rownames_to_column("barcodes") %>% dplyr::select(barcodes, Type, CellType)


mat <- 
Meta %>% left_join(Nicheput, by = "barcodes") %>% dplyr::select(-barcodes) %>%
  group_by(CellType, Type) %>% summarise_at(vars(colnames), mean) %>% unite("NewCol", Type, CellType) %>% 
  column_to_rownames("NewCol")

score <- 
mat %>% rownames_to_column("Type") %>% 
  pivot_longer(cols = colnames) %>% 
  group_by(name) %>% summarise(mean = mean(value)) %>% 
  top_n(200, mean) %>% pull(name)

#cells <- WhichCells(Sobj_integrated, idents = "LUM_HR-pos", downsample = 1000)
#mat <- column_to_rownames(filter(Nicheput, barcodes %in% cells), "barcodes")

show <- colnames(mat) %>% intersect(receptors_oi$to)

image <- 
pheatmap((mat[,show]), 
         scale = "column", 
         color = viridis::viridis(n = 100), 
         cutree_cols = 1, 
         cluster_rows = F, 
         annotation_col = Ann,
         fontsize_row = 10, 
         fontsize_col = 8, 
         treeheight_col = 10, 
         gaps_row = c(2,4,6,8,10,12,14,16,18), 
         treeheight_row = 2,
         )


Ann <- 
  weighted_networks$lr_sig %>% 
  filter(to %in% show, from %in% secr) %>%
  pivot_wider(names_from = "from", values_from =  "weight") %>%
  column_to_rownames("to") 

Ann[!is.na(Ann)] <- "TRUE"

Ann

pheatmap(t(Ann), scale = "column", cluster_cols = F)



as.data.frame(colnames((mat[,receptors_oi$Receptor]))) %>% 
  dplyr::select(to = 1) %>% left_join(Ann, by =  "to") %>% arrange(to)








Ann <- 
weighted_networks$lr_sig %>% 
  filter(from %in% inters) %>%
  group_by(from) %>% 
  top_n(10, weight) %>% mutate(weight = to) %>% 
  pivot_wider(names_from = "from", values_from =  "to") %>%
  column_to_rownames("weight")


Ann[!is.na(Ann)] <- "YES"
Ann[is.na(Ann)] <- "NO"



save2("NicheNet_NewSet_SecrAnn", scale = 1.5)


save2 <- function(name, scale) {
  
  ggsave(plot   = image,
         path     = file.path("Figures/4.15_NicheNet_SCENIC_Refined/PNG/"),
         filename = paste0(name,".png"),
         width    = 16, 
         height   = 9,
         scale    = scale,
         dpi      = 300)
}




