PercFiltered <- readRDS("2021_ReRun/Percent_GeneExpression_CellType_Type.rds")

Idents(Sobj) <- "CellType"

for (i in 1:length(levels(Sobj))) {
  
  Subset <- subset(Sobj, idents = levels(Sobj)[i])
  genes <- PercFiltered %>% filter(CellType == levels(Sobj)[i], CF > 5 | TM > 5) %>% pull(Gene) %>% unique() 
  Subset@assays$RNA@data[genes, ] %>% as.matrix() %>% t() %>% write.csv(paste0("Matrices_Grnboost/Matrix_", levels(Sobj)[i], ".csv"))
  
}






