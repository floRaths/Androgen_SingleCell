rm(list = ls())
dev.off(dev.list()["RStudioGD"])


library(tidyverse)
setwd("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/")


make_narrow <- function() {
  A <- 
    Motifs_Narrow_TM_JASPAR %>% 
    select(idxATAC:last_col()) %>% 
    filter(GeneRna %in% perc_Type_filtered5$Gene) %>%
    pivot_longer(names_to = "Motif", values_to = "value", starts_with("Motif")) %>% 
    filter(!is.na(value)) %>% 
    separate(Motif, into = c("A", "B"), sep = "\\.DistanceTo") %>%
    separate(A, into = c("Motif", "C"), sep = "\\.Sig") %>% 
    unite(X, C, B) %>% 
    mutate(dtbs = "JASPAR") %>%
    pivot_wider(names_from = X, values_from = value) %>% 
    dplyr::rename(DistanceToSummit = "NA_Summit", Signal = "nal_NA") %>% 
    left_join(Motif_Join, by = c("Motif", "dtbs")) %>% 
    group_by(TF, idxATAC) %>% 
    mutate(logFC = log2(Signal/signalBigWig), Type = "Trans")
  
  B <- 
    Motifs_Narrow_CF_JASPAR %>% 
    select(idxATAC:last_col()) %>% 
    filter(GeneRna %in% perc_Type_filtered5$Gene) %>%
    pivot_longer(names_to = "Motif", values_to = "value", starts_with("Motif")) %>% 
    filter(!is.na(value)) %>% 
    separate(Motif, into = c("A", "B"), sep = "\\.DistanceTo") %>%
    separate(A, into = c("Motif", "C"), sep = "\\.Sig") %>% 
    unite(X, C, B) %>% 
    mutate(dtbs = "JASPAR") %>%
    pivot_wider(names_from = X, values_from = value) %>% 
    dplyr::rename(DistanceToSummit = "NA_Summit", Signal = "nal_NA") %>% 
    left_join(Motif_Join, by = c("Motif", "dtbs")) %>% 
    group_by(TF, idxATAC) %>% 
    mutate(logFC = log2(Signal/signalBigWig), Type = "Cis")
  
  
  
  
  C <- 
    Motifs_Narrow_TM_CISPB %>% 
    select(idxATAC:last_col()) %>% 
    filter(GeneRna %in% perc_Type_filtered5$Gene) %>%
    pivot_longer(names_to = "Motif", values_to = "value", starts_with("Motif")) %>% 
    filter(!is.na(value)) %>% 
    separate(Motif, into = c("A", "B"), sep = "\\.DistanceTo") %>%
    separate(A, into = c("Motif", "C"), sep = "\\.Sig") %>% 
    unite(X, C, B) %>% 
    mutate(dtbs = "CISBP") %>%
    pivot_wider(names_from = X, values_from = value) %>% 
    dplyr::rename(DistanceToSummit = "NA_Summit", Signal = "nal_NA") %>% 
    left_join(Motif_Join, by = c("Motif", "dtbs")) %>% 
    group_by(TF, idxATAC) %>% 
    mutate(logFC = log2(Signal/signalBigWig), Type = "Trans")
  
  D <- 
    Motifs_Narrow_CF_CISPB %>% 
    select(idxATAC:last_col()) %>% 
    filter(GeneRna %in% perc_Type_filtered5$Gene) %>%
    pivot_longer(names_to = "Motif", values_to = "value", starts_with("Motif")) %>% 
    filter(!is.na(value)) %>% 
    separate(Motif, into = c("A", "B"), sep = "\\.DistanceTo") %>%
    separate(A, into = c("Motif", "C"), sep = "\\.Sig") %>% 
    unite(X, C, B) %>% 
    mutate(dtbs = "CISBP") %>%
    pivot_wider(names_from = X, values_from = value) %>% 
    dplyr::rename(DistanceToSummit = "NA_Summit", Signal = "nal_NA") %>% 
    left_join(Motif_Join, by = c("Motif", "dtbs")) %>% 
    group_by(TF, idxATAC) %>% 
    mutate(logFC = log2(Signal/signalBigWig), Type = "Cis")
  
  Motifs_Narrow <<- bind_rows(A, B, C, D)
  
  rm(A, B, C, D)
  
}


celltype <- "LUM_HR-pos"

Treatment_Response <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Seurat_Objects/Harmony/Treatment_Response_by_CellType.rds") %>% filter(Cluster == celltype, p_val_adj <= 0.05) ### treatment response DE genes
Motif_Join <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/ATAC_Integration/Motif_Summary_Complete_Unfiltered.rds") %>% filter(CellType == celltype) %>% select(Motif, TF, avg_logFC, dtbs)

LUMexpr            <- do.call(rbind.data.frame, readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Seurat_Objects/Reboot/Expression_list.rds")) %>% 
  ungroup() %>% filter(expression == "expressed", CellType == celltype) ### percentage of cell expressing all genes

cluster.averages   <- readRDS("Output/Reboot/TM-CF_cluster_averages_across_CellTypes.rds") %>% select(-avg_logFC, -p_val_adj) %>% 
  filter(Cluster == celltype) %>% left_join(Treatment_Response, by = c("Gene")) %>% select(-p_val, -pct.1, -pct.2) ### average gene expression across TM/CF



### build a table that combines percent expression across cells and average expression and potential fold change
perc_Type <- 
  LUMexpr %>% 
  dplyr::select(1, Gene = 2, 6) %>% 
  pivot_wider(names_from = Type, values_from = perc) %>% 
  dplyr::select(1, perc_CF = 2, perc_TM = 3) %>% 
  left_join(cluster.averages) %>% 
  dplyr::select(1:3, avg_CF = "CF", avg_TM = "TM", 7, 8)

rm(LUMexpr, cluster.averages)

perc_Type_filtered5 <- perc_Type %>% filter(perc_CF >= 5 | perc_TM >= 5) 




Motifs_Narrow_TM_CISPB  <- read_tsv(paste0("ATAC_Integration/NarrowPeak_Files/", celltype, "_Trans_Motif.cisBP_narrowAnnotatedPeakFile.tsv"))
Motifs_Narrow_CF_CISPB  <- read_tsv(paste0("ATAC_Integration/NarrowPeak_Files/", celltype, "_Cis_Motif.cisBP_narrowAnnotatedPeakFile.tsv"))
Motifs_Narrow_TM_JASPAR <- read_tsv(paste0("ATAC_Integration/NarrowPeak_Files/", celltype, "_Trans_Motif_narrowAnnotatedPeakFile.tsv"))
Motifs_Narrow_CF_JASPAR <- read_tsv(paste0("ATAC_Integration/NarrowPeak_Files/", celltype, "_Cis_Motif_narrowAnnotatedPeakFile.tsv"))


make_narrow()

rm(Motifs_Narrow_CF_CISPB, Motifs_Narrow_TM_CISPB, Motifs_Narrow_CF_JASPAR, Motifs_Narrow_TM_JASPAR, Motif_Join)





Motifs_HiCon <-
  Motifs_Narrow %>% ungroup() %>% select(-Is.DE.in.celltype) %>%
  left_join(Treatment_Response, by = c("GeneRna" = "Gene")) %>% 
  select(-p_val, - pct.1, -pct.2, -Cluster) %>%
  unite(Peak, GeneRna, idxATAC, remove = F) %>%
  filter(logFC >= 0.25) %>% 
  filter(DistanceToSummit <= 500, DistanceToSummit >= -500) %>% 
  dplyr::rename(TF.log_FC = "avg_logFC.x", avg_logFC = "avg_logFC.y")


Motifs_Clean <- 
  Motifs_HiCon %>% 
  select(-Motif, -dtbs, -DistanceToSummit, -Signal, -logFC, -FDR, -Correlation, -signalBigWig) %>% 
  mutate(True = "TRUE") %>% unique() %>%
  pivot_wider(names_from = Type, values_from = True) %>% 
  unite(Type, Trans, Cis) %>%
  mutate(Type = recode(Type, "TRUE_TRUE" = "TM+CF", "TRUE_NA" = "TM_Only", "NA_TRUE" = "CF_Only")) %>%
  group_by(TF, GeneRna) %>%
  mutate(peak.TF = n())



#mat <- 
  Motifs_Clean %>% ungroup() %>% filter(GeneRna %in% c("AREG", "EREG")) %>% 
  select(GeneRna, TF, avg_logFC, Type, peak.TF) %>% mutate(occur = 1) %>%
  pivot_wider(names_from = TF, values_from = occur) %>% mutate(sum = AR + PGR) %>% 
    filter(sum == 2) %>% arrange(Peak) %>% print(n = 500)
    #column_to_rownames("Peak")




Motifs_Narrow %>% ungroup() %>% select(-Is.DE.in.celltype) %>%
  left_join(Treatment_Response, by = c("GeneRna" = "Gene")) %>% 
  select(-p_val, - pct.1, -pct.2, -Cluster) %>%
  unite(Peak, GeneRna, idxATAC, remove = F) %>%
  filter(TF %in% c("AR", "PGR", "KLF6", "BATF", "JUN", "CUX2")) %>%
  #filter(logFC >= 0.25) %>% 
  #filter(DistanceToSummit <= 500, DistanceToSummit >= -500) %>%
  #filter(!is.na(avg_logFC)) %>%
  
  ggplot(aes(x = logFC, y = DistanceToSummit)) + 
  geom_point(alpha = 0.25) + 
  geom_hline(yintercept = 0, color = "darkgreen", size = 1, alpha = 0.5) + 
  geom_vline(xintercept = 0, color = "darkgreen", size = 1, alpha = 0.5) #+ 
  facet_wrap(vars(TF))





Test %>% ungroup() %>% filter(!is.na(avg_logFC)) %>%
  group_by(TF, GeneRna) %>% mutate(n = n()) %>% 
  filter(TF == "ZBTB7A") %>% 
  ggplot(aes(x = avg_logFC, y = n, color = Type)) + 
  geom_jitter(size = 2) 






mat <-
  Motifs_HiCon %>% #filter(TF %in% c("AR", "PGR", "NR3C1", "NR3C2")) %>%
  select(Peak, TF, TF.log_FC) %>% 
  #filter(!is.na(TF.log_FC)) %>% 
  distinct() %>%
  mutate(occup = 1) %>% 
  ungroup() %>% 
  select(TF, occup, Peak) %>% 
  distinct() %>% 
  pivot_wider(names_from = TF, values_from = occup) %>% 
  column_to_rownames("Peak")


mat[is.na(mat)] <- 0


pheatmap(t(mat[,c("AR", "PGR", "NR3C1", "NR3C2")]), 
         color = viridis::magma(n = 100, begin = 0.2, end = 0.9), 
         fontsize_col = 0.00001,
         treeheight_col = 0)



res1 <- cor(mat)
round(res1, 2)

tfs <- c("AR", "PGR", "NR3C1", "FOXA1", "FOXN3", "JUN", "JUND", "JUNB", "FOSL2", "SMARCC1")

sel <- res1 %>% as.data.frame() %>% 
  rownames_to_column("A") %>% 
  pivot_longer(colnames(res1)) %>%
  filter(value < 1) %>% 
  group_by(A) %>% 
  mutate(max = max(value)) %>%
  filter(max > 0.3) %>% 
  pull(A) %>% 
  unique()





pheatmap(res1[tfs, tfs], 
         border_color = NA ,
         color = viridis::magma(n = 100, end = 0.9), 
         main = celltype, 
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation", 
         clustering_method = "ward.D2", 
         treeheight_row = 0,
         treeheight_col = 20,
         fontsize_row = 20, fontsize_col = 20,
         ) %>% 
  as.ggplot()






library(VennDiagram)


set1 <- AR_Targets
set2 <- PGR_Targets


venn.diagram(
  x = list(set1, set2),
  category.names = c("AR" , "PGR"),
  filename = "Venn.png",
  output=F
)






