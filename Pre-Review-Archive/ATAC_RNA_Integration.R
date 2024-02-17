library(tidyverse)
setwd("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/")

make_z = function(x){ return ( (x - mean(x)) / sd(x) ) }

celltype <- "Adipocyte"



#Motifs <- do.call(cbind.data.frame, readRDS(paste0("ATAC_Integration/drive-download-20200830T214337Z-001/", celltype, "_MotifAndEncodeAndAREnrichment_noFiltering.RDS"))@listData) %>% dplyr::select(-Comparison, -ATAC.seqnames, -ATAC.start, -ATAC.end, -Is.DE.in.celltype)
Motifs_Pre <- do.call(cbind.data.frame, readRDS("ATAC_Integration/peakTsvFiles/allCellTypes_QuantitativeMotifEnrichment_noFiltering_signalScore.RDS")@listData) %>% rowid_to_column()

### load resources:
Treatment_Response <- do.call(rbind.data.frame, readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Output/Reboot/Master_DE_List_0.1cutoff.rds")) %>% filter(Cluster == celltype, p_val_adj <= 0.05) ### treatment response DE genes

LUMexpr            <- do.call(rbind.data.frame, readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Seurat_Objects/Reboot/Expression_list.rds")) %>% 
                      ungroup() %>% filter(expression == "expressed", CellType == celltype) ### percentage of cell expressing all genes

cluster.averages   <- readRDS("Output/Reboot/TM-CF_cluster_averages_across_CellTypes.rds") %>% select(-avg_logFC, -p_val_adj) %>% 
                      filter(Cluster == celltype) %>% left_join(Treatment_Response, by = c("Gene", "Cluster")) %>% select(-p_val, -pct.1, -pct.2) ### average gene expression across TM/CF



scale_ATAC_signal <- function(Motifs_Pre) {
ToScale <- Motifs_Pre %>% dplyr::select(starts_with("Signal")) %>% colnames()

Scaled_Signals <- 
  Motifs_Pre %>%
  dplyr::select(1:20) %>% 
  #dplyr::select(1:19) %>% 
  #filter(GeneRna == "TP63") %>% 
  pivot_longer(all_of(ToScale)) %>% 
  separate(name, into = c("Signal", "CellType", "Gender"), sep = "\\.") %>%
  group_by(GeneRna) %>% 
  mutate(Plus1 = value + 1) %>%
  mutate(log10 = log10(Plus1)) %>%
  mutate(Scaled = make_z(log10)) %>% 
  ungroup()
}
Scaled_Signals <- scale_ATAC_signal(Motifs_Pre)


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



build_protein_links <- function() {
  
  TFs <- read_tsv("~/Documents/SCENIC/pySCENIC_Full/scenicdata/Resources/hs_hgnc_tfs.txt", col_names = "Gene") %>% mutate(TF = "YES")
  protein_links_full <- read_table2("utilities/9606.protein.links.full.v11.0.txt")
  protein_info       <- read_table2("utilities/9606.protein.info.v11.0.txt")
  
  
  TF_expression <- 
    perc_Type %>% 
    filter(Gene %in% TFs$Gene, 
           perc_CF >= 5  | perc_TM >= 5, 
           avg_CF >= 0.5 | avg_TM >= 0.5) %>% 
    mutate(TF = "TRUE") %>% 
    dplyr::select(Gene, TF)
  
  
  
  protein_links <- 
    protein_links_full %>%
    left_join(protein_info,  by = c("protein1" = "protein_external_id")) %>%
    left_join(protein_info,  by = c("protein2" = "protein_external_id"), suffix = c("_1", "_2")) %>%
    left_join(TF_expression, by = c("preferred_name_1" = "Gene")) %>%
    left_join(TF_expression, by = c("preferred_name_2" = "Gene"), suffix = c("_1", "_2")) %>%
    dplyr::rename(protein_1 = preferred_name_1, protein_2 = preferred_name_2, ENS1 = protein1, ENS2 = protein2) %>%
    relocate(c(protein_1, annotation_1, TF_1, protein_2, annotation_2, TF_2, combined_score), .after = ENS2) %>%
    relocate(c(ENS1, ENS2, protein_size_1, protein_size_2), .after = textmining_transferred) %>%
    mutate(TF_1 = replace_na(TF_1, "FALSE"), TF_2 = replace_na(TF_2, "FALSE")) %>%
    filter(experiments > 0 | experiments_transferred > 0)
  
  rm(protein_links_full, protein_info)
  

  protein_links <- 
    protein_links %>% 
    left_join(Treatment_Response, by = c("protein_1" = "Gene")) %>% 
    dplyr::select(-p_val, -pct.1, -pct.2, -p_val_adj, -Cluster) %>%
    left_join(Treatment_Response, by = c("protein_2" = "Gene")) %>%
    dplyr::select(-p_val, -pct.1, -pct.2, -p_val_adj, -Cluster) %>%
    dplyr::relocate(avg_logFC.x, .after = annotation_1) %>%
    dplyr::relocate(avg_logFC.y, .after = annotation_2) %>%
    dplyr::rename(pro1.Log2FC = avg_logFC.x, pro2.Log2FC = avg_logFC.y)
  
}
protein_links <- build_protein_links()


interact_select <- 
  protein_links %>% 
  filter(protein_1 %in% perc_Type_filtered5$Gene) %>%
  inner_join(perc_Type_filtered5, by = c("protein_2" = "Gene")) %>%
  dplyr::select(1:9)


build_motif_summary2 <- function() {
  
  # Build Searchable Motif Summary
  names  <- colnames(dplyr::select(Motifs_Pre, 21:length(colnames(Motifs_Pre))))
  
  ## pull out motif source categories
  
  ar_chip <- Motifs_Pre %>% dplyr::select(starts_with("AR_")) %>% colnames()
  enco <-    Motifs_Pre %>% dplyr::select(starts_with("ENCODE")) %>% colnames()
  cisb <-    Motifs_Pre %>% dplyr::select(starts_with("Motif.cisBP"))  %>% colnames()
  dbs  <-    c(ar_chip, enco, cisb)
  
  jasp <- Motifs_Pre %>% dplyr::select(-dbs, -c(1:20)) %>% colnames()
  #ar_chip <- Motifs %>% dplyr::select(-dbs, -c(1:15)) %>% colnames()
  
  
  
  Test <- perc_Type %>% dplyr::select("TF" = 1) ### just a gene list
  
  
  ### take all JASPAR motifs, pull the related TFs and then merge with perc_Type to allow filtering
  jasp <- 
    jasp %>% 
    as_tibble() %>% 
    separate(value, into = c("PREFIX", "B", "C", "D", "NUM"), remove = F) %>% 
    left_join(Test, by = c("B" = "TF"), keep = T) %>% 
    left_join(Test, by = c("C" = "TF"), keep = T) %>% 
    left_join(Test, by = c("D" = "TF"), keep = T) %>% 
    dplyr::select(-c(2:6)) %>% 
    filter(!is.na(TF.x) | !is.na(TF.y) | !is.na(TF)) %>% 
    #head(30) %>% 
    pivot_longer(c("TF.x", "TF.y", "TF"), 
                 names_to = "XYZ", 
                 names_repair = "minimal", 
                 values_drop_na = T) %>% 
    dplyr::select(-XYZ, Motif = 1,  TF = 3)  %>%
    inner_join(perc_Type, by = c("TF" = "Gene")) %>% 
    mutate(dtbs = "JASPAR")
  
  
  enco <- 
    enco %>% 
    as_tibble() %>% 
    separate(value, into = c("PREFIX", "B", "C", "D", "NUM"), remove = F) %>% 
    mutate(C = toupper(C), D = toupper(D)) %>% 
    mutate(C = recode(C, "ERALPHA" = "ESR1", "P300" = "EP300")) %>%
    left_join(Test, by = c("NUM" = "TF"), keep = T) %>% 
    left_join(Test, by = c("C" = "TF"), keep = T) %>% 
    left_join(Test, by = c("D" = "TF"), keep = T) %>% 
    dplyr::select(-c(2:6)) %>% 
    filter(!is.na(TF.x) | !is.na(TF.y) | !is.na(TF)) %>% 
    #head(30) %>% 
    pivot_longer(c("TF.x", "TF.y", "TF"), 
                 names_to = "XYZ", 
                 names_repair = "minimal", 
                 values_drop_na = T) %>% 
    dplyr::select(-XYZ, Motif = 1,  TF = 3) %>% 
    inner_join(perc_Type, by = c("TF" = "Gene")) %>% 
    mutate(dtbs = "ENCODE")
  
  
  cisb <- 
    cisb %>% 
    as_tibble() %>% 
    separate(value, into = c("PREFIX", "B", "C", "D", "NUM"), remove = F) %>% 
    left_join(Test, by = c("B" = "TF"), keep = T) %>% 
    left_join(Test, by = c("C" = "TF"), keep = T) %>% 
    left_join(Test, by = c("D" = "TF"), keep = T) %>% 
    dplyr::select(-c(2:6)) %>% 
    filter(!is.na(TF.x) | !is.na(TF.y) | !is.na(TF)) %>% 
    #head(30) %>% 
    pivot_longer(c("TF.x", "TF.y", "TF"), 
                 names_to = "XYZ", 
                 names_repair = "minimal", 
                 values_drop_na = T) %>% 
    dplyr::select(-XYZ, Motif = 1,  TF = 3) %>% 
    inner_join(perc_Type, by = c("TF" = "Gene")) %>% 
    mutate(dtbs = "CISBP")
  
  
  ar_chip <- 
    ar_chip %>% 
    as_tibble() %>% 
    mutate(TF = "AR") %>% 
    dplyr::select(Motif = 1, TF) %>% 
    inner_join(perc_Type, by = c("TF" = "Gene")) %>% 
    mutate(dtbs = "ChIP")
  
  
  Motif_Summary_Unfiltered <- rbind(jasp, enco, cisb, ar_chip) %>% dplyr::rename(TF.log_FC = avg_logFC, TF.p_val_adj = p_val_adj)#%>% saveRDS("ATAC_Integration/Motif_Summary_TFs_and_Expression.rds")
  
  #rm(Test, jasp, enco, cisb, ar_chip, dbs, names)
  
  
  
  # Filter and Create Relevant Peak Motif Table
  
  ### only keep motifs that come from TFs which are expressed in at least 5% of TM or CF cells
  Motif_Summary <- 
    Motif_Summary_Unfiltered %>% filter(perc_CF >= 5 | perc_TM >= 5, avg_CF >= 0.5 | avg_TM >= 0.5)
  
  #rm(Motif_Summary_Unfiltered)
  
}
Motif_Summary <- build_motif_summary2()



#build_motif_summary <- function() {

# Build Searchable Motif Summary
names  <- colnames(dplyr::select(Motifs, 16:length(colnames(Motifs))))

## pull out motif source categories
jasp <- Motifs %>% dplyr::select(starts_with("JASPAR")) %>% colnames()
enco <- Motifs %>% dplyr::select(starts_with("ENCODE")) %>% colnames()
cisb <- Motifs %>% dplyr::select(starts_with("CisBP"))  %>% colnames()
dbs  <- c(jasp, enco, cisb)

ar_chip <- Motifs %>% dplyr::select(-dbs, -c(1:15)) %>% colnames()



Test <- perc_Type %>% dplyr::select("TF" = 1) ### just a gene list


### take all JASPAR motifs, pull the related TFs and then merge with perc_Type to allow filtering
jasp <- 
  jasp %>% 
  as_tibble() %>% 
  separate(value, into = c("PREFIX", "B", "C", "D", "NUM"), remove = F) %>% 
  left_join(Test, by = c("B" = "TF"), keep = T) %>% 
  left_join(Test, by = c("C" = "TF"), keep = T) %>% 
  left_join(Test, by = c("D" = "TF"), keep = T) %>% 
  dplyr::select(-c(2:6)) %>% 
  filter(!is.na(TF.x) | !is.na(TF.y) | !is.na(TF)) %>% 
  #head(30) %>% 
  pivot_longer(c("TF.x", "TF.y", "TF"), 
               names_to = "XYZ", 
               names_repair = "minimal", 
               values_drop_na = T) %>% 
    dplyr::select(-XYZ, Motif = 1,  TF = 3)  %>%
  inner_join(perc_Type, by = c("TF" = "Gene")) %>% 
  mutate(dtbs = "JASPAR")


enco <- 
  enco %>% 
  as_tibble() %>% 
  separate(value, into = c("PREFIX", "B", "C", "D", "NUM"), remove = F) %>% 
  mutate(C = toupper(C), D = toupper(D)) %>% 
  mutate(C = recode(C, "ERALPHA" = "ESR1", "P300" = "EP300")) %>%
  left_join(Test, by = c("NUM" = "TF"), keep = T) %>% 
  left_join(Test, by = c("C" = "TF"), keep = T) %>% 
  left_join(Test, by = c("D" = "TF"), keep = T) %>% 
  dplyr::select(-c(2:6)) %>% 
  filter(!is.na(TF.x) | !is.na(TF.y) | !is.na(TF)) %>% 
  #head(30) %>% 
  pivot_longer(c("TF.x", "TF.y", "TF"), 
               names_to = "XYZ", 
               names_repair = "minimal", 
               values_drop_na = T) %>% 
  dplyr::select(-XYZ, Motif = 1,  TF = 3) %>% 
  inner_join(perc_Type, by = c("TF" = "Gene")) %>% 
  mutate(dtbs = "ENCODE")


cisb <- 
  cisb %>% 
  as_tibble() %>% 
  separate(value, into = c("PREFIX", "B", "C", "D", "NUM"), remove = F) %>% 
  left_join(Test, by = c("B" = "TF"), keep = T) %>% 
  left_join(Test, by = c("C" = "TF"), keep = T) %>% 
  left_join(Test, by = c("D" = "TF"), keep = T) %>% 
  dplyr::select(-c(2:6)) %>% 
  filter(!is.na(TF.x) | !is.na(TF.y) | !is.na(TF)) %>% 
  #head(30) %>% 
  pivot_longer(c("TF.x", "TF.y", "TF"), 
               names_to = "XYZ", 
               names_repair = "minimal", 
               values_drop_na = T) %>% 
  dplyr::select(-XYZ, Motif = 1,  TF = 3) %>% 
  inner_join(perc_Type, by = c("TF" = "Gene")) %>% 
  mutate(dtbs = "CISBP")


ar_chip <- 
  ar_chip %>% 
  as_tibble() %>% 
  mutate(TF = "AR") %>% 
  dplyr::select(Motif = 1, TF) %>% 
  inner_join(perc_Type, by = c("TF" = "Gene")) %>% 
  mutate(dtbs = "ChIP")


Motif_Summary_Unfiltered <- rbind(jasp, enco, cisb, ar_chip) %>% dplyr::rename(TF.log_FC = avg_logFC, TF.p_val_adj = p_val_adj)#%>% saveRDS("ATAC_Integration/Motif_Summary_TFs_and_Expression.rds")

#rm(Test, jasp, enco, cisb, ar_chip, dbs, names)



# Filter and Create Relevant Peak Motif Table

### only keep motifs that come from TFs which are expressed in at least 5% of TM or CF cells
Motif_Summary <- 
  Motif_Summary_Unfiltered %>% filter(perc_CF >= 5 | perc_TM >= 5, avg_CF >= 0.5 | avg_TM >= 0.5)

#rm(Motif_Summary_Unfiltered)

}
#Motif_Summary <- build_motif_summary()




Motifs_Acc <-   
  Scaled_Signals %>%
  filter(Scaled > 0, CellType == celltype) %>%
  dplyr::select(rowid, GeneRna, CellType) %>%
  left_join(dplyr::select(Motifs_Pre,
                          idxATAC,
                          rowid, 
                          GeneRna, 
                          unique(Motif_Summary$Motif)), 
            by = c("rowid", "GeneRna")) %>%
  pivot_longer(unique(Motif_Summary$Motif), names_to = "Motif") %>%
  filter(value > 0, GeneRna %in% perc_Type_filtered5$Gene) %>%
  dplyr::rename(Target.Gene = GeneRna) %>%                            
  #dplyr::select(Target.Gene, idxATAC, ATAC.Log2FC, Motif) %>%         # keep relevant columns (can be exapnded)
  left_join(Motif_Summary, by = "Motif") %>%
  dplyr::select(-perc_CF, -perc_TM, -avg_CF, -avg_TM, -TF.p_val_adj) %>%
  left_join(Treatment_Response, by = c("Target.Gene" = "Gene")) %>%
  dplyr::rename(Target.Log2FC = avg_logFC) %>%
  dplyr::select(-p_val, -pct.1, -pct.2, -p_val_adj, -Cluster) %>%
  dplyr::relocate(c(Target.Gene, Target.Log2FC, Motif, idxATAC, value, dtbs), 
                  .after = TF.log_FC) %>% unique() %>%
  #filter(!is.na(TF.log_FC), !is.na(Target.Log2FC)) %>% #### REMOVE!!!!
  group_by(TF, Target.Gene) %>% 
  mutate(peak.TF = n()) 



all_up <- Motifs_Acc %>% filter(Target.Log2FC > 0, TF.log_FC > 0) %>% ungroup() %>% unique()
all_do <- Motifs_Acc %>% filter(Target.Log2FC < 0, TF.log_FC < 0) %>% ungroup() %>% unique()


Motifs_inter <- bind_rows(all_up, all_do)
# list of all genes expressed in 5% of cells <- nothing else considered anywhere



# Take our relevant motifs from above (filered by 5% and min avg_expression)
# keep only significant peaks and expressed genes
# add the TF info from the relevant motif summary to match motifs with TFs
# then also add the DE info for the Target genes
# then calculate peaks/TF in each gene. i.e. TF-X has a motif in 14 peaks of Gene-Y

TF_Target_Motif <- 
  Motifs %>% 
  dplyr::select(1:15, unique(Motif_Summary$Motif)) %>%                # only relevant motifs
  pivot_longer(unique(Motif_Summary$Motif), names_to = "Motif") %>%   # long format
  filter(value == "TRUE", ATAC.FDR <= 0.05, GeneRna %in% perc_Type_filtered5$Gene) %>% # only significant peaks and expressed targets        
  dplyr::rename(Target.Gene = GeneRna) %>%                            
  dplyr::select(Target.Gene, idxATAC, ATAC.Log2FC, Motif) %>%         # keep relevant columns (can be exapnded)
  left_join(Motif_Summary, by = "Motif") %>%
  dplyr::select(-perc_CF, -perc_TM, -avg_CF, -avg_TM, -TF.p_val_adj) %>%
  left_join(Treatment_Response, by = c("Target.Gene" = "Gene")) %>%
  dplyr::rename(Target.Log2FC = avg_logFC) %>%
  dplyr::select(-p_val, -pct.1, -pct.2, -p_val_adj, -Cluster) %>%
  dplyr::relocate(c(Target.Gene, Target.Log2FC, Motif, idxATAC, ATAC.Log2FC, dtbs), 
                  .after = TF.log_FC) %>% 
  #filter(!is.na(TF.log_FC), !is.na(Target.Log2FC)) %>% #### REMOVE!!!!
  group_by(TF, Target.Gene) %>% 
  mutate(peak.TF = n()) 



#protein_select <- protein_links %>% group_by(protein_1) %>% top_n(50, combined_score) # top 25 interaction partners







total_up <- TF_Target_Motif %>% filter(Target.Log2FC > 0, ATAC.Log2FC > 0) %>% ungroup() %>% unique()
total_do <- TF_Target_Motif %>% filter(Target.Log2FC < 0, ATAC.Log2FC < 0) %>% ungroup() %>% unique()

total_up_all <- TF_Target_Motif %>% filter(Target.Log2FC > 0, ATAC.Log2FC > 0, TF.log_FC > 0) %>% ungroup() %>% unique()
total_do_all <- TF_Target_Motif %>% filter(Target.Log2FC < 0, ATAC.Log2FC < 0, TF.log_FC < 0) %>% ungroup() %>% unique()

total_up_all <- TF_Target_Motif %>% filter(ATAC.Log2FC > 0, TF.log_FC > 0) %>% ungroup() %>% unique()
total_do_all <- TF_Target_Motif %>% filter(ATAC.Log2FC < 0, TF.log_FC < 0) %>% ungroup() %>% unique()



Motifs_Agree <- 
  rbind(total_do, total_up) %>% 
  dplyr::select(-peak.TF)


All_Agree <- 
  rbind(total_do_all, total_up_all)# %>% 
  #dplyr::select(-peak.TF) #%>%
  #group_by(TF, Target.Gene) %>% 
  #mutate(peak.TF = n()) %>% 
  #relocate(c(peak.TF), .after = TF)

#rm(total_ar, total_do, total_up)
  


#peak.TF <- 
#  Motifs_inter %>% 
#  dplyr::select(1:8) %>% 
#  unique() %>%
#  group_by(TF, Target.Gene) %>% 
#  mutate(peak.TF = n()) %>% 
#  ungroup() %>%
#  dplyr::select(TF, Target.Gene, peak.TF) %>% 
#  unique()


#Motifs_inter <- 
#  Motifs_inter %>% left_join(peak.TF, by = c("TF", "Target.Gene")) %>% relocate(c(peak.TF), .after = dtbs)


#rm(total_ar, total_do, total_up)














# merge with DE gene list
Motifs_Join <- Motifs %>% 
  dplyr::select(1:16, unique(Motif_Summary$Motif)) %>%
  inner_join(Treatment_Response, by = c("GeneRna" = "Gene")) %>%  
  pivot_longer(unique(Motif_Summary$Motif), names_to = "Motif") %>% 
  inner_join(Motif_Summary, by = "Motif") %>%
  relocate(c(GeneRna, ATAC.Log2FC, avg_logFC, 
             TF.log_FC, Motif, TF, value, 
             perc_CF, perc_TM, avg_CF, 
             avg_TM, dtbs), 
           .after = idxATAC) %>%
  dplyr::rename(Target.Gene = GeneRna, Target.Log2FC = avg_logFC)


### create agreeing peaks and motifs and targets
total_up <- Motifs_Join %>% filter(value == "TRUE", ATAC.FDR <= 0.05, ATAC.Log2FC > 0, Target.Log2FC > 0, TF.log_FC > 0) %>% unique()
total_do <- Motifs_Join %>% filter(value == "TRUE", ATAC.FDR <= 0.05, ATAC.Log2FC < 0, Target.Log2FC < 0, TF.log_FC < 0) %>% unique()
total_ar <- Motifs_Join %>% filter(value == "TRUE", ATAC.FDR <= 0.05, TF == "AR") %>% unique()

Motifs_inter <- 
  rbind(total_ar, total_do, total_up) %>% 
  group_by(TF, Target.Gene) %>% 
  mutate(peak.TF = n()) %>% 
  dplyr::select(-value) %>%
  relocate(c(peak.TF), .after = TF)

rm(total_ar, total_do, total_up)








Motif_Cor <- 
  Motifs_Join %>% 
  filter(value == "TRUE", ATAC.FDR <= 0.05) %>% 
  dplyr::select(1:6) %>% 
  group_by(Motif) %>% 
  mutate(peak.motif = n(), log_cor = cor(ATAC.Log2FC, avg_logFC)) %>% 
  left_join(Motif_Summary, by = "Motif") %>% 
  unique() %>% group_by(TF, GeneRna) %>% mutate(peak.TF = n())


Motif_Cor <- 
  Motifs_Join %>% 
  filter(value == "TRUE", ATAC.FDR <= 0.05) %>% 
  group_by(Motif) %>% 
  mutate(peak.motif = n(), log_cor = cor(ATAC.Log2FC, avg_logFC)) %>% 
  left_join(Motif_Summary, by = "Motif") %>% filter(avg_logFC > 0 & TF) 
  unique() %>% group_by(TF, GeneRna) %>% mutate(peak.TF = n())



Motif_Cor %>% 
  dplyr::select(6:14) %>% 
  unique() %>% 
  filter(peak.motif >= 4, 
         log_cor >= 0.25, 
         avg_CF  >= 0.5 | avg_TM >= 0.5) %>% 
  arrange(-peak.motif)



ar_genes <- Motif_Cor %>% filter(TF == "AR", ATAC.Log2FC > 0, avg_logFC > 0) %>% pull(GeneRna) %>% unique()




### overview of all likely regulators within the DE genes of the celltype
Motif_Cor %>% 
  dplyr::select(6:14) %>% 
  unique() %>% 
  filter(peak.motif >= 4, log_cor >= 0.25, avg_CF >= 0.5 | avg_TM >= 0.5) %>% 
  ggplot(aes(x = reorder(TF, -peak.motif), y = peak.motif)) + 
  geom_point(aes(size = log_cor, alpha = 0.5, stroke = 0)) + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) + 
  geom_text(hjust = -0.25,
            vjust = -0.25, 
            size = 3, 
            aes(label = ifelse(!is.na(TF.log_FC), TF, "")))



### single gene correlation plot
#plot1 <- 
Motifs_Join %>% #filter(str_detect(Motif, "PBX3")) %>%
  #filter(peak.motif >= 4, log_cor >= 0.25, avg_CF >= 0.5 | avg_TM >= 0.5) %>%
  filter(value == "TRUE", ATAC.FDR <= 0.05) %>% dplyr::select(-Motif) %>% unique() %>%
  ggplot(aes(x = avg_logFC, y = ATAC.Log2FC)) + 
  geom_point(alpha = 0.5, stroke = 0, size = 2.5) +#, 
  #aes(color = count, size = count)) +
  stat_smooth(method = "lm")+ 
  geom_text(hjust = -0.25,
            vjust = -0.25, 
            size = 3, 
            aes(label = ifelse(avg_logFC > 0.75 | avg_logFC < -0.75 , GeneRna, ""))) +
  geom_text(hjust = -0.25,
            vjust = -0.25, 
            size = 3, 
            aes(label = ifelse(ATAC.Log2FC > 2.5 | ATAC.Log2FC < -2.5 , GeneRna, ""))) +
  
  ggtitle(label = "Peaks Containing (CisBP_PBX3_518 or JASPAR_PBX3_309) Motif n = 75") + 
  theme(plot.title = element_text(size=15)) + 
  labs(color = "n.motifs", size = "n.motifs") + 
  scale_color_viridis_c(option = "magma", end = 0.9)


plot1 + plot2



### TF regulatory score calculated on all significant peaks
mat <- 
Motif_Cor %>% filter(avg_logFC > 0, peak.motif >= 4, log_cor >= 0.25, avg_CF >= 0.5 | avg_TM >= 0.5) %>% 
  group_by(TF, Motif, GeneRna) %>% 
  arrange(GeneRna) %>% 
  mutate(motif.gene = n()) %>% 
  group_by(TF, GeneRna) %>% 
  mutate(TF.gene = n(), 
         reg.score = ifelse(avg_logFC > 0, TF.gene * avg_TM, TF.gene * avg_CF)) %>%#,
         #reg.score_CF = TF.gene * avg_CF, 
         #reg.score_TM = TF.gene * avg_TM) %>% 
  dplyr::select(GeneRna, TF, reg.score) %>% 
  unique() %>% 
  pivot_wider(names_from = TF, values_from = reg.score) %>%
  column_to_rownames("GeneRna")


mat <- 
Motif_Cor %>% unite(Motif, TF, Motif, remove = F) %>%
  filter(avg_logFC < 0, Motif %notin% ctcf, peak.motif >= 10, log_cor >= 0.25, avg_CF >= 1 | avg_TM >= 1) %>%
  group_by(Motif) %>% filter(perc_CF == max(perc_CF)) %>%
  group_by(Motif, GeneRna) %>% 
  arrange(GeneRna) %>% 
  mutate(motif.gene = n(), 
         reg.score = ifelse(avg_logFC > 0, motif.gene * avg_TM, motif.gene * avg_CF)) %>% 
  
  dplyr::select(GeneRna, Motif, motif.gene) %>% 
  unique() %>% 
  pivot_wider(names_from = Motif, values_from = motif.gene) %>%
  column_to_rownames("GeneRna")




mat <- 
Motif_Cor %>% 
  filter(GeneRna %in% ar_genes, 
         TF %in% relevant_TFs, TF %notin% "CTCF", avg_CF >= 1 | avg_TM >= 1) %>% 
  group_by(GeneRna, TF) %>% 
  mutate(count = n(), n = count * avg_TM) %>% 
  dplyr::select(GeneRna, TF, n) %>% 
  unique() %>% 
  pivot_wider(names_from = TF, values_from = n) %>% 
  column_to_rownames("GeneRna")


mat <- 
percentile %>% 
  filter(target %in% ar_genes, TF %in% ar_tfs) %>% 
  dplyr::select(-percentile_rank) %>%  
  pivot_wider(names_from = TF, values_from = importance) %>%
  column_to_rownames("target")












Ann <- Treatment_Response %>% mutate(Direct = ifelse(avg_logFC > 0, "Up", "Down")) %>% select(Gene, Direct) %>% column_to_rownames("Gene")

mat <- 
  Motifs_Acc %>%
  filter(Target.Log2FC > 0.15 | Target.Log2FC < -0.15, !is.na(TF.log_FC), peak.TF > 0) %>% 
  distinct(TF, Target.Gene, peak.TF) %>% 
  pivot_wider(names_from = Target.Gene, values_from = peak.TF) %>% 
  column_to_rownames("TF")


mat <- 
  Test %>% filter(peak.TF > 0) %>%
  group_by(TF) %>% 
  #mutate(Plus1 = value + 1) %>%
  mutate(log10 = log10(peak.TF)) %>%
  mutate(Scaled = make_z(log10)) %>% 
  ungroup() %>% 
  select(-peak.TF, -log10) %>%
  pivot_wider(names_from = Target.Gene, values_from = Scaled) %>% 
  column_to_rownames("TF")


mat <- 
  Treatment_Response %>% 
  filter(avg_logFC > 0.2 | avg_logFC < -0.2) %>%
  group_by(Gene) %>% 
  mutate(n = n()) %>% 
  filter(n > 2) %>% 
  select(Gene, Cluster, avg_logFC) %>% 
  pivot_wider(names_from = Gene, values_from = avg_logFC) %>% 
  column_to_rownames("Cluster")



Ann <- Motifs_Acc %>% ungroup() %>% select(TF, TF.log_FC) %>% unique() %>% mutate(label = ifelse(TF.log_FC > 0, "UP",  "DOWN")) %>% mutate(label = replace_na(label, "SAME")) %>% column_to_rownames("TF") %>% select(-TF.log_FC)

#xxx <- Acc_Basal %>% ungroup() %>% distinct(TF) %>% sample_n(50) %>% pull(TF)

mat <- 
  Motifs_Acc %>% ungroup() %>% 
  select(-Motif, -dtbs, -rowid, -value) %>% 
  unique() %>% filter(TF.log_FC >= 0.16 | TF.log_FC <= -0.16) %>%
  #filter(TF %in% c(xxx, "TP63", "KLF6")) %>%
  group_by(TF) %>% 
  mutate(n = n()) %>% 
  #filter(n < 4000) %>%
  select(idxATAC, TF) %>% 
  mutate(True = 1) %>% 
  unique() %>% 
  pivot_wider(names_from = TF, values_from = "True") %>% 
  sample_n(6000) %>% column_to_rownames("idxATAC")




mat <- bind_rows(pos, neg, bas, fib) %>% select(-TF.log_FC) %>% filter(peak.TF > 2, TF %notin% "CTCF") %>% pivot_wider(names_from = Cell, values_from = peak.TF) %>% column_to_rownames("TF")


mat[is.na(mat)] <- 0


#image <- 
  pheatmap(mat, 
           scale = "column", 
           fontsize = 10, 
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           #color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)),
           #annotation_row = Ann,
           #annotation_colors = ann_colors, 
           treeheight_row = 0, treeheight_col = 0, 
           #cutree_rows = 1, 
           #cutree_cols = 2,
           fontsize_row = 10,
           fontsize_col = 0.0010, 
           labels_col = NULL, 
           color = viridis::magma(n = 20, begin = 0, end = 0.8), 
           main = "Main"
  )




mat <-   
  cluster.averages %>% 
  filter(Gene %in% c("ZBTB7A", "KLF6")) %>% 
  select(Gene, CF, TM, Cluster) %>% 
  pivot_longer(c("CF", "TM")) %>% 
  unite(Cluster, name, Cluster) %>% 
  pivot_wider(values_from = value, names_from = Cluster) %>% 
  column_to_rownames("Gene")  
  
  
pheatmap(t(mat), 
         scale = "column", 
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "euclidean")




xxx <- 
  Motifs_Acc %>% 
  filter(TF %in% c("KLF6", "ZBTB7A")) %>% 
  select(-Motif, -dtbs, -value) %>% 
  group_by(TF, Target.Gene) %>% 
  unique() %>%
  mutate(peak.TF = n()) %>% filter(peak.TF > 4) %>%
  select(TF, Target.Gene, Target.Log2FC, peak.TF) %>%
  unique() 


cluster.averages %>% filter(Cluster == "LUM_HR-pos", Gene %in% xxx)



### Table with AR-motifs only
AR_Peaks <-
  Motifs %>% 
  dplyr::select(1:16, all_of(dbs)) %>% 
  inner_join(Treatment_Response, by = c("GeneRna" = "Gene")) %>% 
  pivot_longer(dbs, names_to = "Motif") %>% 
  #filter(value == "TRUE") %>% 
  group_by(GeneRna, idxATAC) %>% 
  mutate(count = n()) %>%
  relocate(c(GeneRna, ATAC.Log2FC, avg_logFC, Motif, value, count), .after = idxRNA)




AR_Peaks <-
  Motifs %>% 
  select(1:16, all_of(ar_chip)) %>% 
  inner_join(Treatment_Response, by = c("GeneRna" = "Gene")) %>% 
  pivot_longer(ar_chip, names_to = "Motif") %>% 
  filter(value == "TRUE") %>% 
  group_by(GeneRna, idxATAC) %>% 
  mutate(count = n()) %>%
  relocate(c(GeneRna, ATAC.Log2FC, avg_logFC, Motif, value, count), .after = idxRNA)








Motifs_3 %>% 
  filter(value == "TRUE", ATAC.Pval < 0.05, str_detect(Motif, "JASPAR")) %>% 
  select(-Motif) %>% unique() %>%
  ggplot(aes(x = avg_logFC, y = ATAC.Log2FC, color = ATAC.Pval)) + 
  geom_point() + stat_smooth(method = "lm") + 
  geom_text(hjust = -0.25,
            vjust = -0.25, 
            size = 3, 
            aes(label = ifelse(avg_logFC > 0.75 | avg_logFC < -0.75 , GeneRna, ""))) + 
  geom_hline(yintercept = 0, color = "grey")







Marker_list    <- vector("list", length = length(Vari_Tfs2))

for (i in 1:length(Vari_Tfs2)) {
  
  Marker_list[[i]] <- Motifs_Expr %>% 
    left_join(oeaks_per_genes) %>% 
    filter(str_detect(Motif, query)) %>% 
    select(-Motif, -value) %>% unique() %>% 
    group_by(GeneRna) %>% 
    mutate(count = n()) %>% 
    top_n(50, count) %>%
    ggplot(aes(x = avg_logFC, y = ATAC.Log2FC)) + 
    geom_point(alpha = 0.5, stroke = 0,
               aes(size = count, color = ifelse(ATAC.Pval <= 0.05, "s", "n.s"))) + 
    stat_smooth(method = "lm")+ 
    geom_text(hjust = -0.25,
              vjust = -0.25, 
              size = 3, 
              aes(label = ifelse(avg_logFC > 0.75 | avg_logFC < -0.75 , GeneRna, ""))) +
    geom_text(hjust = -0.25,
              vjust = -0.25, 
              size = 3, 
              aes(label = ifelse(ATAC.Log2FC > 2 | ATAC.Log2FC < -2 , GeneRna, ""))) +
    
    ggtitle(label = query) + 
    theme(plot.title = element_text(size=15)) + 
    labs(color = "Significance") + 
    scale_color_manual(values = c("black", "red"))
  
}



query <- "BATF"

plot4 <- Motifs_Expr %>% filter(GeneRna %in% topCUX, str_detect(Motif, query)) %>% 
  ggplot(aes(x = avg_logFC, y = ATAC.Log2FC)) + 
  geom_point(size = 1, alpha = 0.5, 
             aes(color = ifelse(ATAC.Pval <= 0.05, "s", "n.s"))) + 
  stat_smooth(method = "lm")+ 
  geom_text(hjust = -0.25,
            vjust = -0.25, 
            size = 3, 
            aes(label = ifelse(avg_logFC > 0.1 | avg_logFC < -0.1 , GeneRna, ""))) +
  geom_text(hjust = -0.25,
            vjust = -0.25, 
            size = 3, 
            aes(label = ifelse(ATAC.Log2FC > 2 | ATAC.Log2FC < -2 , GeneRna, ""))) +
  
  ggtitle(label = query) + 
  theme(plot.title = element_text(size=15)) + 
  labs(color = "Significance") + 
  scale_color_manual(values = c("black", "red"))


plot1 + plot2 + plot3 + plot4




query <- "NFIC"

Motifs_Expr %>% 
  left_join(oeaks_per_genes) %>% 
  filter(str_detect(Motif, query)) %>% 
  select(-Motif, -value) %>% unique() %>% 
  group_by(GeneRna) %>% 
  mutate(count = n()) %>% 
  top_n(10, count) %>%
  ggplot(aes(x = avg_logFC, y = ATAC.Log2FC)) + 
  geom_point(alpha = 0.5, stroke = 0,
             aes(size = count, color = ifelse(ATAC.Pval <= 0.05, "s", "n.s"))) + 
  stat_smooth(method = "lm")+ 
  geom_text(hjust = -0.25,
            vjust = -0.25, 
            size = 3, 
            aes(label = ifelse(avg_logFC > 0.75 | avg_logFC < -0.75 , GeneRna, ""))) +
  geom_text(hjust = -0.25,
            vjust = -0.25, 
            size = 3, 
            aes(label = ifelse(ATAC.Log2FC > 2 | ATAC.Log2FC < -2 , GeneRna, ""))) +
  
  ggtitle(label = query) + 
  theme(plot.title = element_text(size=15)) + 
  labs(color = "Significance") + 
  scale_color_manual(values = c("black", "red"))














Marker_list    <- vector("list", length = length(genes))

for (i in 1:length(genes)) {
  Marker_list[[i]] <- Motifs_4 %>% filter(str_detect(Motif, genes[i]))
  
}




build_protein_links <- function() {
  
  TFs <- read_tsv("~/Documents/SCENIC/pySCENIC_Full/scenicdata/Resources/hs_hgnc_tfs.txt", col_names = "Gene") %>% mutate(TF = "YES")
  protein_links_full <- read_table2("utilities/9606.protein.links.full.v11.0.txt")
  protein_info       <- read_table2("utilities/9606.protein.info.v11.0.txt")


TF_expression <- 
  perc_Type %>% 
  filter(Gene %in% TFs$Gene, 
         perc_CF >= 5  | perc_TM >= 5, 
         avg_CF >= 0.5 | avg_TM >= 0.5) %>% 
  mutate(TF = "TRUE") %>% 
  dplyr::select(Gene, TF)



protein_links <- 
  protein_links_full %>%
  left_join(protein_info,  by = c("protein1" = "protein_external_id")) %>%
  left_join(protein_info,  by = c("protein2" = "protein_external_id"), suffix = c("_1", "_2")) %>%
  left_join(TF_expression, by = c("preferred_name_1" = "Gene")) %>%
  left_join(TF_expression, by = c("preferred_name_2" = "Gene"), suffix = c("_1", "_2")) %>%
  dplyr::rename(protein_1 = preferred_name_1, protein_2 = preferred_name_2, ENS1 = protein1, ENS2 = protein2) %>%
  relocate(c(protein_1, annotation_1, TF_1, protein_2, annotation_2, TF_2, combined_score), .after = ENS2) %>%
  relocate(c(ENS1, ENS2, protein_size_1, protein_size_2), .after = textmining_transferred) %>%
  mutate(TF_1 = replace_na(TF_1, "FALSE"), TF_2 = replace_na(TF_2, "FALSE")) %>%
  filter(experiments > 0 | experiments_transferred > 0)
  
rm(protein_links_full, protein_info)



protein_links <- 
  protein_links %>% 
  left_join(Treatment_Response, by = c("protein_1" = "Gene")) %>% 
  dplyr::select(-p_val, -pct.1, -pct.2, -p_val_adj, -Cluster) %>%
  left_join(Treatment_Response, by = c("protein_2" = "Gene")) %>%
  dplyr::select(-p_val, -pct.1, -pct.2, -p_val_adj, -Cluster) %>%
  dplyr::relocate(avg_logFC.x, .after = annotation_1) %>%
  dplyr::relocate(avg_logFC.y, .after = annotation_2) %>%
  dplyr::rename(pro1.Log2FC = avg_logFC.x, pro2.Log2FC = avg_logFC.y)

}


protein_links <- build_protein_links()







