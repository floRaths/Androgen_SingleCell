library(tidygraph)
library(ggraph)
library(scales)


# on filtered ATAC data ---------------------------------------------------


#####################################

create_graph <- function(x, y) {
  
  TFs <- read_tsv("~/Documents/SCENIC/pySCENIC_Full/scenicdata/Resources/hs_hgnc_tfs.txt", col_names = "Gene") %>% mutate(TF = "YES")
  
  tfs <- x %>% distinct(TF)          %>% rename(label = TF)
  tgs <- x %>% distinct(Target.Gene) %>% rename(label = Target.Gene)
  pr1 <- y %>% distinct(protein_1)   %>% rename(label = protein_1)
  pr2 <- y %>% distinct(protein_2)   %>% rename(label = protein_2)
  #int <- x %>% distinct(Intrct)      %>% rename(label = Intrct)
  
  nodes <<- 
    full_join(tfs, tgs, by = "label") %>% 
    full_join(pr1, by = "label") %>%
    full_join(pr2, by = "label") %>%
    rowid_to_column("id") %>% 
    left_join(Treatment_Response, by = c("label" = "Gene")) %>% 
    left_join(TFs, by = c("label" = "Gene")) %>% 
    select(1,2,4,9) %>% 
    mutate(TF = ifelse(is.na(TF), "Target", "TF"), 
           avg_logFC = replace_na(avg_logFC, "0"), 
           avg_logFC = as.numeric(avg_logFC))
  
  
  

  per_route1 <- x %>% select(TF, 
                            target = "Target.Gene", 
                            ATAC.Peaks = "peak.TF") %>% 
    unique() %>% mutate(ATAC.Peaks = rescale(ATAC.Peaks, to = c(1,100)))
  
  
  per_route2 <- 
    y %>% select(protein_1, protein_2, combined_score) %>% 
    unique() %>% 
    group_by(protein_1) %>% 
    mutate(combined_score = rescale(combined_score, to = c(1,35))) %>% 
    ungroup()
  
  
  
  ##########################
  
  edges <- per_route1 %>% 
    left_join(nodes, by = c("TF" = "label")) %>% 
    rename(from = id)
  
  
  edges <- edges %>% 
    left_join(nodes, by = c("target" = "label")) %>% 
    rename(to = id)
  
  
  edges1 <- select(edges, from, to, weight = ATAC.Peaks) %>% mutate(Source = "ATAC")
  
  ##########################
  
  edges <- per_route2 %>% 
    left_join(nodes, by = c("protein_1" = "label")) %>% 
    rename(from = id)
  
  
  edges <- edges %>% 
    left_join(nodes, by = c("protein_2" = "label")) %>% 
    rename(to = id)
  
  
  edges2 <- select(edges, from, to, weight = combined_score) %>% mutate(Source = "STRINGII")
  
  edges <- bind_rows(edges1,edges2)
  
  ##########################
  
  routes_tidy <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  
  
  routes_tidy %>% 
    activate(edges) %>% 
    arrange(desc(weight))
  
}


#diff_TFs <- TF_Target_Motif %>% filter(!is.na(TF.log_FC)) %>% pull(TF) %>% unique()

query_TF   <- c()

add_TFs <- Motifs_Acc %>% filter(TF %in% query_TF, Target.Log2FC < 0) %>% ungroup() %>% unique()


#x <- TF_Target_Motif %>% filter(peak.TF > 0) %>% 
#  group_by(TF) %>% 
#  distinct(Target.Gene, peak.TF) %>% 
#  top_n(5, peak.TF) %>% 
#  ungroup() %>% unique()

#Motifs_inter <- Motifs_Acc %>% filter(Target.Gene %in% c("JUNB", "JUN", "JUND", "YY1", "ZBTB7A"), TF %notin% c("CTCF"), !is.na(Target.Log2FC))



limit = 0.25
min.mot = 5

x <- Motifs_inter %>% 
  filter(Target.Log2FC > limit | Target.Log2FC < -limit) %>% 
  filter(TF.log_FC > limit | TF.log_FC < -limit) %>% 
  bind_rows(add_TFs) %>%
  group_by(TF, Target.Gene) %>% 
  mutate(peak.TF = n()) %>% 
  filter(peak.TF > min.mot) %>% 
  group_by(TF) %>% 
  distinct(Target.Gene, peak.TF) %>% 
  top_n(50, peak.TF) %>% 
  ungroup() %>% unique()



y <- interact_select %>% group_by(protein_1) %>% top_n(20, combined_score) %>%  
  filter(protein_1 %in% c(), #%>% ungroup() %>% unique() 
                                         TF_2 == "TRUE" | !is.na(pro2.Log2FC )) %>% ungroup() %>% unique()

#"PIK3R1", "FN1", "IGFBP5", "CACNA1D", "MYBPC1", "FHL1", "EGFL7", "EGFR", "ERBB4"


routes_tidy <- create_graph(x, y)




ggraph(routes_tidy, layout = "stress") + 
  #geom_edge_link   (aes(width  = weight), alpha = 0.5, color = "gray60") +
  #geom_edge_link   (aes(width  = Interaction), alpha = 0.5, color = "powderblue") +
  geom_edge_fan2   (aes(color = Source, width = weight, linetype = Source), alpha = 0.5, arrow = arrow(angle = 10, ends = "last", type = "closed")) +
  scale_edge_width (range = c(0.2, 2)) +
  
  geom_node_point  (aes(colour = avg_logFC, 
                        size = avg_logFC, shape = TF), cex = 7) + 
  
  scale_color_gradient2(low = "orange", mid = "grey70", high = "purple") +
  scale_edge_colour_manual(values = c("gray60", "lightseagreen")) +
  
  geom_node_text (aes(label = label), 
                 size = 4,
                 repel = TRUE,
                 color = "black") + 
  #labs(edge_width = "motifs") +
  scale_size(range = c(1,6))+
  theme_graph(background = "white", 
              base_family = "Lato Bold") + 
  #ggtitle("ATAC_FC, Target_FC and TF_FC must agree (AR exempt) ... min. 1 peaks/connection")#
  #ggtitle(paste(celltype, sep=", ", collapse=", ") )
  ggtitle(label =  paste(celltype, " - Connections on accessible regions. FC of RNA and TF, must agree"),
          subtitle = paste("exempt TFs:", paste(query_TF, 
                                                paste("  min. target FC =  +/-", limit),
                                                paste("  min.motifs =", min.mot),
                                 
                                 sep = ", ", collapse = ", ")))






#-------------------------------------------------


make_z = function(x){ return ( (x - mean(x)) / sd(x) ) }



Motifs_Pre <- do.call(cbind.data.frame, readRDS("ATAC_Integration/peakTsvFiles/allCellTypes_QuantitativeMotifEnrichment_noFiltering_signalScore.RDS")@listData) %>% rowid_to_column()

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





dplyr::select(1, GeneRna, name, Scaled) %>%
  unite(Gene_Peak, 2, 1) %>%
  pivot_wider(names_from = name, values_from = Scaled) %>%
  column_to_rownames("Gene_Peak")



pheatmap(t(sample_n(mat, 1000)), labels_col = F, 
         color = viridis::viridis(n = 100))


test <- 
read_tsv("ATAC_Integration/peakTsvFiles/Fibroblast_Cis_female.tsv") %>% 
  full_join(read_tsv("ATAC_Integration/peakTsvFiles/Fibroblast_Trans_male.tsv"), 
            by = c("seqnames", "idx", "start", "end")) %>% 
  dplyr::select(-Log2FC.x, -FDR.x, -MeanDiff.x, -Log2FC.y, -FDR.y, -MeanDiff.y) %>%
  mutate(Access.Fibroblast = "TRUE")



Motifs <- 
Motifs_Pre %>% 
  dplyr::select(-(12:19)) %>% 
  as_tibble() %>%
  left_join(test, by = c("ATAC.seqnames" = "seqnames", 
                         "ATAC.start" = "start", 
                         "ATAC.end" = "end")) %>% 
  dplyr::relocate(c(idx, Access.Fibroblast), .after = ATAC.end) %>%
  filter(!is.na(Access.Fibroblast)) #%>% print(n = 1000)






build_motif_summary2 <- function() {
  
  # Build Searchable Motif Summary
  names  <- colnames(dplyr::select(Motifs, 14:length(colnames(Motifs))))
  
  ## pull out motif source categories
  
  ar_chip <- Motifs %>% dplyr::select(starts_with("AR_")) %>% colnames()
  enco <- Motifs %>% dplyr::select(starts_with("ENCODE")) %>% colnames()
  cisb <- Motifs %>% dplyr::select(starts_with("Motif.cisBP"))  %>% colnames()
  dbs  <- c(ar_chip, enco, cisb)
  
  jasp <- Motifs %>% dplyr::select(-dbs, -c(1:13)) %>% colnames()
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


#TF_Target_Motif <- 
  Motifs %>% 
  dplyr::select(1,8, unique(Motif_Summary$Motif)) %>% unique() %>%                # only relevant motifs
  pivot_longer(unique(Motif_Summary$Motif), names_to = "Motif") %>%   # long format
  filter(value > 0, GeneRna %in% perc_Type_filtered5$Gene) %>% # only significant peaks and expressed targets        
  dplyr::rename(Target.Gene = GeneRna) %>%                            
  #dplyr::select(Target.Gene, idxATAC, ATAC.Log2FC, Motif) %>%         # keep relevant columns (can be exapnded)
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
