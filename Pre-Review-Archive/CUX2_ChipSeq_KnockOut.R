
library(readxl)

CUX2_Adeno <- read_excel("External_Data/CUX2_Adeno-siRNA.xlsx", sheet = "CUX2-Adeno", skip = 5, )

colnames(CUX2_Adeno) <- gsub(" ", "_", colnames(CUX2_Adeno))
colnames(CUX2_Adeno) <- gsub("-", "_", colnames(CUX2_Adeno))

CUX2_Adeno <- CUX2_Adeno %>% dplyr::select(Gene = 4, 9, 10) 
CUX2_Adeno <- CUX2_Adeno %>% dplyr::select(Gene = 4, 10) 

CUX2_Adeno$Gene <- toupper(CUX2_Adeno$Gene)



CUX2_ChIP_Seq <- read_excel("External_Data/CUX2_ChIP-Seq.xlsx", 
                            sheet = "CUX2-ChIP", skip = 4)

colnames(CUX2_ChIP_Seq) <- gsub(" ", "_", colnames(CUX2_ChIP_Seq))
colnames(CUX2_ChIP_Seq) <- gsub("-", "_", colnames(CUX2_ChIP_Seq))

CUX2_ChIP_Seq$CUX2_Target_Gene_Symbol <- toupper(CUX2_ChIP_Seq$CUX2_Target_Gene_Symbol)


CHIP_Overlap <- 
CUX2_ChIP_Seq %>% dplyr::select(Gene = 33, 36) %>% unique() %>% inner_join(LUM_pos_Treatment_Response) %>% arrange(-avg_logFC)







test <- 
  c("PLP1", "AMACR", "ACSL3", "SC4MOL", "HMGCS2", "IDI1", "FASN", "SCD", "ELOVL7", "MBOAT2", "GNE", "CYP1A1", "SC5D", "NSDHL", "ETNK1", "HMGCS1", "EBP", "FADS1", "ACSS2", "LSS", "NANS", "FDPS", "OLAH", "CHKA", "DHCR7", "DHCR24", "LPCAT3", "DPM3", "CMAS", "CDS1", "ELOVL5", "FADS3")

test <- 
c("KLK3", "MSMB", "TMEPAI", "CUX2", "LOC131368", "SCRG1", "MGC29937", "NUDT10", "TRPM8", "CACNG4", "ADAMTS2", "CSE-C", "EIF6", "CDC45L", "GOLGA3", "ACPP", "TMEM59", "SORD", "SCGB1A1", "CYP4F2", "SMOC1", "GMPR", "GUCY2D", "MYBPC1", "MGC40574", "GMCL1L", "SHAX3", "FLJ36004", "SLC5A10", "RAB26", "NDRG1", "MED28", "PRO1848", "KIAA1661", "PDS5A", "TTLL1", "SCGB3A1", "GPR157", "CREB3L4", "IFI44L", "SH3GL3", "RNASEL")






CUX2_Adeno$Effect_Broad <- 
  recode(CUX2_Adeno$Effect_of_Adeno_Cux2, 
       "induced in female" = "induced",
        "induced in male" = "induced",
        "induced in male & female" = "induced",
        "induced in male, suppressed in female" = "induced",
        "suppressed in female" = "suppressed",
        "suppressed in male" = "suppressed",
        "suppressed in male & female" = "suppressed")



mat <- 
CUX2_Adeno %>% dplyr::select(1,3) %>% unique() %>%
  inner_join(LUM_pos_Treatment_Response) %>%
  dplyr::select(1,2,4) %>% 
  pivot_wider(names_from = "Effect_of_Adeno_Cux2", values_from = "avg_logFC")



mat[is.na(mat)] <- 0


mat <- mat %>% column_to_rownames("Gene")
pheatmap(mat)
pheatmap(t(mat), fontsize_row = 13,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation", annotation_col = Ann)


#image <- 
  pheatmap(mat, scale = "column", fontsize = 10, 
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "euclidean",
           #color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)),
           #annotation_col = Ann,
           #annotation_colors = ann_colors, 
           treeheight_row = 0, treeheight_col = 0, 
           cutree_rows = 1, 
           cutree_cols = 10,
           fontsize_row = 9,
           fontsize_col = 6, 
           labels_col = NULL
  )


  c("AHNAK", "CHPT1", "EGFR", "MBNL2", "MPP7", "NCAM2", "PARD3B", "THRB")

  Ann <- CUX2_Adeno %>% dplyr::select(1,2, 3) %>% unique() %>%
    +   inner_join(LUM_pos_Treatment_Response) %>% dplyr::select(1,2) %>% unique() %>% filter(Gene %notin% c("AHNAK", "CHPT1", "EGFR", "MBNL2", "MPP7", "NCAM2", "PARD3B", "THRB")) %>% column_to_rownames("Gene")
  
  
  
  
  
  
  
  
  

  
  
de <- clusterProfiler::bitr(filter(CellType_Response_List, Cluster == "Adipocyte", p_val_adj <=0.05)$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)



filter(CellType_Response_List, Cluster == "Adipocyte", p_val_adj <=0.05) %>% inner_join(de, by = c("Gene" = "SYMBOL")) %>%
  mutate(test = ifelse(avg_logFC > 0, "#f8fc03", "#ddccff")) %>% 
  filter(p_val_adj <= 0.05) %>% 
  dplyr::select("#hsa" = 8,9) %>% write_csv("test.csv")






ARinteract <- 
ARinteract %>% 
  separate(Coregulator, into = c("A", "B")) %>% 
  pivot_longer(cols = c("A", "B")) %>% 
  dplyr::select(-2) %>% 
  filter(!is.na(value))
