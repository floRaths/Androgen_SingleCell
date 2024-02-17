library(gprofiler2)

Marks  <-  readRDS("Submission/utilities/DE_Test_Everything_TMvCF_20kDownSample.rds")

query        <- filter(Marks, avg_logFC > 0) %>% rownames()
Results      <- gost(query, organism = "hsapiens", evcodes = T)
query        <- filter(Marks, avg_logFC < -0) %>% rownames()
Result_Down <- gost(query, organism = "hsapiens", evcodes = T)


a <- Results$result[ ,1:14]     %>% filter(source == "REAC") %>% select(-1, -2) %>% mutate(Figure = "S1f - upregulated")
b <- Results$result[ ,1:14]     %>% filter(source != "REAC") %>% select(-1, -2) %>% mutate(Figure = "S1f - upregulated")
c <- Result_Down$result[ ,1:14] %>% filter(source == "REAC") %>% select(-1, -2) %>% mutate(Figure = "S1f - downregulated")
d <- Result_Down$result[ ,1:14] %>% filter(source != "REAC") %>% select(-1, -2) %>% mutate(Figure = "S1f - downregulated")

data_1 <- bind_rows(a, b, c, d)




percentile <- 
  read_tsv(paste0("Output/Reboot/SCENIC/expr_mat.adjacencies_all_cells.tsv")) %>% 
  group_by(TF) %>% 
  mutate(percentile_rank = ntile(importance, 100))

clust = "PPARG"

query <- percentile %>% filter(TF == clust, percentile_rank > 95) %>% pull(target)

Results_PPARG <- gost(query, organism = "hsapiens")

a <- Results_PPARG$result[ ,1:14] %>% filter(source == "GO:BP") %>% select(-1, -2) %>% mutate(Figure = "S5f - PPARG module")
b <- Results_PPARG$result[ ,1:14] %>% filter(source != "GO:BP") %>% select(-1, -2) %>% mutate(Figure = "S5f - PPARG module")

data_2 <- bind_rows(a, b)




Treatment_Response      <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Seurat_Objects/Harmony/Treatment_Response_by_CellType.rds") ### treatment response DE genes
cabello_aguilar_LR_PMID <- read_csv("~/Downloads/cabello_aguilar_LR_PMID.csv")
cabello <- cabello_aguilar_LR_PMID %>% filter(!is.na(PMIDs))

ligands_de <- Treatment_Response %>% filter(Gene %in% cabello$ligand) %>% pull(Gene) %>% unique()
receptors <- cabello %>% filter(ligand %in% ligands_de) %>% pull(receptor)

receptors_de <- Treatment_Response %>% filter(Gene %in% receptors) %>% pull(Gene) %>% unique()

query <- c(ligands_de, receptors_de)


library(gprofiler2)
library(viridis)

clust = "Combined Signaling"
Results_SIGNAL <- gost(query, organism = "hsapiens", evcodes = TRUE)

a <- Results_SIGNAL$result[ ,1:14] %>% filter(source == "KEGG") %>% select(-1, -2) %>% mutate(Figure = "S7a - Differential Signaling")
b <- Results_SIGNAL$result[ ,1:14] %>% filter(source != "KEGG") %>% select(-1, -2) %>% mutate(Figure = "S7a - Differential Signaling")

data_3 <- bind_rows(a, b)








percentile <- 
  read_tsv(paste0("Output/Reboot/SCENIC/expr_mat.adjacencies_all_cells.tsv")) %>% 
  group_by(TF) %>% 
  mutate(percentile_rank = ntile(importance, 100))

clust = "GPAM"

query <- percentile %>% filter(TF == clust, percentile_rank > 95) %>% pull(target)

Results_GPAM <- gost(query, organism = "hsapiens")

a <- Results_GPAM$result[ ,1:14] %>% filter(source == "WP") %>% select(-1, -2) %>% mutate(Figure = "S7f - GPAM module")
b <- Results_GPAM$result[ ,1:14] %>% filter(source != "WP") %>% select(-1, -2) %>% mutate(Figure = "S7f - GPAM module")

data_4 <- bind_rows(a, b)







library(xlsx)
write.xlsx(select(data_1, 13, 1:11), file = "Submission/Figures_and_Legends/Supp_Tables/Gprofiler2_FullResults.xlsx", sheetName = "Figure - S1f", row.names = F)
write.xlsx(select(data_2, 13, 1:11), file = "Submission/Figures_and_Legends/Supp_Tables/Gprofiler2_FullResults.xlsx", sheetName = "Figure - S5f", row.names = F, append = T)
write.xlsx(select(data_3, 13, 1:11), file = "Submission/Figures_and_Legends/Supp_Tables/Gprofiler2_FullResults.xlsx", sheetName = "Figure - S7a", row.names = F, append = T)
write.xlsx(select(data_4, 13, 1:11), file = "Submission/Figures_and_Legends/Supp_Tables/Gprofiler2_FullResults.xlsx", sheetName = "Figure - S7f", row.names = F, append = T)










Bind %>% ungroup() %>%
  mutate(factor = group1 * FC) %>% 
  filter(group1 > 0.1, FC > 0, Wilcox <= 0.05,
         str_detect(Module, "REACTOME|BIOCARTA|WP_|PID|KEGG|HALLMARK")) %>% 
  
  select(Subcluster = "CLuster", Module, Avg.Score = "group1", logFC = "FC", Wilcox) %>% 
  arrange(Subcluster, Wilcox) %>% left_join(gene_counts)





