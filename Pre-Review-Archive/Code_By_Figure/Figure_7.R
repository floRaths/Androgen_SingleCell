PI3K_TF_Heatmap <- function() {

Bind %>% saveRDS("Output/PI3K_TF_logFC.rds")
Bind  <-  readRDS("Output/PI3K_TF_logFC.rds")

mat <- 
Bind %>% 
  filter(Gene %in% c(tfs), p_val_adj <= 0.05, abs(avg_logFC) > 0.1) %>% 
  select(Gene, Cluster, avg_logFC) %>% 
  group_by(Gene) %>% 
  mutate(n = n()) %>% 
  filter(n > 4) %>% 
  pivot_wider(names_from = Cluster, 
              values_from = avg_logFC) %>% 
  select(-n) %>% 
  column_to_rownames("Gene")


mat[is.na(mat)] <- 0

my.breaks <- c(seq((min(mat)), -0.00001, length.out = 50), seq(0.00001, (max(mat)), length.out = 50))

p <- 
pheatmap(mat, breaks = my.breaks, fontsize = 20, 
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "euclidean",
         color = colorRampPalette(rev(brewer.pal(11, "BrBG")))(100))
}

p <- PI3K_TF_Heatmap()

p %>% save_x(data = ., name = "PI3K_TF_Heatmap", 1, 16, 10, svg = T)





PI3K_combined_Signaling <- function() {

cabello_aguilar_LR_PMID <- read_csv("~/Downloads/cabello_aguilar_LR_PMID.csv")
cabello <- cabello_aguilar_LR_PMID %>% filter(!is.na(PMIDs))

ligands_de <- Treatment_Response %>% filter(Gene %in% cabello$ligand) %>% pull(Gene) %>% unique()
receptors <- cabello %>% filter(ligand %in% ligands_de) %>% pull(receptor)

receptors_de <- Treatment_Response %>% filter(Gene %in% receptors) %>% pull(Gene) %>% unique()

query <- c(ligands_de, receptors_de)


library(gprofiler2)
library(viridis)

clust = "Combined Signaling"
Results <- gost(query, organism = "hsapiens", evcodes = TRUE)




p <- 
  ggplot(head(filter(Results$result, source == "KEGG"), 5), 
         aes(fct_reorder(term_name, p_value, .desc = T), intersection_size, fill = p_value)) + 
  geom_col() +
  scale_fill_viridis(begin = 0.4, end = 0.7, direction = -1) +
  theme(axis.title.y = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 22, face = "bold")) +
  ylab(label = "Number of Genes in Pathway") +
  ggtitle(paste(unique(head(filter(Results$result, source == "KEGG"), 15)$source), "- Cluster", clust, "Markers")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
  coord_flip()

}

p <- PI3K_combined_Signaling()

p %>% save_x(data = ., name = "PI3K_combined_Signaling", 1, 16, 10, svg = T)
