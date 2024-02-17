
library(gprofiler2)
library(viridis)

clust = "0"
#query <- Treatment_Response_Subclus2er %>% filter(Cluster == clust, avg_logFC > 0.25, Gene %notin% common) %>% top_n(100, avg_logFC) %>% pull(Gene)
#query <- against_all %>% filter(p_val_adj <= 0.05, pct.1 > 0.3) %>% top_n(150, avg_logFC) %>% pull(gene)
#query <- mrks %>% filter(p_val_adj <= 0.05, cluster == clust, pct.1 > 0.5) %>% top_n(100, avg_logFC) %>% pull(gene)
query <- percentile %>% filter(TF == clust, percentile_rank > 90) %>% pull(target)
#query <- cutree(image$tree_col, k = 4) %>% as.data.frame() %>% rename(Clust = ".") %>% filter(Clust == 1) %>% rownames_to_column("Gene") %>% pull(Gene)
#query <- `Markers_LUM_HR-neg` %>% filter(cluster == clust, p_val_adj <= 0.05) %>% top_n(100, avg_logFC) %>% pull(gene)
#query <- Treatment_Response %>% filter(Cluster == "Fibroblast", Gene %notin% common, avg_logFC < -0.25, pct.2 > 0.3) %>% pull(Gene)
Results <- gost(query, organism = "hsapiens", evcodes = T)

 
Results$result %>% dplyr::select(2:6, 10,11) #%>% head(filter(Results$result, source == "GO:BP")) %>% head(15)



#p4 <- 
ggplot(head(filter(Results$result, source == "GO:BP"), 15), 
       aes(fct_reorder(term_name, p_value, .desc = T), intersection_size, fill = p_value)) + 
  geom_col() +
  scale_fill_viridis(begin = 0.25, trans = "reverse", direction = -1) +
  theme(axis.title.y = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 22, face = "bold")) +
  ylab(label = "Number of Genes in Pathway") +
  ggtitle(paste(unique(head(filter(Results$result, source == "GO:BP"), 15)$source), "-", clust, "Module")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
  coord_flip()

#p2 <- 
ggplot(head(filter(Results$result, source == "REAC"), 15), 
       aes(fct_reorder(term_name, p_value, .desc = T), intersection_size, fill = p_value)) + 
  geom_col() +
  scale_fill_viridis(begin = 0.25, trans = "reverse", direction = -1) +
  theme(axis.title.y = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 22, face = "bold")) +
  ylab(label = "Number of Genes in Pathway") +
  ggtitle(paste(unique(head(filter(Results$result, source == "REAC"), 15)$source), "- Cluster", clust, "Markers")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
  coord_flip() + scale_y_sqrt()


#p %>% save_x(data = ., name = "KEGG_Signaling_Enrichment", 1, 10, 5, svg = T)


#p <- 
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

  
  ggplot(head(filter(Results$result, source == "WP"), 8), 
         aes(fct_reorder(term_name, p_value, .desc = T), intersection_size, fill = p_value)) + 
  geom_col() +
    scale_fill_viridis(begin = 0.4, end = 0.7, direction = -1) +
  #scale_fill_viridis(begin = 0.25, trans = "reverse", direction = -1) +
  theme(axis.title.y = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 22, face = "bold")) +
  ylab(label = "Number of Genes in Pathway") +
  ggtitle(paste(unique(head(filter(Results$result, source == "WP"), 15)$source), "- Cluster", clust, "Markers")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
  coord_flip()


p %>% save_x(data = ., name = paste0("GPAM_Enrichment"), 1, 10, 5, svg = T)
