celltype = "Fibroblast"


meta <- Sobj_integrated@meta.data %>% rownames_to_column("ID") %>% select(ID, CellType, Subcluster, Sample, Type)
query_expr <- LUMexpr %>% filter(CellType == celltype, perc >= 5) %>% pull(TF) %>% unique()
query_tfs  <- cluster.averages %>% left_join(TFs) %>% filter(Cluster == celltype, CF >= 0.5 | TM >= 0.5, TF == "YES") %>% pull(Gene)

query <- TFs_by_Subcluster %>% filter(CellType == celltype) %>% filter(avg_logFC > 0.25 | avg_logFC < -0.25) %>% pull(Gene) %>% unique()
query_module <- query %>% paste0("_95th_1")



Freq_Table <- table(Sobj_integrated$CellType, 
                    Sobj_integrated$Sample) %>% as.data.frame() %>% arrange(-Freq)

enough <- Freq_Table %>% 
  unite(CT_Sample, Var1, Var2, remove = F) %>%
  filter(Freq >= 50) %>% pull(CT_Sample)


data <- readRDS("Output/Reboot/SCENIC/Module_Scores_3.rds")
data <- data %>% rownames_to_column("ID") %>% select(ID, contains("_1")) %>% left_join(meta)

Module_Averages <- 
  data %>%
  group_by(Subcluster, Type) %>% 
  mutate_at(vars(contains("_1")), mean) %>% 
  pivot_longer(contains("_1"), names_to = "Module") %>% 
  select(-ID) %>% 
  unique()
  

Module_Averages  <-   readRDS("Output/Reboot/SCENIC/Module_Averages_95th_Uniqued.rds")  

wilcox <- Module_Averages %>% ungroup() %>% 
  filter(CellType == celltype) %>% 
  pivot_longer(contains("_1"), names_to = "Module") %>%
  filter(Module == query_module) %>% 
  #unite(CT_Sample, CellType, Sample, remove = F) %>%
  #filter(CT_Sample %in% enough) %>% 
  group_by(Subcluster, Module) %>%
  do(w = wilcox.test(value~Type, data = ., paired=FALSE)) %>% 
  summarise(Subcluster, Module, Wilcox = w$p.value) #%>% filter(Wilcox <= 0.05) %>% arrange(Wilcox) %>% add_column(CellType = input_vargenes)



Sig_Module <- wilcox %>% filter(Wilcox <= 0.05) %>% 
  filter(CellType == celltype) %>% 
  separate(Module, into = c("Gene", "B"), sep = "_9", remove = F) %>% 
  filter(Gene %in% query_expr)  %>% arrange(Wilcox) %>% pull(Module)



Module_Avg_Sig <- 
  Module_Averages %>% select(Sample, CellType, Type, Sig_Module) %>% 
  unite(CT_Sample, CellType, Sample, remove = F) %>%
  filter(CellType == celltype, CT_Sample %in% enough) %>%
  #group_by(Sample, CellType) %>% 
  #mutate_at(vars(contains("_1")), mean) %>% 
  #unique() %>%
  pivot_longer(Sig_Module, names_to = "Module") %>% unique()



Modules_High <- 
Module_Avg_Sig %>% group_by(Type, Module) %>% 
  mutate(avg = mean(value)) %>% 
  select(-CT_Sample, -Sample, -value) %>% 
  unique() %>% group_by(Type) %>% 
  filter(avg > 0) %>% pull(Module) %>% unique()
  
  ggplot(aes(x = reorder(Module, -avg), y = avg, fill = Type, color = Type)) + geom_point() +
  
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=14,face="bold"), 
        axis.text.x = element_text(size=10,angle = 90, vjust = 0.5, hjust=1)
  )





Module_Avg_Sig %>% filter(Module %in% Modules_High) %>%
ggplot(aes(x = reorder(Module, -value), y = value, fill = Type)) +
  
  geom_violin(alpha = 0.65, trim = F, scale = "width" , draw_quantiles = T) +
  
  geom_boxplot(aes(color = Type), fill = "white", width = 0.2, alpha = 0.99, 
               position=position_dodge(0.9), outlier.alpha = 0) +
  
  geom_point(aes(color = Type), 
             position = position_jitterdodge(jitter.width = 0.1, 
                                             jitter.height = 0,
                                             dodge.width = 0.9), 
             alpha = 0.99,
             size = 3) +
  
  
  stat_compare_means(paired = F, method = "wilcox", 
                     label.x.npc = 0.25, vjust = -5) +
  
  scale_fill_viridis_d(begin = 0.3, end = 0.8) +
  scale_color_manual(values =  c("gray20", "gray20")) +
  
  ggtitle(label = "Average Module Score per Sample") +
  
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=14,face="bold"), 
        axis.text.x = element_text(size=10,angle = 90, vjust = 0.5, hjust=1)
        ) #+
  #scale_y_sqrt() +
  #facet_wrap(vars(Module), scales = "free_y")

  
  
Module_Avg_Sig %>% filter(Module %in% c("PRRX1_95th_1")) %>%
  ggplot(aes(x = Type, y = value, fill = Type)) +
  
  geom_violin(alpha = 0.65, trim = F, scale = "width" , draw_quantiles = T) +
  
  geom_boxplot(aes(color = Type), fill = "white", width = 0.2, alpha = 0.99, 
               position=position_dodge(0.9), outlier.alpha = 0) +
  
  geom_point(aes(color = Type), 
             position = position_jitterdodge(jitter.width = 0.1, 
                                             jitter.height = 0,
                                             dodge.width = 0.9), 
             alpha = 0.99,
             size = 3) +
  
  
  stat_compare_means(paired = F, method = "wilcox", 
                     label.x.npc = 0.25, vjust = -5) +
  
  scale_fill_viridis_d(begin = 0.3, end = 0.8) +
  scale_color_manual(values =  c("gray20", "gray20")) +
  
  ggtitle(label = "Average Module Score per Sample") +
  
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=14,face="bold"), 
        axis.text.x = element_text(size=10,angle = 90, vjust = 0.5, hjust=1)
  ) +
#scale_y_sqrt() +
facet_wrap(vars(Module), scales = "free_y", ncol = 1)
  
  
  

  
  
  
  
  
  
check <- wilcox %>% filter(Wilcox <= 0.05) %>% 
  filter(CellType == celltype) %>% 
  separate(Module, into = c("Gene", "B"), sep = "_9", remove = F) %>% 
  filter(Module %in% Modules_High) %>% pull(Gene)
  
  
Marker_list    <- vector("list", length = length(check))  

for (i in 1:length(check)) {
  
  clust = check[i]
  query <- percentile %>% filter(TF == clust, percentile_rank > 95) %>% pull(target)
  
  Results <- gost(query, organism = "hsapiens")
  
  
  Marker_list[[i]] <- Results$result %>% dplyr::select(2:6, 10,11) %>% filter(source == "GO:BP") %>%
    mutate(Query = clust) %>% head(10) #%>% head(filter(Results$result, source == "GO:BP")) %>% head(15)
  
}
  


check <- All_Basal %>% distinct(cluster) %>% pull(cluster) %>% as.character()

Marker_list    <- vector("list", length = length(check))  

for (i in 1:length(check)) {
  
  clust = check[i]
  query <- All_Basal %>% filter(cluster == clust, p_val_adj <= 0.05) %>% top_n(75, avg_logFC) %>% pull(gene)
  
  Results <- gost(query, organism = "hsapiens")
  
  
  Marker_list[[i]] <- Results$result %>% dplyr::select(2:6, 10,11) %>% filter(source == "GO:BP") %>%
    mutate(Query = clust) %>% head(10) #%>% head(filter(Results$result, source == "GO:BP")) %>% head(15)
  
}


Bind <- do.call(rbind.data.frame, Marker_list)


  
  #p1 <- 
  ggplot(filter(Bind), 
         aes(fct_reorder(term_name, p_value, .desc = T), intersection_size, fill = p_value)) + 
    geom_col() +
    scale_fill_viridis(begin = 0.25, trans = "reverse", direction = -1) +
    theme(axis.title.y = element_blank(), 
          axis.text = element_text(size = 12), 
          plot.title = element_text(size = 22, face = "bold")) +
    ylab(label = "Number of Genes in Pathway") +
    ggtitle(paste(unique(head(filter(Results$result, source == "GO:BP"), 15)$source), "- Module", clust, "Markers")) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
    coord_flip() + 
    facet_wrap(vars(Query), scales = "free")
  