Sobj_integrated@meta.data[16:length(colnames(Sobj_integrated@meta.data))] <- NULL


Marker_list    <- vector("list", length = length(levels))

for (i in 1:length(levels)) {
  
celltype = levels[i]

Marker_list[[i]] <- 
  read_tsv(paste0("Output/Reboot/SCENIC/CellType_SCENIC/Output/expr_mat.adjacencies_" , levels[i], ".tsv")) %>% 
  group_by(TF) %>% 
  mutate(percentile_rank = ntile(importance, 100), CellType = levels[i])


}

cluster <- Treatment_Response_Subcluster %>% filter(CellType %in% celltype) %>% pull(Cluster) %>% unique()


Adipo_Unique2 <- "SMARCA4"

up <- Treatment_Response_Subcluster %>% filter(Cluster %in% cluster) %>% left_join(TFs) %>% filter(TF == "YES") %>% group_by(Cluster) %>% top_n(15, avg_logFC)
do <- Treatment_Response_Subcluster %>% filter(Cluster %in% cluster) %>% left_join(TFs) %>% filter(TF == "YES") %>% group_by(Cluster) %>% top_n(-15, avg_logFC)
Adipo_Unique2 <- c(up$Gene, do$Gene, "NR3C1") %>% unique()

for (i in 1:length(Adipo_Unique2)) {
  
  
  genes <- percentile %>% filter(TF == Adipo_Unique2[i], percentile_rank > 95) %>% pull(target)
  
  Sobj_integrated <- AddModuleScore  (Sobj_integrated, 
                                      assay = "RNA", 
                                      #list(filter(Corr_Bind, str_detect(Clust.Name, "^KLK") | str_detect(Clust.Name, "^AR"))$Gene), 
                                      list(genes), 
                                      nbin = 25, name = paste0(Adipo_Unique2[i],"_Module_"))
  
}


wilcox_GRN <- 
Sobj_integrated@meta.data %>% 
  rownames_to_column("ID") %>% 
  select(ID, Type, Subcluster, paste0(Adipo_Unique2, "_Module_1")) %>% 
  pivot_longer(paste0(Adipo_Unique2, "_Module_1"), names_to = "Module") %>% 
  group_by(Subcluster, Module) %>% 
  do(w = wilcox.test(value~Type, data = ., paired = F)) %>% 
  summarise(Subcluster, Module, Wilcox = w$p.value)

avg_vals <- 
Sobj_integrated@meta.data %>% 
  rownames_to_column("ID") %>% 
  select(ID, Type, Subcluster, paste0(Adipo_Unique2, "_Module_1")) %>% 
  pivot_longer(paste0(Adipo_Unique2, "_Module_1"), names_to = "Module") %>% 
  group_by(Subcluster, Module, Type) %>% 
  mutate(avg = mean(value)) %>% 
  select(-ID, -value) %>% 
  unique()



avg_vals %>% 
  filter(Subcluster == "1") %>% 
  pivot_wider(names_from = Type, values_from = avg) %>%
  mutate(FC = log2(TM/CF)) %>% 
  left_join(wilcox_GRN) %>%
  filter(Wilcox <= 0.05, TM > 0.1) %>% 
  arrange(-FC)



mods <- avg_vals %>% 
  filter(Subcluster == "LAMB1+_Matrix-F") %>% 
  pivot_wider(names_from = Type, values_from = avg) %>%
  mutate(FC = log2(TM/CF)) %>% 
  left_join(wilcox_GRN) %>%
  filter(Wilcox <= 0.05, TM > 0.2) %>% 
  arrange(-FC) %>% head(20) %>% pull(Module)



query <- mods

#query <-  c( "NFIC_Module_1", "NR3C1_Module_1")

plot_modules(query, 5)

plot_modules <- function(query, ncol) {
Sobj_integrated@meta.data %>% 
  rownames_to_column("ID") %>% 
  #select(ID, Type, Subcluster, paste0(Adipo_Unique2, "_1")) %>% 
  #pivot_longer(paste0(Adipo_Unique2, "_1"), names_to = "Module") %>% 
  select(ID, Type, Subcluster, query) %>% 
  pivot_longer(query, names_to = "Module") %>% filter(Subcluster != "7") %>% 
  
  #filter(Module == "PPARG_Module") %>% 
  
  ggplot(aes(x = Subcluster, y = value, fill = Type)) +
  
  geom_split_violin(scale = "width", 
                    alpha = 0.75, 
                    trim = F,
                    lwd = 1.25, 
                    color = "grey30") + 
  
  
  geom_boxplot(aes(color = Type), 
               fill = "grey30", 
               width = 0.1, alpha = 1, 
               lwd = 1.25,
               position=position_dodge(0.7), outlier.alpha = 0) +
  
  scale_fill_viridis_d(option = "C", begin = 0.4, end = 0.8) +
  scale_color_manual(values = c("grey30", "grey30")) +
  
  stat_summary(fun = "median", 
               geom = "point", size = 1.5,
               position = position_dodge(0.7),
               color = "floralwhite") +
  
  
  facet_wrap(vars(Module), scales = "free_y", ncol = ncol) +
  
  #scale_y_sqrt() +
  
  theme(text = element_text(family = "Lato"),
        axis.text.x = element_text(face = "bold", size = 15, angle = 25, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "gray96"),
        strip.text = element_text(size = 13, face = "bold", hjust = 0), 
        legend.key.size = unit(1,"cm"), legend.justification = c(1,1))

}





avg_vals %>% 
  filter(Subcluster == "Macrophage") %>% 
  pivot_wider(names_from = Type, values_from = avg) %>%
  mutate(FC = log2(TM/CF)) %>% 
  left_join(wilcox_GRN) %>%
  filter(Wilcox <= 0.05, CF > 0.) %>% 
  arrange(FC) %>% ggplot(aes(x = FC, y = CF)) + geom_point() + geom_text(aes(label = Module))







gmt <- read.gmt("~/Downloads/h.all.v7.2.symbols.gmt")
gmt2 <- read.gmt("~/Downloads/c2.all.v7.2.symbols.gmt")


gmt <- read.gmt("scNuclei-Libraries/Analysis/Reboot/MsigDB/h.all.v7.2.symbols.gmt")
gmt2 <- read.gmt("scNuclei-Libraries/Analysis/Reboot/MsigDB/c2.all.v7.2.symbols.gmt")
gmt <- c(gmt1, gmt2)
 



gmt_names <- names(gmt)


for (i in 1:length(gmt_names)) {
  
  
  genes <- gmt[[gmt_names[i]]]
  
  Sobj_integrated <- AddModuleScore  (Sobj_integrated, 
                                      assay = "RNA", 
                                      #list(filter(Corr_Bind, str_detect(Clust.Name, "^KLK") | str_detect(Clust.Name, "^AR"))$Gene), 
                                      list(genes), 
                                      nbin = 25, name = paste0(gmt_names[i],"_"))
  
}


long_celltype <- Sobj_integrated@meta.data %>% select(Subcluster, Type, contains("_1")) %>% pivot_longer(contains("_1")) %>% group_by(Subcluster, Type, name) %>%  mutate(avg = mean(value)) %>% select(-value) %>% unique()
abs_diff_celltype <- long_celltype %>% pivot_wider(names_from = Type, values_from = avg) %>% group_by(Subcluster, name) %>% mutate(diff = abs(CF - TM))
querym <- abs_diff_celltype %>% filter(Subcluster == "LAMB1+_Matrix-F") %>% arrange(-diff) %>% head(20) %>% pull(name) %>% unique()



Module_Averages <- 
  #Sobj_integrated@meta.data %>% 
  #select(Sample, Type, CellType, Subcluster, contains("_1")) %>% 
  data %>% group_by(Sample, CellType) %>% 
  mutate_at(vars(contains("_1")), mean) %>% 
  unique()


wilcox <- Module_Averages %>% 
  pivot_longer(contains("_1"), names_to = "Module") %>%
  group_by(CellType, Module) %>%
  do(w = wilcox.test(value~Type, data = ., paired=FALSE)) %>% 
  summarise(CellType, Module, Wilcox = w$p.value) #%>% filter(Wilcox <= 0.05) %>% arrange(Wilcox) %>% add_column(CellType = input_vargenes)



Sig_Module <- wilcox %>% filter(Wilcox <= 0.05) %>% 
  filter(Subcluster == celltype) %>% 
  #separate(Module, into = c("Gene", "B"), sep = "_", remove = F) %>% 
  #filter(Gene %in% query_expr)  %>% 
  arrange(Wilcox) %>% pull(Module)



Module_Avg_Sig <- 
  Module_Averages %>% select(Sample, CellType, Type, modules) %>% 
  #unite(CT_Sample, CellType, Sample, remove = F) %>%
  #filter(CellType == celltype, CT_Sample %in% enough) %>%
  #group_by(Sample, Subcluster) %>% 
  #mutate_at(vars(contains("_1")), mean) %>% 
  #unique() %>%
  pivot_longer(contains("_1"), names_to = "Module") %>% unique()


#Modules_High <- 
  Module_Avg_Sig %>% group_by(Type, Module) %>% 
  mutate(avg = mean(value)) %>% 
  select(-Sample, -value) %>% 
  unique() %>% group_by(Type) %>% 
  filter(avg < 0) %>% #pull(Module) %>% unique()

ggplot(aes(x = reorder(Module, -avg), y = avg, fill = Type, color = Type)) + geom_point() +
  
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=14,face="bold"), 
        axis.text.x = element_text(size=10,angle = 90, vjust = 0.5, hjust=1)
  )





Module_Avg_Sig %>% #filter(Subcluster %in% "Lymphatic_EC", Module %in% Modules_High) %>%
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
  ) +
#scale_y_sqrt() +
facet_wrap(vars(Subcluster), scales = "free_y")





