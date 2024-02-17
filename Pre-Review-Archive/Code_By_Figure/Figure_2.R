Motif.cisBP_chromvar_zscores <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/ATAC_Integration/INSR_Story/ChromVar/Motif.cisBP_chromvar_zscores.RDS")
correlation_tf_expression_motif_activity <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/ATAC_Integration/INSR_Story/ChromVar/correlation_tf_expression_motif_activity.rds")

z <- Motif.cisBP_chromvar_zscores@assays@data@listData$z %>% as.matrix() %>% t()
meta <- Motif.cisBP_chromvar_zscores@colData

meta_brief <- meta %>% as.data.frame() %>%  
  rownames_to_column("ID") %>% 
  select(ID, Sample, CellType = "predictedCellTypeGroup", Type = "SampleType", predictedGroup)

z_df <- z %>% as.data.frame() %>% rownames_to_column("ID")




df <- meta_brief %>% left_join(select(z_df, ID, contains("_"))) %>% 
  pivot_longer(contains("_"), names_to = "Motif", values_to = "Score")


df <- df %>% separate(Motif, into = c("TF", "num"), remove = F)


#p <- 
  df %>% filter(CellType == "LUM_HR-pos", TF %in% c("AR")) %>% 
  
  ggplot(aes(x = CellType, y = Score, fill = Type)) + 
  geom_violin(scale = "width") +
  geom_boxplot(aes(color = Type), size = 1,
               fill = "grey30", 
               width = 0.2, 
               position = position_dodge(width = 0.9), outlier.alpha = 0) + 
  
  stat_summary(fun = "median", 
               geom = "point", size = 3,
               position = position_dodge(0.9),
               color = "floralwhite") +
  
  scale_color_manual(values = c("grey30", "grey30")) + 
  #scale_fill_manual(values = cols_type) + 
  stat_compare_means(paired = F, method = "wilcox", 
                       label.x.npc = 0.25, vjust = 1) +
    
  facet_wrap(vars(TF), scales = "free_y", ncol = 3) + 
  
  theme(text = element_text(family = "Lato", size = 40), 
        axis.title = element_blank(), 
        panel.grid = element_line(colour = "grey90"),
        panel.background = element_blank(), 
        legend.position = "none",
        panel.border = element_rect(size = 2, color = "grey30", fill = NA))


p %>% save_x(data = ., name = "ChromVar_AR_LUM_HR-pos", 1, 7, 10, svg = T)  






dat <- 
Sobj_integrated@assays$RNA@counts[c("AR", "ESR1", "PGR"), ] %>% 
  as.matrix() %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("ID") %>% 
  left_join(select(rownames_to_column(Sobj_integrated@meta.data, "ID"), ID, Type, Sample)) %>% 
  column_to_rownames("ID") %>% mutate(AR = ifelse(AR == 0, 0, 1), ESR1 = ifelse(ESR1 == 0, 0, 1), PGR = ifelse(PGR == 0, 0, 1))





#p2 <- 
dat %>% 
  group_by(Type, Sample, AR, ESR1, PGR) %>% 
  summarise(n = n()) %>% 
  group_by(Sample) %>% 
  mutate(total = sum(n), perc = (n/total)*100) %>% 
  
  mutate(A_ESR1      = ifelse((AR == 0 & ESR1 == 1 & PGR == 0), "Y", "N"),
         B_AR        = ifelse((AR == 1 & ESR1 == 0 & PGR == 0), "Y", "N"),
         C_PGR       = ifelse((AR == 0 & ESR1 == 0 & PGR == 1), "Y", "N"),
         
         D_ESR1_AR   = ifelse((AR == 1 & ESR1 == 1 & PGR == 0), "Y", "N"),
         E_ESR1_PGR  = ifelse((AR == 0 & ESR1 == 1 & PGR == 1), "Y", "N"),
         F_AR_PGR    = ifelse((AR == 1 & ESR1 == 0 & PGR == 1), "Y", "N"),
         
         G_All       = ifelse((AR == 1 & ESR1 == 1 & PGR == 1), "Y", "N"),
         H_None      = ifelse((AR == 0 & ESR1 == 0 & PGR == 0), "Y", "N")
         )  %>% select(-H_None) %>% 
  
  ungroup() %>% 
  pivot_longer(9:15, names_to = "Intersect") %>% filter(value == "Y") %>% 
  
  ggplot(aes(x = Intersect, y = perc, fill = Type)) + 
  geom_boxplot(outlier.alpha = 0, size = 1.5, color = "grey30", width = 0.65) + 
  
  geom_point(aes(color = factor(Sample, levels = levels_samp)), 
             position = position_dodge(width = 0.75), 
             size = 4, 
             alpha = 1) + 
  
  scale_fill_manual(values = c("#9d3396", "#eb933e")) + 
  scale_color_manual(values = cols_samp) + 
  
  stat_compare_means(paired = F, method = "t.test", 
                     label.x.npc = 0.25, vjust = 1) + 
  
  theme(text = element_text(family = "Lato", size = 20), 
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "grey30", size = 2, fill = NA),
        panel.grid = element_blank(), 
        legend.position = "none"
  )


p <- as.ggplot(p1) + p2 + patchwork::plot_layout(ncol = 1) 

p %>% save_x(data = ., name = "LUM_HR-pos_UPSET_HR", 1, 16, 9, svg = T)  










plot_distance <- function() {
  mot_dist <- NarrowPeaks_DistanceToSummit200_AllCellTypes_Nov16 %>% 
    filter(CellType == celltype, Motif == motif) %>% 
    select(idxATAC, Mot = "Motif", DistMot = "DistanceToSummit")
  
  df <- 
    NarrowPeaks_DistanceToSummit200_AllCellTypes_Nov16 %>% 
    filter(CellType == celltype) %>% 
    left_join(mot_dist) %>% 
    select(idxATAC, GeneRna, Motif, Mot, Signal, Type, DistanceToSummit, DistMot) %>% 
    filter(!is.na(Mot)) %>% 
    distinct() %>% 
    mutate(diff = abs(DistMot - DistanceToSummit)) %>% 
    group_by(Motif) %>% 
    summarise(n = n(), mean = mean(diff)) %>% 
    arrange(-n)
  
  df <- df %>% mutate(avg = median(mean))
  avg <- df %>% pull(avg) %>% unique()
  
  filter(df, Motif != motif) %>% 
    ggplot(aes(x = n, y = mean)) + 
    
    geom_hline(aes(yintercept = avg), color = "tomato", alpha = 0.35) +
    
    geom_point(color = "lightseagreen", size = 2) + 
    geom_text_repel(data = filter(df, Motif != motif, 
                                  n > 2*(mean(df$n)) | mean < 90), 
                    aes(label = Motif), color = "grey20") +
    
    ggtitle(label = paste(motif, "   avg. distance from co-occuring motifs in", celltype)) +
    
    ylab(label = "avg. motif distance") + 
    xlab(label = "# co-occurrences") + 
    
    annotate(geom = "text", x = max(filter(df, Motif != motif)$n), y = (avg)-3, label = "median", color = "tomato", alpha = 0.35) +
    
    theme(text = element_text(family = "Lato", size = 20), 
          title = element_text(family = "Lato", face = "bold", size = 15), 
          panel.grid = element_blank(), 
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "grey30", size = 2)
    )
  
}

celltype = "LUM_HR-pos"
motif = filter(mot_list, TF == "AR") %>% pull(Motif)

p <- plot_distance()

p %>% save_x(data = ., name = "LUM_HR-pos_AR-Distance_Plot", 1, 16, 9, svg = T)  




fil_data <- tissue_enrichment_data_AllSamples_gender_split %>% separate(tissue, into = c("Tissue", "Gender"), sep = "_") %>% filter(pval <= 0.05)

p <- 
tissue_enrichment_data_AllSamples_gender_split %>% 
  separate(tissue, into = c("Tissue", "Gender"), sep = "_") %>%
  ggplot(aes(x = fc_vs_median, y = -log10(pval), color = Gender)) + 
  geom_vline(xintercept = 0, color = "grey70") + 
  geom_hline(yintercept = -log10(0.05), color = "grey70", linetype = "dashed") + 
  geom_point(size = 5) + 
  geom_text_repel(data = fil_data, aes(label = Tissue), color = "grey30", size = 5) +
  scale_color_manual(values = c("#9d3396", "#eb933e")) +
  
  theme(text = element_text(family = "Lato", size = 20), 
                       title = element_text(family = "Lato", face = "bold", size = 15), 
                       panel.grid = element_blank(), 
                       panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "grey30", size = 2)
                       )

p %>% save_x(data = ., name = "LUM_HR-pos_AR-tissue_association", 1, 12, 9, svg = T)  








Sobj_integrated <- AddMetaData(Sobj_integrated, readRDS("Output/LUM-Pos_AR_Status_Metadata.rds"))

JUNFOS_Split <- function(){

query <- c("JUN", "FOSB")


p <- 
  Sobj_integrated@assays$RNA@data[query, ] %>% as.matrix() %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("ID") %>% 
  left_join(select(rownames_to_column(Sobj_integrated@meta.data, "ID"), ID, CellType, Type, Subcluster, ARtat)) %>% unite(ARstat, Type, ARtat, remove = F) %>% 
  #filter(CellType == "Fibroblast") %>%
  #filter(Subcluster %in% levels(Sobj_integrated)[1:3]) %>% 
  pivot_longer(query, names_to = "Gene") %>% filter(Gene != "AR") %>% 
  #mutate(value = value + 1) %>% 
  
  ggplot(aes(x = Type, y = value, fill = ARtat)) + 
  geom_violin(size = 1, color = "grey20", scale = "width") + 
  
  scale_fill_manual(values = c("#9d3396", "#eb933e")) + 
  
  #scale_y_log10() + 
  #scale_y_sqrt() + 
  facet_wrap(facets = c("Gene", "ARtat"), scales = "free_y", ncol = 4) + 
  
  ylab("expression") + 
  
  theme(text = element_text(family = "Lato", size = 20), 
        legend.position = "top", 
        #axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(color = "grey30", size = 2, fill = NA),
        panel.grid = element_blank())

}

p %>% save_x(data = ., name = paste0("LUM-HRpos_JUN_FOS_Split"), 1, 10, 5, svg = T)





linearModel <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/ATAC_Integration/LUM_HR-pos_Story/linearModelNuclearReceptorsExpressionMotifActivity_all.RDS")

mat <- linearModel %>% 
  select(Gene, Motif, effectSize) %>% #filter(Gene != "ESR1", Motif != "ESR1_661") %>%  
  pivot_wider(names_from = Gene, values_from = "effectSize") %>% 
  column_to_rownames("Motif")


my.breaks <- c(seq((min(mat)), -0.00001, length.out = 50), seq(0.00001, (max(mat)), length.out = 50))

p <- 
pheatmap(mat[c(-7,-8) , c(-7, -8)], breaks = my.breaks,
         color = colorRampPalette(brewer.pal(11, "BrBG"))(100))

p %>% save_x(data = ., name = paste0("LUM-HRpos_NucReceptors_linearModel"), 1, 10, 5, svg = T)




permutation <- function(){

up <- read_csv("~/Downloads/tissue_enrichment_data_MaleSamples_FC0.20.csv") %>% mutate(Set = "up")
do <- read_csv("~/Downloads/tissue_enrichment_data_MaleSamples_FC-0.20.csv") %>% mutate(Set = "down")
AR_expression_MaleSamples <- read_csv("~/Downloads/AR_expression_MaleSamples.csv") %>% select(tissue = "Tissue", AR_median = "value")


data <- 
bind_rows(up, do) %>% 
  left_join(AR_expression_MaleSamples, by = c("tissue")) %>% 
  mutate(fc_vs_median = replace_na(fc_vs_median, 0)) %>% 
  mutate(label = ifelse(pval <= 0.05, tissue, "n.s."))
  

#p <-   
data %>% 
ggplot(aes(x = fc_vs_median, y = AR_median)) + 
  geom_vline(xintercept = 0, color = "grey90") +
  geom_text_repel(aes(label = tissue), size = 5) + 
  geom_point(data = filter(data, pval <= 0.05), aes(color = tissue), size = 5) + 
  geom_point(data = filter(data, label == "n.s."), color = "grey80", size = 5) + 
  geom_smooth(method = "lm", se = F) + 
  scale_color_d3() +
  facet_wrap(facets = "Set", ncol = 1) +   
  theme(text = element_text(family = "Lato", size = 20), 
                                     legend.position = "top", 
                                     #axis.text.x = element_blank(), 
                                     #axis.ticks.x = element_blank(), 
                                     axis.title.x = element_blank(),
                                     #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                                     panel.background = element_blank(),
                                     panel.border = element_rect(color = "grey30", size = 2, fill = NA),
                                     panel.grid = element_blank())



}

p <- permutation()

p %>% save_x(data = ., name = paste0("LUM-HRpos_permutation"), 1, 10, 10, svg = T)





itgam %>% filter(Expression == "Present") %>%
  separate(Replicate, into = c("Group", "Rest"), sep = "_", remove = F) %>%
  ggplot(aes(x = Replicate, y = Freq, fill = Group)) + 
  geom_col() + 
  ylim(c(0,1)) +
  scale_fill_brewer(palette = "Dark2")
  
  
