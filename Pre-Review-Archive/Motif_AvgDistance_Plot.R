mot_list <- NarrowPeaks_DistanceToSummit200_AllCellTypes_Nov16 %>% pull(Motif) %>% unique() %>% as_tibble() %>% separate(value, into = c("TF", "B"), remove = F) %>% select(-B, TF, Motif= "value")

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

plot_distance()






topn <- 
NarrowPeaks_DistanceToSummit200_AllCellTypes_Nov16 %>% 
  filter(CellType == "Basal", GeneRna %in% REAC_Smooth_Muscle_Contraction) %>% 
  select(idxATAC, GeneRna,Motif, Signal, Type) %>% 
  group_by(Motif) %>% 
  summarise(n = n()) %>% 
  arrange(-n) %>% 
  separate(Motif, into = c("Gene", "B"), remove = F) %>% 
  select(-B) %>% 
  left_join(filter(cluster.averages, Cluster == "Basal")) %>% 
  filter(CF > 2, avg_logFC < 0) %>% 
  top_n(10, n)



#diff <- 
  NarrowPeaks_DistanceToSummit200_AllCellTypes_Nov16 %>% 
  filter(CellType == "Basal", GeneRna %in% REAC_Smooth_Muscle_Contraction) %>% 
  select(idxATAC, GeneRna,Motif, Signal, Type) %>%  
  pivot_wider(names_from = Type, values_from = Signal) %>% 
  mutate(FC = log(Trans/Cis)) %>% 
  filter(Motif %in% topn$Motif, Cis > 20, !is.na(FC)) %>% 
  group_by(Motif) %>% 
  summarise(avgFC = mean(FC)) %>% 
  arrange(avgFC) %>% 
    
    ggplot(aes(y = reorder(Motif, - avgFC), x = avgFC)) + 
    geom_col(fill = "lightseagreen") + 
    #geom_text(data = topn, aes(label = n)) +
    
    ggtitle(label = paste("Accessibility loss among 'Smooth Muscle' peaks")) +
    
    ylab(label = "Motif") + 
    xlab(label = "avg loss in ATAC signal") + 
    
    theme(text = element_text(family = "Lato", size = 20), 
        title = element_text(family = "Lato", face = "bold", size = 15), 
        panel.grid = element_line(color = "grey92"), 
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey30", size = 2)
        )
  
  
  
  
  
  
  NarrowPeaks_DistanceToSummit200_AllCellTypes_Nov16 %>% 
    filter(CellType == "Basal", GeneRna %in% REAC_Smooth_Muscle_Contraction) %>% 
    select(idxATAC, GeneRna,Motif, Signal, Type) %>%  
    #pivot_wider(names_from = Type, values_from = Signal) %>% 
    filter(Motif %in% topn$Motif) %>% 
    ggplot(aes(x = Motif, y = Signal, fill = Type)) + 
    geom_violin() + 
    geom_boxplot() + 
    scale_y_log10()