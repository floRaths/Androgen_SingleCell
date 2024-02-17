EnrichmentOf <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/ATAC_Integration/Motif_Enrichment_RDS_Files/EnrichmentOfExpressedTfsOrMotifsInDifferentCellTypes_Oct25.RDS")
Motif_Summary <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/ATAC_Integration/Motif_Summary_Complete_Unfiltered.rds")


Enrichment <- 
  EnrichmentOf %>% 
  mutate(name = "Motif.cisBP") %>% 
  unite(Motif.Long, name, TF, sep = ".", remove = F) %>% 
  separate(Comparison, into = c("CellType", "B"), sep = " ", extra = "merge") %>% 
  separate(B, into = c("Type", "D"), sep = " ", extra = "merge") %>% 
  dplyr::select(-D, -name) %>% 
  dplyr::rename(Motif = "TF") %>% 
  left_join(Motif_Summary, by = c("Motif.Long" = "Motif", "CellType")) %>% 
  as_tibble()



plot_enrichment <- function() {
  
  motifs <- 
    Enrichment %>% 
    mutate(Enrichment_Div = ifelse(Type == "Trans", Enrichment, Enrichment * -1), 
           Avg_Div = ifelse(Type == "Trans", avg_TM, avg_CF * -1), FracPeak = (nPeaksWithMotif/nPeaks) * 100) %>% 
    filter(CellType == celltype, 
           abs(Avg_Div) > 0.5,
           FracPeak >= 10) %>% 
    group_by(Type) %>% top_n(-8, rank) %>% pull(Motif) %>% unique()
  
  
  p <-   
    Enrichment %>% filter(CellType == celltype, Motif %in% motifs) %>% 
    mutate(Enrichment_Div = ifelse(Type == "Trans", Enrichment, Enrichment * -1), 
           Avg_Div = ifelse(Type == "Trans", avg_TM, avg_CF * -1)) %>% 
    
    ggplot(aes(y = reorder(Motif, Enrichment_Div), x = Enrichment_Div, fill = Type)) + 
    
    geom_col(width = 0.85, alpha = 1) + 
    geom_text(aes(label = nPeaks)) +
    #geom_col(aes(x = Avg_Div), width = 0.85, alpha = 0.1) + 
    
    #scale_fill_gradient2(low = "#9d3396", mid = "floralwhite", high = "#eb933e") +
    scale_fill_manual(values = c("#9d3396", "#eb933e")) +
    
    ggtitle(label = paste0(celltype, " Motif Enrichment")) + 
    scale_y_discrete(position = "right")+
    #ylim(-xlim, xlim) +
    
    geom_vline(xintercept = 0, lwd = 3, color = "white") +
    
    theme(text = element_text(family = "Lato", size = 25),
          legend.position = "none", 
          
          title = element_text(size = 20, colour = "grey20"),
          
          axis.text.x = element_text(),
          axis.text.y = element_text(),
          
          panel.grid.major.x  = element_blank(),
          panel.grid.minor.x  = element_blank(),
          panel.grid.minor.y  = element_blank(),
          panel.grid.major.y  = element_line(color = "grey70", 
                                             lineend = "round",
                                             linetype = "dotted", 
                                             size = 1),
          
          panel.background = element_rect(fill = "white", 
                                          color = "grey70", 
                                          size = 2),
          
          axis.title.y = element_blank(),
          axis.title.x = element_blank())
  
  
  
  
  
  
  q <- 
    Enrichment %>% filter(CellType == celltype, Motif %in% motifs) %>% 
    mutate(FracPeak = (nPeaksWithMotif/nPeaks) * 100, Enrichment_Div = ifelse(Type == "Trans", Enrichment, Enrichment * -1), 
           Avg_Div = ifelse(Type == "Trans", avg_TM, avg_CF * -1)) %>% 
    
    ggplot(aes(y = reorder(Motif, Enrichment_Div), x = FracPeak, fill = Type)) + 
    geom_col(width = 0.85, alpha = 1, position = position_dodge()) + 
    scale_fill_manual(values = c("#9d3396", "#eb933e")) +
    scale_x_reverse() + theme(legend.position = "none", 
                              axis.text.y = element_blank(), 
                              axis.ticks.y = element_blank(), 
                              axis.title.y = element_blank(), 
                              axis.title.x = element_blank(), 
                              panel.background = element_rect(fill = "white", color = "grey70", size = 2), 
                              panel.grid.major = element_line(color = "grey90"),
                              text = element_text(family = "Lato", size = 25))
  
  
  q + p + patchwork::plot_layout(widths = c(1, 4))
  
}

celltype = "LUM_HR-pos"

p <- plot_enrichment()


p %>% save_x(data = ., name = paste0("TF_Enrichment_", celltype), 1, 7, 7, svg = T)    



