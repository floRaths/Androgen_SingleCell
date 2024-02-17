query <- c("PTEN", "PIK3R1", "AKT3")


p <- 
Sobj_integrated@assays$RNA@data[query, ] %>% as.matrix() %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("ID") %>% 
  left_join(select(rownames_to_column(Sobj_integrated@meta.data, "ID"), ID, CellType, Type)) %>% 
  #filter( == "Fibroblast") %>%
  #filter(Subcluster %in% levels(Sobj_integrated)[1:4]) %>% 
  pivot_longer(query, names_to = "Gene") %>% filter(Gene != "AR") %>% 
  #mutate(value = value + 1) %>% 
  
  ggplot(aes(x = CellType, y = value, fill = Type)) + 
  geom_violin(size = 1, color = "grey20", scale = "width") + 
  
  scale_fill_manual(values = c("#9d3396", "#eb933e")) + 
  
  #scale_y_log10() + 
  #scale_y_sqrt() + 
  facet_wrap(facets = c("Gene"), scales = "free_y", ncol = 4) + 
  
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




#P2 <- 
adipocyte_features_2020_10_06_v1 %>% #filter(case_label %notin% c("CF-787-821", "TM-7714")) %>% 
  group_by(case_label, type) %>% 
  summarise(mean = mean(area)) %>% 
  ggplot(aes(x = type, y = mean, fill = type)) + 
  
  geom_boxplot(#aes(color = type), 
               color = "grey20",
               width = 0.5, alpha = 1, 
               lwd = 1,
               position=position_dodge(0.7), outlier.alpha = 0.8) +
  
  
  geom_point(color = "grey60", size = 4, alpha = 0.65, 
             position = position_jitterdodge(jitter.width = 0.2,
                                             jitter.height = 0,
                                             dodge.width = 0.9)) +
  
  
  scale_fill_manual(values = c("#9d3396", "#eb933e")) +
  scale_color_manual(values = c(cols, "grey")) +
  
  stat_compare_means(paired = F, method = "t.test", 
                     label.x.npc = 0.25) +
  
    ggtitle("avg. adipocyte area") +
  
  #coord_flip() + 
  
  theme(text = element_text(family = "Lato", size = 20), 
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(color = "grey30", size = 2, fill = NA),
        panel.grid = element_line(color = "grey90"))




p <- P1 + P2 + patchwork::plot_layout(ncol = 1, heights = c(2.5,1))


p %>% save_x(data = ., name = "Lipolysis_Markers_And_Area", 1.5, 5, 5, svg = T)
