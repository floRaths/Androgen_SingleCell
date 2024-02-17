cols <- c("#9d3396", "#eb933e")

adipocyte_features_2020_10_06_v1 <- read_csv("~/Downloads/adipocyte_features_2020_10_06_v1.csv")

#p2 <- 
adipocyte_features_2020_10_06_v1 %>% 
  ggplot(aes(x = type, y = log10_area)) + 
  
  geom_violin(aes(fill = type),
              scale = "width", 
              alpha = 0.75, 
              trim = F,
              lwd = 1.25, 
              color = "grey30") + 

  geom_boxplot(aes(color = type), 
             fill = "grey30", 
             width = 0.1, alpha = 1, 
             lwd = 1.25,
             position=position_dodge(0.7), outlier.alpha = 0.8) +
  
  scale_fill_viridis_d(option = "C", begin = 0.4, end = 0.8) +
  scale_color_manual(values = c("grey30", "grey30")) +
  
  scale_y_log10() + 
  
  stat_summary(fun = "median", 
               geom = "point", size = 1.5,
               position = position_dodge(0.7),
               color = "floralwhite") +
  
  theme(text = element_text(family = "Lato"),
        axis.text.x = element_text(face = "bold", size = 15, angle = -25, hjust = 0, vjust = 1),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "gray96"),
        strip.text = element_text(size = 13, face = "bold", hjust = 0), 
        legend.key.size = unit(1,"cm"), legend.justification = c(1,1))




#p1 <- 
adipocyte_features_2020_10_06_v1 %>% 
  ggplot(aes(x = case_label, y = log10_area)) + 
  
  geom_violin(aes(fill = type),
              scale = "width", 
              alpha = 0.75, 
              trim = F,
              lwd = 1.25, 
              color = "grey30") + 
  
  geom_boxplot(aes(color = type), 
               fill = "grey30", 
               width = 0.1, alpha = 1, 
               lwd = 1.25,
               position=position_dodge(0.7), outlier.alpha = 0.8) +
  
  scale_fill_viridis_d(option = "C", begin = 0.4, end = 0.8) +
  scale_color_manual(values = c("grey30", "grey30")) +
  
  stat_summary(fun = "median", 
               geom = "point", size = 1.5,
               position = position_dodge(0.7),
               color = "floralwhite") +
  
  theme(text = element_text(family = "Lato"),
        axis.text.x = element_text(face = "bold", size = 15, angle = -25, hjust = 0, vjust = 1),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "gray96"),
        strip.text = element_text(size = 13, face = "bold", hjust = 0), 
        legend.key.size = unit(1,"cm"), legend.justification = c(1,1))





#p3 <- 
adipocyte_features_2020_10_06_v1 %>% group_by(case_label, type) %>% summarise(mean = mean(area)) %>% 
  ggplot(aes(x = type, y = mean, fill = type)) + 
  
  geom_boxplot(aes(color = type), 
               #fill = "grey30", 
               width = 0.3, alpha = 1, 
               lwd = 2,
               position=position_dodge(0.7), outlier.alpha = 0.8) +
  
  
  geom_point(size = 5, alpha = 0.65, color = "grey30",
             position = position_jitterdodge(jitter.width = 0.2,
                                             jitter.height = 0,
                                             dodge.width = 0.9)) +
  
  #scale_fill_viridis_d(option = "C", begin = 0.4, end = 0.8) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = c("grey30", "grey30")) +
  
  stat_compare_means(paired = F, method = "t.test", 
                     label.x.npc = 0.25, vjust = -1) +
  #stat_summary(fun = "median", 
  #             geom = "point", size = 1.5,
  #             position = position_dodge(0.7),
  #             color = "floralwhite") +
  
  theme(text = element_text(family = "Lato"),
        axis.text.x = element_text(face = "bold", size = 15, angle = -25, hjust = 0, vjust = 1),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "gray96"),
        strip.text = element_text(size = 13, face = "bold", hjust = 0), 
        legend.key.size = unit(1,"cm"), legend.justification = c(1,1))






image <- 
  p1 + p2 + p3 + patchwork::plot_layout(widths = c(4,1,1))



save(name = paste0("Adipocyte_Area_Vln"), 1, 16, 6, svg = F)
