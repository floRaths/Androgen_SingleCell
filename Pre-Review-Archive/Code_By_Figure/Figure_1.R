data_UMI <- Sobj_integrated@meta.data
data_UMI$CellType <- factor(data_UMI$CellType, levels = levels)

data_NUC <- AtacSeqEmbeddingBasedOnPeakMatrix_20201007
data_NUC$predictedCellTypeGroup <- factor(data_NUC$predictedCellTypeGroup, levels = levels)

obs <- read_csv("~/Downloads/obsRaw.csv")
obs$CellType <- factor(obs$CellType, levels = levels)


cols_type <- c("#A16BAD", "#EFA763")
cols_type_out <- c("#5B2158", "#89582B")


look <-     theme(text = element_text(family = "Lato"),
                  #axis.text.x = element_text(face = "bold", size = 15, angle = -25, hjust = 0, vjust = 1),
                  #axis.title.x = element_blank(),
                  axis.title.y = element_text(face = "bold", size = 18),
                  axis.text.y  = element_text(size = 15), legend.position = "none",
                  
                  axis.title.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_line(color = "grey90"),
                  
                  panel.background = element_rect(fill = "white", color = "grey50", size = 3),
                  
                  strip.text = element_text(size = 13, face = "bold", hjust = 0), 
                  legend.key.size = unit(1,"cm"), legend.justification = c(1,1))





plot_vln_UMI <- function() {
  
  p1 <- data_UMI %>% 
    ggplot(aes(x = Type, y = nCount_RNA)) + 
    
    geom_violin(aes(fill = Type),
                scale = "width", 
                alpha = 0.75, 
                trim = F,
                lwd = 1.25, 
                color = "grey30") + 
    
    geom_boxplot(aes(color = Type), 
                 fill = "grey30", 
                 width = 0.1, alpha = 1, 
                 lwd = 1.25,
                 position=position_dodge(0.7), outlier.alpha = 0.8) +
    
    stat_summary(fun = "median", 
                 geom = "point", size = 1.5,
                 position = position_dodge(0.7),
                 color = "floralwhite") +
    
    
    scale_fill_viridis_d(option = "C", begin = 0.4, end = 0.8) +
    scale_color_manual(values = c("grey30", "grey30")) +
    
    
    
    scale_y_log10(limits = limits_UMI) + 
    
    look
    


  p2 <-  
    data_UMI %>% 
    ggplot(aes(x = CellType, y = nCount_RNA, fill = Type)) + 
    
    geom_boxplot(aes(color = Type), 
                 #fill = "grey30", 
                 width = 0.5, alpha = 0.75, 
                 lwd = 1.25,
                 position=position_dodge(0.7), outlier.alpha = 0.5) +
    
    stat_compare_means(paired = F, method = "wilcox", 
                       label.x.npc = 0.25, vjust = 1) +
    
    stat_summary(fun = "median", 
                 geom = "point", size = 1.5,
                 position = position_dodge(0.7),
                 color = "floralwhite") +
    
    
    scale_fill_viridis_d(option = "C", begin = 0.4, end = 0.8) +
    scale_color_manual(values = c("grey30", "grey30")) +
    
    
    scale_y_log10(limits = limits_UMI) + 
    
    look + theme(axis.title.y = element_blank(), 
                 axis.text.y  = element_blank(), 
                 axis.ticks.y = element_blank())
  
  
  




  
  p3 <- data_NUC %>% 
    ggplot(aes(x = SampleType, y = NucleosomeRatio)) + 
    
    geom_violin(aes(fill = SampleType),
                scale = "width", 
                alpha = 0.75, 
                trim = F,
                lwd = 1.25, 
                color = "grey30") + 
    
    geom_boxplot(aes(color = SampleType), 
                 fill = "grey30", 
                 width = 0.1, alpha = 1, 
                 lwd = 1.25,
                 position=position_dodge(0.7), outlier.alpha = 0.8) +
    
    stat_summary(fun = "median", 
                 geom = "point", size = 1.5,
                 position = position_dodge(0.7),
                 color = "floralwhite") +
    
    
    scale_fill_viridis_d(option = "C", begin = 0.4, end = 0.8) +
    scale_color_manual(values = c("grey30", "grey30")) +
    
    scale_y_log10(limits = limits_NUC) + 
    
    look
  
  
  
  p4 <-  
    data_NUC %>% 
    ggplot(aes(x = predictedCellTypeGroup, y = NucleosomeRatio, fill = SampleType)) + 
    
    geom_boxplot(aes(color = SampleType), 
                 #fill = "grey30", 
                 width = 0.5, alpha = 0.75, 
                 lwd = 1.25,
                 position=position_dodge(0.7), outlier.alpha = 0.5) +
    
    stat_summary(fun = "median", 
                 geom = "point", size = 1.5,
                 position = position_dodge(0.7),
                 color = "floralwhite") +
    
    
    scale_fill_viridis_d(option = "C", begin = 0.4, end = 0.8) +
    scale_color_manual(values = c("grey30", "grey30")) +
    
    
    scale_y_log10(limits = limits_NUC) + 
    
    look + theme(axis.title.y = element_blank(), 
                 axis.text.y  = element_blank(), 
                 axis.ticks.y = element_blank())
  
  
  p5 <- 
    obs %>% 
    select(CellType, Type, ratio = `spliced%`) %>% 
    group_by(Type) %>% 
    mutate(avg_ratio = mean(ratio)) %>% 
    select(-ratio, -CellType) %>% distinct() %>% 
    
    ggplot(aes(x = Type, y = avg_ratio, fill = Type, color = Type)) + 
    
    geom_col(position = position_dodge(), width = 0.8, size = 1) + 
    
    geom_text(aes(label = round(avg_ratio, digits = 1)), 
              position = position_dodge(width = 1)) +
    
    scale_fill_manual(values = cols_type) + 
    scale_color_manual(values = cols_type_out) + 
    
    
    scale_y_continuous(limits = limits_SPL) + 
    
    look
  
  
  
  p6 <- 
    obs %>% 
    select(CellType, Type, ratio = `spliced%`) %>% 
    group_by(CellType, Type) %>% 
    mutate(avg_ratio = mean(ratio)) %>% 
    select(-ratio) %>% distinct() %>% 
    
    ggplot(aes(x = CellType, y = avg_ratio, fill = Type, color = Type)) + 
    
    geom_col(position = position_dodge(), width = 0.75, size = 1) + 
    
    #geom_text(aes(label = round(avg_ratio, digits = 1)), 
    #          position = position_dodge(width = 1)) +
    
    scale_fill_manual(values = cols_type) + 
    scale_color_manual(values = cols_type_out) + 
    
    scale_y_continuous(limits = limits_SPL) + 
    
    look + theme(axis.title.y = element_blank(), 
                 axis.text.y  = element_blank(), 
                 axis.ticks.y = element_blank())
  
  
  
  
  p1 + p2 + p3 + p4 + p5 + p6 + patchwork::plot_layout(ncol = 2, widths = c(1,5))
  
}

limits_UMI = c(400, 30000)
limits_NUC = c(0.3, 1.5)
limits_SPL = c(0, 30.1)

plot_vln_UMI()

plot %>% save_x(data = ., name = paste0("UMI_Nucleosome_Splicing_Combo"), 1, 12, 12, svg = T)









