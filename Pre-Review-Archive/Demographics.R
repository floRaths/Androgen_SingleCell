Cohort_demographics <- read_csv("utilities/Cohort_demographics.csv")

metadata <- Sobj_integrated@meta.data %>% 
  rownames_to_column("ID") %>% 
  left_join(Cohort_demographics, by = "Sample") %>% 
  column_to_rownames("ID")

Sobj_integrated <- AddMetaData(Sobj_integrated, metadata = metadata)



table(Sobj_integrated$Mens.stat) %>% as.data.frame() %>% arrange(-Freq)
Idents(Sobj_integrated) <- "Mens.stat"
pre_Subsample  <- subset(Sobj_integrated, idents = c("na", "pre"),  downsample = 24369)
post_Subsample <- subset(Sobj_integrated, idents = c("na", "post"), downsample = 24369)
Idents(Sobj_integrated) <- "Type"
comb_subset <- subset(Sobj_integrated, idents = c("TM", "CF"), downsample = 24369)

DimPlot(Sobj_integrated, group.by =  "Mens.stat", pt.size = 0.5, split.by = "Type")


Subsample@meta.data %>% 
  select(Sample, Age, Mens.stat, Subcluster, Type) %>% 
  as_tibble() %>% 
  group_by(Subcluster) %>% 
  mutate(total = n()) %>% 
  group_by(Subcluster, Mens.stat) %>% 
  mutate(count = n(), perc = (count/total)*100) %>% 
  select(Mens.stat, Subcluster, perc) %>% 
  distinct() %>% 
  
  ggplot(aes(x = Subcluster, y = perc, fill = Mens.stat)) + 
  geom_col(position = position_dodge()) + 
  scale_fill_d3()



mens.response <- 
Treatment_Response %>% 
  mutate(avg_log2FC = avg_logFC, Mens.stat = "both") %>% 
  select(-avg_logFC) %>% 
  bind_rows(Bind)

mens.response <-  readRDS("Output/TreatResponse_Mens.Status_Downsampled24K_each.rds")
mens.response <-  readRDS("Output/Mens.Status_Separated_Downsampled_Response.rds")
#mens.response <-  readRDS("Output/Mens.Status_Separated_EqualCellType_Downsampled_Response.rds")



comparison_heatmap <- function(x) {
  
  Clusters <- c("LUM_HR-pos", "LUM_HR-neg", "Basal", "Fibroblast", "Adipocyte", "Blood_EC", "Lymphatic_EC", "Vasc_Accessory", "Myeloid", "Lymphoid")
  #Clusters <- Treatment_Response %>% pull(Cluster) %>% unique()

  Marker_list    <- vector("list", length = length(Clusters))


for (i in seq_along(Marker_list)) {
  
  topgenes <- 
    Treatment_Response %>% 
    filter(Cluster == Clusters[i], abs(avg_logFC) > x, p_val_adj <= 0.05) %>% 
    arrange(-avg_logFC) %>% 
    pull(Gene)
  
  
  mat <- 
    mens.response %>% filter(Cluster == Clusters[i], Gene %in% topgenes, p_val_adj <= 0.05) %>%
    select(Gene, avg_log2FC, Mens.stat) %>% 
    pivot_wider(names_from = Mens.stat, values_from = avg_log2FC) %>% 
    column_to_rownames("Gene")
  
  mat <- mat[, c("both", "pre", "post")]
  
  mat[is.na(mat)] <- 0
  my.breaks <- c(seq((min(mat)), -0.00001, length.out = 50), seq(0.00001, (max(mat)), length.out = 50))
  my.colors <- c(colorRampPalette(colors = c("#002c2d", "#00494b", "#076769", "#208288", "#3ea8a6", "#76c7be", "#bbe4d1", "#F2F2F2"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("#F2F2F2", "#fed693", "#f5ad52", "#dd8629", "#bd651a", "#9c4511", "#79260b", "#580000"))(length(my.breaks)/2))
  
  Marker_list[[i]] <- 
    pheatmap(mat, border_color = NA, treeheight_row = 0, treeheight_col = 0,
             show_rownames = F, cluster_cols = F,
             color = my.colors,
             #annotation_row = Ann,
             breaks = my.breaks, main = Clusters[i]) %>% 
    as.ggplot()
}

extract_plots <- function(){
p1 <- Marker_list[[1]]
p2 <- Marker_list[[2]]
p3 <- Marker_list[[3]]
p4 <- Marker_list[[4]]
p5 <- Marker_list[[5]]
p6 <- Marker_list[[6]]
p7 <- Marker_list[[7]]
p8 <- Marker_list[[8]]
p9 <- Marker_list[[9]]
p10 <- Marker_list[[10]]

p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + patchwork::plot_layout(ncol = 5)

}
extract_plots()

}
comparison_heatmap_downsampled <- function(x) {
  
  #Clusters <- mens.response %>% pull(Cluster) %>% unique()
  Clusters <- c("LUM_HR-pos", "LUM_HR-neg", "Basal", "Fibroblast", "Adipocyte", "Blood_EC", "Lymphatic_EC", "Vasc_Accessory", "Myeloid", "Lymphoid")
  Marker_list    <- vector("list", length = length(Clusters))
  
  
  for (i in seq_along(Marker_list)) {
    
    topgenes <- 
      mens.response %>% 
      filter(Mens.stat == "both", Cluster == Clusters[i], abs(avg_log2FC) > x, p_val_adj <= 0.05) %>% 
      arrange(-avg_log2FC) %>% 
      pull(Gene)
    
    
    mat <- 
      mens.response %>% filter(Cluster == Clusters[i], Gene %in% topgenes, p_val_adj <= 0.05) %>%
      select(Gene, avg_log2FC, Mens.stat) %>% 
      pivot_wider(names_from = Mens.stat, values_from = avg_log2FC) %>% 
      column_to_rownames("Gene")
    
    mat <- mat[, c("both", "pre", "post")]
    
    mat[is.na(mat)] <- 0
    my.breaks <- c(seq((min(mat)), -0.00001, length.out = 50), seq(0.00001, (max(mat)), length.out = 50))
    my.colors <- c(colorRampPalette(colors = c("#002c2d", "#00494b", "#076769", "#208288", "#3ea8a6", "#76c7be", "#bbe4d1", "#F2F2F2"))(length(my.breaks)/2), 
                   colorRampPalette(colors = c("#F2F2F2", "#fed693", "#f5ad52", "#dd8629", "#bd651a", "#9c4511", "#79260b", "#580000"))(length(my.breaks)/2))
    
    Marker_list[[i]] <- 
      pheatmap(mat, border_color = NA, treeheight_row = 0, treeheight_col = 0,
               show_rownames = F, cluster_cols = F,
               color = my.colors,
               #annotation_row = Ann,
               breaks = my.breaks, main = Clusters[i]) %>% 
      as.ggplot()
  }
  
  extract_plots <- function(){
    p1 <- Marker_list[[1]]
    p2 <- Marker_list[[2]]
    p3 <- Marker_list[[3]]
    p4 <- Marker_list[[4]]
    p5 <- Marker_list[[5]]
    p6 <- Marker_list[[6]]
    p7 <- Marker_list[[7]]
    p8 <- Marker_list[[8]]
    p9 <- Marker_list[[9]]
    p10 <- Marker_list[[10]]
    
    p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + patchwork::plot_layout(ncol = 5)
    
  }
  extract_plots()
  
}

comparison_heatmap(0.25)
p <- comparison_heatmap_downsampled(0.5)


save_x(p, name = "Mens.Status_DE_Impact", scale = 1, w = 12, h = 7, svg = T)





$p <- 
comb_subset@meta.data %>% 
  select(Sample, Mens.stat, CellType, Type) %>% 
  group_by(Type, CellType, Mens.stat) %>% 
  summarise(n = n()) %>% 
  arrange(CellType) %>% group_by(CellType) %>% 
  mutate(total = sum(n)) %>% filter(Mens.stat != "na") %>% 
  
  ggplot(aes(x = factor(CellType, levels = levels), y = n, fill = factor(Mens.stat, levels = c("na", "pre", "post")))) + 
  geom_col(position = position_fill()) + 
  scale_fill_d3(name = "Mens.status") + 
  theme(text = element_text(size = 15))

save_x(p, name = "Mens.Status_ComparisonProportion", scale = 1, w = 12, h = 7, svg = T)


p <- Sobj_integrated@meta.data %>% 
  select(Sample, Age, Mens.stat) %>% 
  distinct() %>% 
  as_tibble() %>% 
  group_by(Mens.stat) %>% 
  mutate(median = median(Age)) %>% 
  ggplot(aes(x = Mens.stat, y = Age, fill = factor(Mens.stat, levels = c("na", "pre", "post")))) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_point(size = 4, alpha = 0.5) + 
  scale_fill_d3(name = "Mens.status") + 
  theme(text = element_text(size = 15))


save_x(p, name = "Mens.Status_Age_Boxplot", scale = 1, w = 5, h = 7, svg = T)








Sobj@meta.data %>% 
  rownames_to_column("SAMPID") %>% 
  as_tibble() %>% 
  select(Gender, SAMPID, Fibroblast_1) %>% filter(between(Fibroblast_1, 0.025, 0.2)) %>% 
  ggplot(aes(fill = Gender, x = reorder(SAMPID, -Fibroblast_1 ), y = Fibroblast_1)) + 
  geom_col() + facet_wrap(facets = "Gender")




Sobj@meta.data %>% 
  rownames_to_column("SAMPID") %>% 
  as_tibble() %>% 
  select(SAMPID, Gender, Fibroblast_1) %>% 
  filter(between(Fibroblast_1, 0.1, 0.5)) %>% 
  ggplot(aes(fill = Gender, x = Gender, y = Fibroblast_1)) + 
  geom_boxplot()

















Sobj_integrated@meta.data %>% 
  select(Sample, Subcluster, Type, Mens.stat, Age) %>% 
  #group_by(Subcluster) %>% 
  mutate(total = n()) %>% 
  group_by(Sample, Subcluster) %>% 
  mutate(n = n(), perc = (n/total)*100) %>% 
  distinct() %>% 
  
  ggplot(aes(x = Subcluster, y = perc, fill = Type)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = position_jitterdodge())





table(Sobj_integrated$Sample, Sobj_integrated$Subcluster) %>% 
  as.data.frame() %>% 
  separate(Var1, into = c("Type", "B"), extra = "merge", remove = F) %>% 
  as_tibble() %>% rename(Sample = "Var1", Subcluster = "Var2")
  




meta <- Sobj_integrated@meta.data %>% 
  select(Sample, Mens.stat, Age) %>% distinct()
  
  
prop.table(table(Sobj_integrated$Subcluster, Sobj_integrated$Sample), margin = 2) %>% 
  as.data.frame() %>% as_tibble() %>% 
  separate(Var2, into = c("Type", "B"), extra = "merge", remove = F) %>% select(-B) %>% 
  rename(Subcluster = "Var1", Sample = "Var2") %>% left_join(meta, by = "Sample") %>% 
  
  ggplot(aes(x = Subcluster, y = Freq, fill = Type)) + 
  geom_boxplot(outlier.alpha = 0) +
  geom_point(aes(color = Mens.stat), position = position_dodge(width = 0.2), size = 2) + 
  scale_color_d3(palette = "category20")








