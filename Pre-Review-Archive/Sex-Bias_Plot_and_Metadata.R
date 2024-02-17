collapse.medians <- function() {
  pheno <- 
    GTEx_Analysis_v8_Annotations_SampleAttributesDS %>% 
    select(1, SMTS, SMTSD) %>% 
    separate(SAMPID, into = c("A", "B", "C"), extra = "merge", remove = F) %>% 
    unite(SAMPLEID, A, B, sep = "-") %>% 
    select(-C) %>% 
    left_join(GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS, by = c("SAMPLEID" = "SUBJID")) %>% 
    mutate(SEX = ifelse(SEX == 1, "M", "F"), SAMPID = str_replace_all(SAMPID, "-", "\\."))
  
  
  Females <- pheno %>% filter(SEX == "F") %>% pull(SAMPID)
  Males   <- pheno %>% filter(SEX == "M") %>% pull(SAMPID)
  
  
  Female_Data <- dat.gct %>% 
    select(Description, intersect(Females, colnames(dat.gct))) %>% 
    filter(Description %in% c(filter(Treatment_Response, Cluster == "LUM_HR-pos", abs(avg_logFC) > 0.2)$Gene, "AR")) %>% 
    pivot_longer(intersect(Females, colnames(dat.gct)), names_to = "SAMPID") %>% 
    left_join(pheno)
  
  
  Male_Data <- dat.gct %>% 
    select(Description, intersect(Males, colnames(dat.gct))) %>% 
    filter(Description %in% c(filter(Treatment_Response, Cluster == "LUM_HR-pos", abs(avg_logFC) > 0.2)$Gene, "AR")) %>% 
    pivot_longer(intersect(Males, colnames(dat.gct)), names_to = "SAMPID") %>% 
    left_join(pheno)
  
  
  collapsed_medians     <<- bind_rows(Female_Data, Male_Data) %>% select(Description, SAMPID, value, SMTS)      %>% distinct() %>% group_by(SMTS, Description) %>% summarise(median = median(value))
  collapsed_medians_SEX <<- bind_rows(Female_Data, Male_Data) %>% select(Description, SAMPID, value, SMTS, SEX) %>% distinct() %>% group_by(SEX, SMTS, Description) %>% summarise(median = median(value))
}
build.annotation <- function() {
  A <- GTEx_Analysis_v8_Annotations_SampleAttributesDS %>% select(SMTS, SMTSD) %>% unique() %>% arrange(SMTSD) %>% 
    filter(SMTSD %notin% c("Bladder", "Cervix - Ectocervix", "Cervix - Endocervix", "Fallopian Tube", "Ovary", "Prostate", "Testis", "Uterus", "Vagina", "Cells - Leukemia cell line (CML)", "Kidney - Medulla"))
  B <- effect_size %>% colnames() %>% as_tibble() %>% filter(value != "gene") %>% arrange(value)
  
  collapse_effsize_columns <<- bind_cols(A, B)
  rm(A, B)
  
  ens.anno <<- dat.gct %>% select(1,2)
  anno.dir <<- Treatment_Response %>% filter(Cluster == "LUM_HR-pos") %>% mutate(Dir = ifelse(avg_logFC > 0, "UP", "DO")) %>% select(Description = "Gene", Dir)
}
build.quantiles  <- function() {
  
  test <- effect_size_collapsed %>% filter(Description %in% c(checkUP, checkDO)) %>% mutate(abs = abs(eff_size))
  
  names <- levels
  
  Marker_list    <- vector("list", length = length(names))
  
  for (i in 1:length(names)) {
    
    Marker_list[[i]] <- 
      test %>% filter(SMTS == names[i]) %>% 
      pull(abs) %>% 
      quantile(probs = seq(0, 1, 0.05), na.rm = T) %>% 
      as.data.frame() %>% 
      rownames_to_column() %>% 
      select(tile = "rowname", tile_value = 2) %>% 
      mutate(SMTS = names[i])
    
  }
  
  quantiles <<- do.call(rbind.data.frame, Marker_list)
  
}
pick_genes <- function() {
  
  motif <- 
    NarrowPeaks_DistanceToSummit200_AllCellTypes_Nov16 %>% 
    filter(CellType == "LUM_HR-pos") %>% 
    select(Motif) %>% 
    distinct() %>% 
    filter(str_detect(Motif, paste0("^", query, "_"))) %>% pull(Motif)
  
  
  idx_FC <- NarrowPeaks_DistanceToSummit200_AllCellTypes_Nov16 %>% filter(CellType == "LUM_HR-pos", Motif == motif) %>% select(idxATAC, Signal, Type) %>% distinct() %>% pivot_wider(names_from = Type, values_from = Signal) %>% mutate(Trans = replace_na(Trans, 1), Cis = replace_na(Cis, 1), FC = log(Trans/Cis)) %>% filter(FC > 0.25) %>% pull(idxATAC)
  
  UP <-
    Treatment_Response %>% filter(Cluster == "LUM_HR-pos", avg_logFC > 0.25) %>% pull(Gene)
  
  checkUP <<- 
    NarrowPeaks_DistanceToSummit200_AllCellTypes_Nov16 %>% 
    filter(GeneRna %in% UP, CellType == "LUM_HR-pos", Motif == motif, idxATAC %in% idx_FC) %>% 
    pull(GeneRna) %>% 
    unique() %>% 
    sort()
  
  DO <-
    Treatment_Response %>% filter(Cluster == "LUM_HR-pos", avg_logFC < -0.25) %>% pull(Gene)
  
  checkDO <<- 
    NarrowPeaks_DistanceToSummit200_AllCellTypes_Nov16 %>% 
    filter(GeneRna %in% DO, CellType == "LUM_HR-pos", Motif == motif, idxATAC %in% idx_FC) %>% 
    pull(GeneRna) %>% 
    unique() %>% 
    sort()
  
}

dat.gct <- read.delim(file="~/Downloads/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", skip=2)



Treatment_Response <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Seurat_Objects/Harmony/Treatment_Response_by_CellType.rds") ### treatment response DE genes
GTEx_Analysis_v8_Annotations_SampleAttributesDS  <- read_delim("~/Downloads/SexBias/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", "\t",  escape_double = FALSE, trim_ws = TRUE)
GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS <- read_delim("~/Downloads/SexBias/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
anno_dir <- Treatment_Response %>% filter(Cluster == "LUM_HR-pos") %>% mutate(Direction = ifelse(avg_logFC > 0, "UP", "DO")) %>% select(Description = "Gene", Direction)
NarrowPeaks_DistanceToSummit200_AllCellTypes_Nov16 <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Output/NarrowPeaks_DistanceToSummit200_AllCellTypes_Nov16.rds")

LFSR <-           read_delim("~/Downloads/GTEx_Analysis_v8_sbgenes 2/LFSR.tsv", "\t",           escape_double = FALSE, trim_ws = TRUE)
effect_size <-    read_delim("~/Downloads/GTEx_Analysis_v8_sbgenes 2/effect_size.tsv", "\t",    escape_double = FALSE, trim_ws = TRUE)
effect_size_se <- read_delim("~/Downloads/GTEx_Analysis_v8_sbgenes 2/effect_size_se.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)


query = "CUX2"

pick_genes()



collapse.medians()
build.annotation()

levels <- collapsed_medians %>% filter(Description == "AR", SMTS %in% collapse_effsize_columns$SMTS) %>% arrange(-median) %>% pull(SMTS)
levels_all <- collapsed_medians %>% filter(Description == "AR") %>% arrange(-median) %>% pull(SMTS)

effect_size_collapsed <- 
  ens.anno %>% 
  right_join(effect_size, by = c("Name" = "gene")) %>% 
  select(-Name) %>% pivot_longer(2:45) %>% 
  left_join(collapse_effsize_columns , by = c("name" = "value")) %>% mutate(value = replace_na(value, 0)) %>% 
  group_by(SMTS, Description) %>% 
  summarise(eff_size = mean(value)) %>% mutate(eff_size = replace_na(eff_size, 0))


build.quantiles()

#effsize_mean_difference <- 
#  effect_size_collapsed %>% filter(Description %in% c(checkUP, checkDO)) %>% left_join(anno.dir) %>% group_by(SMTS, Dir) %>% mutate(mean.eff_size = mean(eff_size)) %>% ungroup() %>% select(SMTS, mean.eff_size, Dir) %>% distinct() %>% pivot_wider(names_from = Dir, values_from = mean.eff_size) %>% 
#    mutate(Diff = (DO - UP)) #%>% select(-UP, -DO)
  
  
effsize_mean_difference <- 
  effect_size_collapsed %>% filter(Description %in% c(checkUP, checkDO)) %>% left_join(anno.dir) %>% mutate(Dir = replace_na(Dir, "NA")) %>% group_by(SMTS, Dir) %>% mutate(mean.eff_size = mean(eff_size)) %>% ungroup() %>% select(SMTS, mean.eff_size, Dir) %>% distinct() %>% 
    pivot_wider(names_from = Dir, values_from = mean.eff_size) %>% pivot_longer(c(UP, DO), names_to = "Direction")


df <- 
effect_size_collapsed %>% filter(Description %in% c(checkUP, checkDO, "AR")) %>% left_join(anno.dir) %>% 
  left_join(dplyr::rename(pivot_wider(filter(collapsed_medians, Description %in% c("CUX2", "AR")), names_from = "Description", values_from = "median"), AR_median = "AR", CUX2_median = "CUX2")) %>% 
  group_by(SMTS) %>% mutate(variance = var(eff_size)) %>% left_join(quantiles) %>% left_join(effsize_mean_difference)



AR <- 
df %>% filter(Description == "AR") %>% 
  select(SMTS, AR_median) %>% distinct() %>% 
  ggplot(aes(y = reorder(SMTS, AR_median), x = AR_median)) + 
  geom_col() + ggtitle("AR Expression") +
  theme(axis.text.y = element_text(size = 20), 
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
                     panel.background = element_blank(), 
                     panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
                     panel.grid = element_line(color = "grey90"))


CUX2 <- 
  df %>% filter(Description == "CUX2") %>% 
  select(SMTS, CUX2_median, AR_median) %>% distinct() %>% 
  ggplot(aes(y = reorder(SMTS, AR_median), x = CUX2_median)) + 
  geom_col() + ggtitle("CUX2 Expression") +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
        panel.grid = element_line(color = "grey90"))


CUX2_sex <- 
df %>% filter(Description == "CUX2") %>% 
  ggplot(aes(y = reorder(SMTS, AR_median), x = eff_size)) + 
  geom_col() + ggtitle("CUX2 Sex Bias") +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
        panel.grid = element_line(color = "grey90"))


Quant <- 
df %>% filter(tile %in% c("25%", "50%", "90%")) %>% 
  ggplot(aes(y = reorder(SMTS, AR_median), x = tile_value, fill = tile)) + 
  geom_col(position = position_dodge()) + ggtitle("Quantiles") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), legend.position = "bottom",
        panel.background = element_blank(), 
        panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
        panel.grid = element_line(color = "grey90"))

Varicance <- 
df %>% select(SMTS, variance, AR_median) %>% 
  distinct() %>% 
  ggplot(aes(y = reorder(SMTS, AR_median), x = variance)) + 
  geom_col() + ggtitle("Variance") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), legend.position = "bottom",
        panel.background = element_blank(), 
        panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
        panel.grid = element_line(color = "grey90"))

Eff_Size <- 
df %>% select(SMTS, eff_size, Dir, AR_median) %>% 
  distinct() %>% 
  ggplot(aes(y = reorder(SMTS, AR_median), x = eff_size, fill = Dir)) + 
  geom_boxplot(outlier.alpha = 0.35) + ggtitle("Effect Size in UP/DOWN") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), legend.position = "bottom",
        panel.background = element_blank(), 
        panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
        panel.grid = element_line(color = "grey90"))


#Diff <- 
#df %>% select(SMTS, Diff, AR_median) %>% 
#  distinct() %>% 
#  ggplot(aes(y = reorder(SMTS, AR_median), x = Diff)) + 
#  geom_col() + ggtitle("abs eff_size mean diff") +
#  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
#        axis.ticks.y = element_blank(), legend.position = "bottom",
#        panel.background = element_blank(), 
#        panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
#        panel.grid = element_line(color = "grey90"))


Dir <- 
  df %>% select(SMTS, Direction, value, AR_median) %>% 
  distinct() %>% 
  ggplot(aes(y = reorder(SMTS, AR_median), x = value, fill = Direction)) + 
  geom_col(position = position_dodge()) + 
  ggtitle("eff_size mean") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), legend.position = "bottom",
        panel.background = element_blank(), 
        panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
        panel.grid = element_line(color = "grey90"))




AR + Quant + Varicance + Dir + Eff_Size + patchwork::plot_layout(ncol = 5, widths = c(1,1,1,1,3))

AR + CUX2 + CUX2_sex + Quant + Varicance + Dir + Eff_Size + patchwork::plot_layout(ncol = 7, widths = c(1,1,1,1,1,1,3))


p %>% save_x(data = ., name = "IFNG_Response_Immune", 1, 16, 9, svg = T)    











#### sex bias heatmaps

filler <- collapsed_medians %>% pull(SMTS) %>% unique() %>% as_tibble() %>% mutate(UP = 0, DO = 0, `25%` = 0,  `50%` = 0,  `90%` = 0) %>% dplyr::rename(SMTS = value)


AR_median <- collapsed_medians %>% filter(Description == "AR") %>% ggplot(aes(y = reorder(SMTS, median), x = median)) + geom_col()

AR_median %>% save_x(data = ., name = "AR_Median", 1, 10, 7, svg = T)

mat <- df %>% filter(tile %in% c("25%", "50%", "90%")) %>% select(SMTS, tile, tile_value) %>% unique() %>% pivot_wider(names_from = tile, values_from = tile_value) %>% arrange(factor(SMTS, levels = levels)) %>% bind_rows(select(filter(filler, SMTS %notin% levels), -UP, -DO)) %>% arrange(factor(SMTS, levels = levels_all)) %>% column_to_rownames("SMTS")

AR_quantiles <- pheatmap(mat, cluster_rows = F, cluster_cols = F, color = viridis::viridis(n = 100))



mat <- df %>% select(SMTS, Direction, value) %>% distinct() %>% pivot_wider(names_from = Direction, values_from = value) %>% bind_rows(select(filter(filler, SMTS %notin% levels), -`25%`, -`50%`, -`90%`)) %>% arrange(factor(SMTS, levels = levels_all)) %>% column_to_rownames("SMTS")
my.breaks <- c(seq((min(mat)), -0.00001, length.out = 50), seq(0.00001, (max(mat)), length.out = 50))
my.colors <- c(colorRampPalette(colors = c("#002c2d", "#00494b", "#076769", "#208288", "#3ea8a6", "#76c7be", "#bbe4d1", "#ffffe0"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("#ffffe0", "#fed693", "#f5ad52", "#dd8629", "#bd651a", "#9c4511", "#79260b", "#580000"))(length(my.breaks)/2))

AR_Bias <- pheatmap(mat, cluster_rows = F, cluster_cols = F, color = my.colors, breaks = my.breaks)




p <- as.ggplot(AR_quantiles) + as.ggplot(AR_Bias)


p %>% save_x(data = ., name = "AR_quant_Effsize", 1, 10, 7, svg = T)



tissue_enrichment_data_AllSamples <- read_csv("~/Downloads/tissue_enrichment_data_AllSamples.csv")

mat <- tissue_enrichment_data_AllSamples %>% select(-count) %>% arrange(factor(tissue, levels = levels_all)) %>% column_to_rownames("tissue")
my.breaks <- c(seq((min(mat)), 0.499, length.out = 50), seq(0.5, (max(mat)), length.out = 50))
my.colors <- c(colorRampPalette(colors = c("#002c2d", "#00494b", "#076769", "#208288", "#3ea8a6", "#76c7be", "#bbe4d1", "#ffffe0"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("#ffffe0", "#fed693", "#f5ad52", "#dd8629", "#bd651a", "#9c4511", "#79260b", "#580000"))(length(my.breaks)/2))

pheatmap(mat, cluster_rows = F, cluster_cols = F, color = my.colors, breaks = my.breaks)



mat <- 
  tissue_enrichment_data_AllSamples_1_ %>% select(-count, -quantile) %>% arrange(factor(tissue, levels = levels_all)) %>% mutate(fc_vs_median = str_replace(fc_vs_median, "Inf", "1"), fc_vs_median = log(as.numeric(fc_vs_median))) %>% column_to_rownames("tissue")

my.breaks <- c(seq((min(mat)), -0.00001, length.out = 50), seq(0.00001, (max(mat)), length.out = 50))
my.colors <- c(colorRampPalette(colors = c("#002c2d", "#00494b", "#076769", "#208288", "#3ea8a6", "#76c7be", "#bbe4d1", "#ffffe0"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("#ffffe0", "#fed693", "#f5ad52", "#dd8629", "#bd651a", "#9c4511", "#79260b", "#580000"))(length(my.breaks)/2))

p <- pheatmap(mat, cluster_rows = F, cluster_cols = F, color = my.colors, breaks = my.breaks)


p %>% save_x(data = ., name = "Tissue_Enrichment_AR", 1, 10, 7, svg = T)









##### CUX2

CUX2_median <- 
  collapsed_medians %>% filter(Description == "CUX2") %>% ggplot(aes(y = factor(SMTS, levels = rev(levels_all)), x = median)) + geom_col()

CUX2_median %>% save_x(data = ., name = "CUX2_Median", 1, 10, 7, svg = T)

mat <- df %>% filter(tile %in% c("25%", "50%", "90%")) %>% select(SMTS, tile, tile_value) %>% unique() %>% pivot_wider(names_from = tile, values_from = tile_value) %>% arrange(factor(SMTS, levels = levels)) %>% bind_rows(select(filter(filler, SMTS %notin% levels), -UP, -DO)) %>% arrange(factor(SMTS, levels = levels_all)) %>% column_to_rownames("SMTS")

CUX2_quantiles <- pheatmap(mat, cluster_rows = F, cluster_cols = F, color = viridis::viridis(n = 100))




mat <- df %>% select(SMTS, Direction, value) %>% distinct() %>% pivot_wider(names_from = Direction, values_from = value) %>% bind_rows(select(filter(filler, SMTS %notin% levels), -`25%`, -`50%`, -`90%`)) %>% arrange(factor(SMTS, levels = levels_all)) %>% column_to_rownames("SMTS")
my.breaks <- c(seq((min(mat)), -0.00001, length.out = 50), seq(0.00001, (max(mat)), length.out = 50))
my.colors <- c(colorRampPalette(colors = c("#002c2d", "#00494b", "#076769", "#208288", "#3ea8a6", "#76c7be", "#bbe4d1", "#ffffe0"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("#ffffe0", "#fed693", "#f5ad52", "#dd8629", "#bd651a", "#9c4511", "#79260b", "#580000"))(length(my.breaks)/2))

CUX2_Bias <- pheatmap(mat, cluster_rows = F, cluster_cols = F, color = my.colors, breaks = my.breaks)







mat <- df %>% filter(Description == "CUX2") %>% select(SMTS, eff_size) %>% unique() %>% bind_rows(select(filter(filler, SMTS %notin% levels), -`25%`, -`50%`, -`90%`)) %>% arrange(factor(SMTS, levels = levels_all)) %>% mutate(eff_size = replace_na(eff_size, 0)) %>% select(-UP, -DO) %>% column_to_rownames("SMTS")

my.breaks <- c(seq((min(mat)), -0.00001, length.out = 50), seq(0.00001, (max(mat)), length.out = 50))
my.colors <- c(colorRampPalette(colors = c("#002c2d", "#00494b", "#076769", "#208288", "#3ea8a6", "#76c7be", "#bbe4d1", "#ffffe0"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("#ffffe0", "#fed693", "#f5ad52", "#dd8629", "#bd651a", "#9c4511", "#79260b", "#580000"))(length(my.breaks)/2))

CUX2_effsize <- pheatmap(mat, cluster_rows = F, cluster_cols = F, color = my.colors, breaks = my.breaks)




p <- as.ggplot(CUX2_quantiles) + as.ggplot(CUX2_Bias) + as.ggplot(CUX2_effsize)


p %>% save_x(data = ., name = "CUX2_quant_Effsize", 1, 15, 7, svg = T)







df1 <- collapsed_medians %>% filter(Description == "AR")
p <- 
  effect_size_collapsed %>% filter(Description == "CUX2") %>% 
  left_join(df1, by = "SMTS") %>% ggplot(aes(x = eff_size, y = median)) + 
  geom_point(aes(size = median)) + 
  geom_text_repel(aes(label = SMTS)) + 
  xlim(c(-1.7, 1.7)) + 
  stat_smooth(method = "lm", size=1)

p %>% save_x(data = ., name = "CUX2_Effsize_vs_AR", 1, 15, 7, svg = T)












mat <- effect_size_collapsed %>% filter(Description %in% checkUP) %>% pivot_wider(names_from = Description, values_from = eff_size) %>% arrange(factor(SMTS, levels = levels)) %>% column_to_rownames("SMTS")

my.breaks <- c(seq((min(mat)), -0.00001, length.out = 50), seq(0.00001, (max(mat)), length.out = 50))
my.colors <- c(colorRampPalette(colors = c("#002c2d", "#00494b", "#076769", "#208288", "#3ea8a6", "#76c7be", "#bbe4d1", "#ffffe0"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("#ffffe0", "#fed693", "#f5ad52", "#dd8629", "#bd651a", "#9c4511", "#79260b", "#580000"))(length(my.breaks)/2))


m1 <- pheatmap(mat, breaks = my.breaks, color = my.colors, show_colnames = F, cluster_rows = F, border_color = NA, 
         clustering_distance_cols = "euclidean", fontsize_row = 20, main = "UP Genes", show_rownames = F)



mat <- effect_size_collapsed %>% filter(Description %in% checkDO) %>% pivot_wider(names_from = Description, values_from = eff_size) %>% arrange(factor(SMTS, levels = levels)) %>% column_to_rownames("SMTS")

my.breaks <- c(seq((min(mat)), -0.00001, length.out = 50), seq(0.00001, (max(mat)), length.out = 50))
my.colors <- c(colorRampPalette(colors = c("#002c2d", "#00494b", "#076769", "#208288", "#3ea8a6", "#76c7be", "#bbe4d1", "#ffffe0"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("#ffffe0", "#fed693", "#f5ad52", "#dd8629", "#bd651a", "#9c4511", "#79260b", "#580000"))(length(my.breaks)/2))


m2 <- pheatmap(mat, breaks = my.breaks, color = my.colors, show_colnames = F, cluster_rows = F, border_color = NA, 
         clustering_distance_cols = "euclidean", fontsize_row = 20, main = "DOWN Genes")


p <- cowplot::plot_grid(as.ggplot(m1), as.ggplot(m2), rel_widths = c(1,1.5), ncol = 2)

p %>% save_x(data = ., name = "Sex_Bias_Heatmap", 1.5, 12, 9, svg = T)    











dflong <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Output/Reboot/GTEX_MeadianTPM_long.rds")
signif_sbgenes <- read_delim("~/Downloads/GTEx_Analysis_v8_sbgenes 2/signif.sbgenes.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
TPM <- dflong %>% pull(Tissue) %>% unique() %>% sort() %>% as_tibble() %>% filter(value %notin% exclude) %>% dplyr::rename(TPM = "value")
anno <- dflong %>% select(ID, hgnc_symbol)

assign <- signif_sbgenes %>% 
  pull(tissue) %>% 
  unique() %>% 
  sort() %>% 
  as_tibble() %>% 
  select(SB = 1) %>% 
  bind_cols(TPM) %>% 
  print(n = 50)



levels <- dflong %>% left_join(assign, by = c("Tissue" = "TPM")) %>% 
  filter(!is.na(SB), !str_detect(SB, "Cells_")) %>% select(-1, -2, -Tissue) %>% 
  filter(hgnc_symbol == "AR") %>% arrange(-value) %>% pull(SB)













query = "AR"

motif <- 
  NarrowPeaks_DistanceToSummit200_AllCellTypes_Nov16 %>% 
  filter(CellType == "LUM_HR-pos") %>% 
  select(Motif) %>% 
  distinct() %>% 
  filter(str_detect(Motif, paste0("^", query, "_"))) %>% pull(Motif)


idx_FC <- NarrowPeaks_DistanceToSummit200_AllCellTypes_Nov16 %>% filter(CellType == "LUM_HR-pos", Motif == motif) %>% select(idxATAC, Signal, Type) %>% distinct() %>% pivot_wider(names_from = Type, values_from = Signal) %>% mutate(Trans = replace_na(Trans, 1), Cis = replace_na(Cis, 1), FC = log(Trans/Cis)) %>% filter(FC > 0.25) %>% pull(idxATAC)

UP <-
  Treatment_Response %>% filter(Cluster == "LUM_HR-pos", avg_logFC > 0.25) %>% pull(Gene)

checkUP <- 
  NarrowPeaks_DistanceToSummit200_AllCellTypes_Nov16 %>% 
  filter(GeneRna %in% UP, CellType == "LUM_HR-pos", Motif == motif, idxATAC %in% idx_FC) %>% 
  pull(GeneRna) %>% 
  unique() %>% 
  sort()

DO <-
  Treatment_Response %>% filter(Cluster == "LUM_HR-pos", avg_logFC < -0.25) %>% pull(Gene)

checkDO <- 
  NarrowPeaks_DistanceToSummit200_AllCellTypes_Nov16 %>% 
  filter(GeneRna %in% DO, CellType == "LUM_HR-pos", Motif == motif, idxATAC %in% idx_FC) %>% 
  pull(GeneRna) %>% 
  unique() %>% 
  sort()








































