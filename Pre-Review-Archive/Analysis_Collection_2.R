# Library Loading ---------------------------------------------------------
library(Seurat)
library(ggsci)
library(scales)
library(tidyverse)
library(ggpubr)
library(EnhancedVolcano)
library(pheatmap)
library(readxl)
library(readr)
library(RColorBrewer)
library(cowplot)
library(ggplotify)

setwd("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/")

col <- pal_d3("category20")(20)
options(tibble.print_max = 150, tibble.print_min = 35)
`%notin%` <- Negate(`%in%`)




# function to save plot named "image"
save     <- function(name, scale) {
  
  dir.create(paste0("Figures/5.13_CellType_Vignettes/", celltype, "/PNG"), recursive = T, showWarnings = F)
  dir.create(paste0("Figures/5.13_CellType_Vignettes/", celltype, "/SVG"), recursive = T, showWarnings = F)
  
  png.path <- paste0("Figures/5.13_CellType_Vignettes/", celltype, "/PNG")
  svg.path <- paste0("Figures/5.13_CellType_Vignettes/", celltype, "/SVG")
  
  ggsave(plot   = image,
         path     = svg.path,
         filename = paste0(name,".svg"),
         width    = 16, 
         height   = 9,
         scale    = scale,
         dpi      = 300)
  
  ggsave(plot   = image,
         path     = png.path,
         filename = paste0(name,".png"),
         width    = 16, 
         height   = 9,
         scale    = scale,
         dpi      = 300)
}
save2    <- function(name, scale) {
  
  ggsave(plot     = quick,
         path     = figure.path,
         filename = paste0(name,".png"),
         width    = 16, 
         height   = 9,
         scale    = scale,
         dpi      = 300)
}

sig.read <- function(n) {
  Signalingpathways_AR <- read_excel("utilities/Signalingpathways_AR.xlsx")  %>% dplyr::select(Family = "Family", Gene = "Gene", Percentile = "Percentile", Disc.Rate = `Discovery Rate`, GMFC = "GMFC") %>% mutate(Family = "AR") %>% head(n)
  SignalingpathwaysESR <- read_excel("utilities/Signalingpathways_ESR.xlsx") %>% dplyr::select(Family = "Family", Gene = "Gene", Percentile = "Percentile", Disc.Rate = `Discovery Rate`, GMFC = "GMFC") %>% mutate(Family = "ESR") %>% head(n)
  SignalingpathwaysPGR <- read_excel("utilities/Signalingpathways_PGR.xlsx") %>% dplyr::select(Family = "Family", Gene = "Gene", Percentile = "Percentile", Disc.Rate = `Discovery Rate`, GMFC = "GMFC") %>% mutate(Family = "PGR") %>% head(n)
  Sign <- rbind(Signalingpathways_AR, SignalingpathwaysESR, SignalingpathwaysPGR)
  
  print(Sign)
  return(Sign)
  
}






# Load DataSets -----------------------------------------------------------


# load Seurat Obejct
Sobj_integrated <- readRDS("Seurat_Objects/Harmony/Sobj_All_Cells_Sep28.rds")
Sobj_integrated <- readRDS("Seurat_Objects/Harmony/CellType/Harmony_Myeloid.rds")
Sobj_integrated <- readRDS("Seurat_Objects/Harmony/Group/Harmony_Batch_Vasculature.rds")

DefaultAssay(Sobj_integrated) <- "RNA"
Idents(Sobj_integrated)       <- "CellType"
Idents(Sobj_integrated)       <- "Subcluster"
levels(Sobj_integrated)       <- c("LUM_HR-pos", "LUM_HR-neg", "Basal", "Fibroblast", "Adipocyte", "Blood_EC", "Lymphatic_EC", "Vasc_Accessory", "Myeloid", "Lymphoid")

cols <- c("#1BA8B3", "#24948C", "#2F6282", "#71B340", "#F6D579", "#DB3B58", "#A23B72", "#F29CA3", "#8459DE", "#8088E6")

levels <- c("LUM_HR-pos", "LUM_HR-neg", "Basal", "Fibroblast", "Adipocyte", "Blood_EC", "Lymphatic_EC", "Vasc_Accessory", "Myeloid", "Lymphoid")
cols   <- c("#44aa99", "#88ccee", "#117733", "#332288", "#ddcc77", "#882255", "#aa4499", "#cc6677", "#999933", "#dddddd")

levels_samp <- c("CF-3920", "CF-19301", "CF-7780", "CF-2797", "CF-318-813", "CF-0404", "CF-428-112", "CF-249-347", "CF-4014", "TM-9469", "TM-1956", "TM-6544", "TM-6477", "TM-8249", "TM-7567", "TM-9817", "TM-2768", "TM-3937")
cols_samp   <- c("#F26C66", "#FFBB78", "#F57F46", "#7FC210", "#2CA02C", "#BCBD22", "#EDBD1F", "#A8486A", "#C49C94", "#9467BD", "#C5B0D5", "#17BECF", "#C5D4D2", "#438ABB", "#FA7AA3", "#F7B6D2", "#FF9896", "#3FB8AF")

cols <- c("#9d3396", "#eb933e")

Treatment_Response <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Seurat_Objects/Harmony/Treatment_Response_by_CellType.rds") ### treatment response DE genes
Treatment_Response_Subcluster <- readRDS("Seurat_Objects/Harmony/Treatment_Response_By_Subcluster.rds")

LUMexpr            <- do.call(rbind.data.frame, readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Seurat_Objects/Reboot/Expression_list.rds")) %>% ungroup() %>% filter(expression == "expressed") ### percentage of cell expressing all genes
cluster.averages   <- readRDS("Output/Reboot/TM-CF_cluster_averages_across_CellTypes.rds") %>% select(-avg_logFC, -p_val_adj) %>% left_join(Treatment_Response, by = c("Gene", "Cluster")) %>% select(-p_val, -pct.1, -pct.2) ### average gene expression across TM/CF
Sample_averages    <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Output/Reboot/Average_Expression_bySample_CellType.rds")
Sample_averages_Subcluster    <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Output/Reboot/Average_Expression_bySample_Subcluster.rds")
TFs                <- read_tsv("~/Documents/SCENIC/pySCENIC_Full/scenicdata/Resources/hs_hgnc_tfs.txt", col_names = "Gene") %>% mutate(TF = "YES")


Commmon_up <- Treatment_Response %>% filter(avg_logFC > 0.25) %>% group_by(Gene) %>% summarise(n = n()) %>% filter(n >= 5) %>% pull(Gene)
Commmon_do <- Treatment_Response %>% filter(avg_logFC < -0.25) %>% group_by(Gene) %>% summarise(n = n()) %>% filter(n >= 5) %>% pull(Gene)
common <- c(Commmon_up, Commmon_do) 
rm(Commmon_up, Commmon_do)




percentile <- 
  read_tsv(paste0("Output/Reboot/SCENIC/expr_mat.adjacencies_all_cells.tsv")) %>% 
  group_by(TF) %>% 
  mutate(percentile_rank = ntile(importance, 100))

percentile <- 
  read_tsv(paste0("Output/Reboot/SCENIC/CellType_SCENIC/Output/expr_mat.adjacencies_Blood_EC.tsv")) %>% 
  group_by(TF) %>% 
  mutate(percentile_rank = ntile(importance, 100))




Adipo_Unique2 <- c("GPAM", "SOX6", "ZCCHC14", "RUFY3", "CTBP2", "BDP1", "KMT2A", "BPTF", "PPARG")
Adipo_Unique2 <- c("leptin")


for (i in 1:length(Adipo_Unique2)) {
  

#genes <- percentile %>% filter(TF == Adipo_Unique2[i], percentile_rank > 95) %>% pull(target)

genes <- gmt[Adipo_Unique2[i]]
  
Sobj_integrated <- AddModuleScore  (Sobj_integrated, 
                                    assay = "RNA", 
                                    #list(filter(Corr_Bind, str_detect(Clust.Name, "^KLK") | str_detect(Clust.Name, "^AR"))$Gene), 
                                    list(genes), 
                                    nbin = 25, name = paste0(Adipo_Unique2[i],"_"))

}

# DimPlot -----------------------------------------------------------------

p1 <- 
DimPlot(Sobj_integrated, 
        #cells = WhichCells(Sobj_integrated, idents = c("TM-3937", "CF-2797")),
        reduction = "umap", 
        group.by  = "Subcluster",
        #group.by  = "Sample", split.by = "Type",
        #cells = WhichCells(Sobj_integrated, downsample = 100), 
        #group.by  = "Sample",
        #split.by = "Type",
        pt.size = 1, 
        label = T,
        repel = T, 
        #ncol = 5,
        label.size = 6,
        #split.by  = "Sample", ncol = 6, 
        #cols = viridis::magma(n = 2, begin = 0.3, end = 0.8),
        #cols = col
        ) + theme(axis.title = element_blank(),
                  legend.position = "bottom", 
                  text =  element_text(family = "Lato"), 
                  #panel.grid = element_line(color = "grey60")
                  ) 
  

plot2 %>% save_x(data = ., name = "test", 1, 10, 10, svg = T)



plot1 + plot2 + patchwork::plot_layout(ncol = 2)


#save2("Immune_Markers", 1.5)


# VlnPLot -----------------------------------------------------------------

query <- c("TLL1", "LNX1", "PCSK5", "RPL13", "RELN", "RADIL", "THSD7B", "CLSTN2", "CACNB2")

#P1 <-  
VlnPlot(Sobj_integrated, log = T, 
            #idents = levels(Sobj_integrated)[1:4],
            features = c("THBS1"),
            group.by = "Subcluster", 
            #split.by = "Type",
            pt.size = 0, 
            ncol = 1, #cols = cols,
            #cols = viridis::magma(n = 2, begin = 0.3, end = 0.8)
        )

p %>% save_x(data = ., name = paste0("Fibroblast_JUNB"), 1, 5, 10, svg = T)


save(name = paste0("VlnPlot_", celltype), 0.8, 5, 10, svg = F)        

# Featureplot -------------------------------------------------------------

features = c("TP63", "TEAD1")

#P2 <- 
FeaturePlot(Sobj_integrated, reduction = "umap",
            #cells = WhichCells(Sobj_integrated, expression = SCORE_1 >= 0.25),
            #features = c("PLXNA4", "HPSE2", "SEMA3C", "EPHA3", "SLC22A3", "ZBTB7C"), 
            features = c("ALDH1A1"),
            pt.size = 1, 
            order = T, 
            #max.cutoff = 4.2,
            cols = viridis::magma(n = 100),
            #cols = c("grey", "Green", "Red")
            ncol = 1
            ) #+ look

save(name = paste0("FeaturePlot2_", celltype), 1, 10, 7, svg = F)


quick1 + quick2 + plot_layout(ncol = 1)



Adipo_Unique2 <- 
  c("CUX2", "ZNF689", "SPDEF", "GATA3", "TBX3", 
    "PGR", "PBX3", "BATF", "GLIS3", "PRDM16", "BCL11A", "RUNX1",
    "GLI3", "EZR", "RFX2", "ATF3", "BACH1", "AHR", "BHLHE40", "HIF1A",
    "ESR1", "AR", "TRPS1", "TFAP2A",
    "ZNF652", "CERS4", "LTF", "ESRRG", "TFAP2B", 
    "NCALD", "RARB", "CREB5", "SOX5", "HES1", "ID2", 
    "RPL35", "RPS4X", 
    "DTL", "ZNF367",
    "NR4A1", "HSPA5")

features3 <- paste0(features3, "_1")




query <- c("PRKAR2B","PLIN1","CIDEA","PPARG","LPL","PNPLA2","ADCY2","FABP4","CAV1","ABHD5", "LIPE-AS1", "PLIN2", "PLIN5")


# Sample_Avg Plots -------------------------------------------------------------


query = c("INSR")
#query <- Treatment_Response %>% filter(Cluster == "LUM_HR-neg", pct.2 > 0.3) %>% top_n(50, avg_logFC) %>% pull(Gene)

#p <-
Sample_averages %>% filter(Gene %in% query, #Sample != "TM-9817",
                           #Cluster %in% c("Adipocyte"), 
                           ) %>% 
  ggplot(aes(x = Cluster, y = avg_Expr, fill = Type)) +
  
  #geom_violin(alpha = 0.65, trim = F, scale = "width" , draw_quantiles = T) +
  
  geom_boxplot(color = "grey30", 
               width = 0.5, 
               alpha = 0.8, 
               position=position_dodge(0.9), 
               lwd = 1,
               outlier.alpha = 0) +
  
  geom_point(aes(color = factor(Sample, levels = levels_samp)),
             #position = position_dodge(width = 0.75), 
             position = position_jitterdodge(jitter.width = 0.000005, 
                                             jitter.height = 0,
                                             dodge.width = 0.5), 
             alpha = 1,
             size  = 3) +
  
  #stat_compare_means(paired = F, method = "wilcox", 
  #                   label.x.npc = 0.25, vjust = 1) +
  
  scale_fill_viridis_d(option = "C", begin = 0.3, end = 0.8) +
  #scale_color_viridis_d(begin = 0.3, end = 0.6) +
  scale_color_manual(values = cols_samp) +
  
  ggtitle(label = "Sample Average Expression") +
  
  theme(text = element_text(family = "Lato", size = 20),
        axis.text = element_text(), #legend.title = element_text("test"),
        #axis.title=element_text(size=14,face="bold"),
        axis.title = element_blank(), #legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey30", size = 2),
        panel.grid = element_blank(),
        strip.text = element_text(size = 15),
        axis.text.x = element_text(size=14, angle = -25, hjust = 0, vjust = 1)
        ) +
  #scale_y_sqrt() +
  #scale_y_log10() +
  facet_wrap(vars(Gene), scales = "free_y", ncol = 2)



p %>% save_x(data = ., name = "PPARG_Average_Expression", 1, 8, 10, svg = T)







Modules <- Sobj_integrated@meta.data %>% filter(CellType == "Adipocyte") %>% select(Sample, Type, query) %>% rownames_to_column("name")
Module_Avg <- Modules %>% group_by(Sample) %>% mutate_at(.vars = query, mean) %>% select(-name) %>% unique() %>% pivot_longer(query)%>% unique()


query = c("GPAM_1")

Module_Avg %>% filter(Sample %in% enough) %>% 
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
  
  ggtitle(label = "Sample Average Expression") +
  
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=14,face="bold"), 
        axis.text.x = element_text(size=10)) +
  #scale_y_sqrt() +
  facet_wrap(vars(name), scales = "free_y")






A <- 
Treatment_Response %>% 
  filter(Cluster == "Adipocyte", 
         !str_detect(Gene, "-AS"), 
         !str_detect(Gene, "\\.[1-5]")) %>% 
  top_n(50, avg_logFC) %>% 
  arrange(-avg_logFC)

B <- 
  Treatment_Response %>% 
  filter(Cluster == "Adipocyte", 
         !str_detect(Gene, "-AS"), 
         !str_detect(Gene, "\\.[1-5]")) %>% 
  top_n(-50, avg_logFC) %>% 
  arrange(-avg_logFC)


bind <- bind_rows(A, B)


DotPlot(Sobj_integrated, features = bind$Gene, group.by = "Sample") + coord_flip()





# Cluster Proportions Boxplot ---------------------------------------------



levels(Sobj_integrated) <- table(Sobj_integrated$Subcluster) %>% as.data.frame() %>% arrange(-Freq) %>% pull(Var1)


Prop <- 
  as.data.frame(prop.table(table(Sobj_integrated$Subcluster, 
                                 Sobj_integrated$Mens.stat), 
                           margin = 2)) %>% 
  dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq") %>% 
  mutate(Type = substring(.$Sample, 1, 2))# %>% arrange(factor(Cluster, levels = levels), desc(Freq))



#Test <- 
  as.data.frame(table(Sobj_integrated$Subcluster, 
                      Sobj_integrated$Sample)) %>% 
  dplyr::select(Cluster = "Var1", Sample = "Var2", Cells = "Freq") %>% 
  mutate(Type = substring(.$Sample, 1, 2)) %>% 
  group_by(Cluster) %>% 
  mutate(Total = sum(Cells)) %>% 
  mutate(Freq = (Cells/Total)*100) %>% 
  group_by(Cluster) %>% 
  mutate(Cutoff_2.5 = ifelse(Freq >= 2.5, "pass", "fail")) %>% 
  #filter(Cluster == 6) %>% arrange(-Freq)
  #mutate(Cutoff_5   = ifelse(Freq >= 5,   "pass", "fail")) %>% 
  
  group_by(Cluster, Cutoff_2.5) %>%
  summarise(n = n()) %>% 
  group_by(Cluster) %>% 
  mutate(Inclusion = ifelse(n >= 6, "pass", "not_pass")) %>% 
  filter(Cutoff_2.5 == "pass") %>% 
  select(RNA_snn_res.0.2 = "Cluster", Contr.Samples = "n", 4) %>% ungroup() %>% 
  right_join(rownames_to_column(Sobj_integrated@meta.data, "ID"), by = "RNA_snn_res.0.2") %>% 
  select(1:4) %>% column_to_rownames("ID")


Sobj_integrated <- MetaData(Sobj_integrated, Test)

filter(Cluster == 4) %>% 
  arrange(-Freq)



#P2 <- 
#Prop %>% filter(!is.na(Freq), Cluster %in% c("Macrophage", "moDC", "Monocyte", "DC")) %>% 
Prop %>% #filter(!is.na(Freq), Cluster %in% c("T-Effector", "CD4_T", "CD8_T", "NK")) %>% 
ggplot(aes(x = Cluster, y = Freq, fill = Type)) + 
  
  #geom_violin(alpha = 0.65, trim = F, scale = "width" , draw_quantiles = T) +
  
  geom_boxplot(color = "grey40", width = 0.5, alpha = 1, lwd = 1.25,
               position=position_dodge(0.75), outlier.alpha = 0) +
  
  geom_point(aes(color = Sample), 
             position = position_jitterdodge(jitter.width = 0.5, 
                                             jitter.height = 0,
                                             dodge.width = 0.4), 
             alpha = 1,
             size  = 4) +
  scale_y_sqrt() + 

  scale_color_manual(values =  cols) +
  
  ggtitle("Fibroblast Proportions") + 
  scale_fill_manual(values = c("#9d3396", "#eb933e")) +
  theme(text = element_text(family = "Lato", size = 25),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        #axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        panel.background = element_blank(), 
        panel.border = element_rect(color = "grey30", size = 2, fill = NA), 
        ) 



p %>% save_x(data = ., name = paste0("Basal_Subcluster_Proportions_Polar"), 1, 16, 9, svg = T)




ggplot(data = Prop, aes(x = "", y = Freq, fill = Cluster )) + 
  geom_bar(stat = "identity", position = position_fill()) +
  #geom_text(aes(label = Freq), position = position_fill(vjust = 0.5)) +
  coord_polar(theta = "y") +
  
  #scale_color_manual(values = cols_cell) +
  
  #facet_wrap(~ Sample, nrow = 3)  #+
  theme(axis.title.x = element_blank(),
      axis.title.y = element_blank()) + 
  theme(legend.position='bottom') + 
  guides(fill=guide_legend(nrow=2, byrow=TRUE))







Prop <- 
  table(Sobj_integrated$Subcluster, Sobj_integrated$Type) %>% 
  as.data.frame() %>% rename(Type = "Var2")
  

p <- 
ggplot(data = Prop, aes(x = "", y = Freq, fill = Type )) + 
  geom_bar(stat = "identity", position = position_fill()) +
  #geom_text(aes(label = Freq), position = position_fill(vjust = 0.5)) +
  coord_polar(theta = "y") +
  
  scale_fill_viridis_d(option = "C", begin = 0.4, end = 0.8) +
  
  facet_wrap(~ Cluster, nrow = 1)  +
  
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        panel.background = element_blank()) #+ 
  
  guides(fill=guide_legend(nrow=2, byrow=TRUE))






  
  
mat <- Bind %>% select(Gene, Cluster, avg_logFC) %>% pivot_wider(names_from = Cluster, values_from = avg_logFC) %>% column_to_rownames("Gene")
  
mat[is.na(mat)] <- 0

image <- 
pheatmap(t(mat), 
         color = viridis::magma(n = 100, begin = 0, end = 0.9),
         main = "logFC - GO:BP Golgi Organization", 
         fontsize_row = 15, fontsize_col = 8,
         treeheight_row = 5, treeheight_col = 15
          )

save(name = "GO:BP_Golgi_Organization", 1, 16, 5)




  
mat <- percentile %>% 
  filter(TF %in% Adipo_Unique, target %in% Adipo_Genes$Gene,  percentile_rank > 90) %>% 
  select(-percentile_rank) %>% 
  pivot_wider(names_from = target, values_from = importance) %>% column_to_rownames("TF")







#image <- 
pheatmap(mat, 
         scale = "column", 
         fontsize = 10, 
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         #color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)),
         #annotation_row = Ann,
         #annotation_colors = ann_colors, 
         treeheight_row = 0, treeheight_col = 0, 
         #cutree_rows = 1, 
         #cutree_cols = 2,
         fontsize_row = 10,
         fontsize_col = 0.0010, 
         labels_col = NULL, 
         color = viridis::magma(n = 100, begin = 0, end = 0.9), 
         main = "Main"
)







Adipo_Sig <- 
Sample_averages %>% 
  filter(Cluster %in% "Adipocyte", Gene %in% genes) %>% 
  group_by(Gene) %>%
  do(w = wilcox.test(avg_Expr~Type, data = ., paired=FALSE)) %>% 
  summarise(Gene, Wilcox = w$p.value) %>% filter(Wilcox <= 0.05) %>% arrange(Wilcox) %>% add_column(CellType = input_vargenes)




#image <- 
EnhancedVolcano(df, df$Gene, 
                x = "avg_logFC",
                y = "p_val_adj", transcriptPointSize = 5, transcriptLabSize = 5,
                #xlim = c(-0.50, 0.5),
                FCcutoff = 0.25, 
                #title = "NR3C1 - High Confidence Targets (Basal)"
                )



save(name = "NR3C1 - High Confidence Targets (Basal)", 1, 12, 9, svg = F)






library(Seurat)
library(tidyverse)
Sobj_integrated <- readRDS("~/scNuclei-Libraries/Analysis/Reboot/Seurat_Objects/Final_Integration/Sobj_GoodSamples_Integrated_Final.rds")

percentile <- 
  read_tsv("expr_mat.adjacencies.tsv") %>% 
  group_by(TF) %>% 
  mutate(percentile_rank = ntile(importance, 100))


TFs <- percentile %>% ungroup() %>% distinct(TF) %>% pull(TF)

for (i in 1:length(TFs)) {
  
  genes <- percentile %>% filter(TF == TFs[i], percentile_rank > 90) %>% pull(target)
  
  Sobj_integrated <- AddModuleScore  (Sobj_integrated, 
                                      assay = "RNA", 
                                      list(genes), 
                                      nbin = 25, name = paste0(TFs[i], "_90th_"))
}



for (i in 1:length(TFs)) {
  
  genes <- percentile %>% filter(TF == TFs[i], percentile_rank > 95) %>% pull(target)
  
  Sobj_integrated <- AddModuleScore  (Sobj_integrated, 
                                      assay = "RNA", 
                                      list(genes), 
                                      nbin = 25, name = paste0(TFs[i], "_95th_"))
}


for (i in 1:length(TFs)) {
  
  genes <- percentile %>% filter(TF == TFs[i], percentile_rank > 99) %>% pull(target)
  
  Sobj_integrated <- AddModuleScore  (Sobj_integrated, 
                                      assay = "RNA", 
                                      list(genes), 
                                      nbin = 25, name = paste0(TFs[i], "_99th_"))
}



for (i in 1:length(TFs)) {
  
  genes <- percentile %>% filter(TF == TFs[i]) %>% top_n(100, percentile_rank) %>% pull(target)
  
  Sobj_integrated <- AddModuleScore  (Sobj_integrated, 
                                      assay = "RNA", 
                                      list(genes), 
                                      nbin = 25, name = paste0(TFs[i], "_top100_"))
}






levels <- allCellTypes_deAnalysisDESeqWithBatchAsCoVar %>% 
  as.data.frame() %>% pull(CellType) %>% unique() %>% as.character()

Marker_list    <- vector("list", length = length(levels))

for (i in 1:length(levels)) {
  

df <- allCellTypes_deAnalysisDESeqWithBatchAsCoVar %>% 
  as.data.frame() %>% 
  filter(CellType == levels[i]) %>% 
  rownames_to_column("Gene") %>% 
  #mutate(Gene = str_sub(Gene, end=-3)) %>% 
  #filter(Gene %in% TFs) %>% 
  select(Gene, avg_logFC = "log2FoldChange", p_val_adj = "FDR")


Marker_list[[i]] <- 
EnhancedVolcano(df, df$Gene, 
                x = "avg_logFC",
                y = "p_val_adj", transcriptLabSize = 3,
                #xlim = c(-0.85, 0.85),
                FCcutoff = 0.25, title = paste("DE genes", levels[i]))


}



EnhancedVolcano(df, df$Gene, 
                x = "avg_logFC",
                y = "p_val_adj", 
                labSize = 3,
                xlim = c(-0.9, 0.9),
                FCcutoff = 0, 
                title = paste("HALLMARK: Fatty Acid Metabolism"))











