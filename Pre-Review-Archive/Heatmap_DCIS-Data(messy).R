library(tidyverse)
library(pheatmap)

setwd("~/Box/Knott_Lab/Flo/Projects/Organoids/Analysis")


DefaultAssay(Hormone.Integrated) <- "RNA"

Pat_A <- AddModuleScore(Pat_A, 
                                    assay = "RNA", 
                                    list(Phgo$gene), 
                                    nbin = 5, name = "PHAGO", search = F)


FeaturePlot(Hormone.Integrated, features = c("PHAGO1"), 
            pt.size = 1, cols = inferno(10, begin = 0.2, end = 1), 
            reduction = "umap", min.cutoff = "q8", order = F, split.by = "Treat")



for (i in length(DCIS)) {
  
}


Idents(Pat_A) <- "Cluster"

VlnPlot(Pat_A, idents = "Basal_Myo", features = "PHAGO1", split.by = "Treat")


Pat_A <- subset(Hormone.Integrated, idents = "Patient_A")



wilcox.test(filter(Pat_A_Phago, Treat == "Control")$Phago.Score, 
            filter(Pat_A_Phago, Treat == "Est+Pro")$Phago.Score, 
            paired = F, alternative = "two.sided")

Label <- unique(Pat_A_Phago$Cluster)


plot_list <- vector("list", length = length(Label))

for (i in seq_along(Label)) {
  plot_list[[i]] <- ggplot(filter(Pat_A_Phago, Cluster == Label[i]), aes(Treat, Phago.Score, fill = Treat)) +
    geom_boxplot(show.legend = F) +
    ggtitle(Label[i]) 
  }

cowplot::plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[7]],plot_list[[5]], plot_list[[6]], plot_list[[8]], plot_list[[4]], ncol = 2)



for (i in seq_along(Label)) {
print(ggplot(filter(filter(Pat_A_Phago, Cluster == Label[i]), Treat == "Control" | Treat == "Est+Pro"), aes(Treat, Phago.Score)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox", label.x.npc = 0.4) +
  ggtitle(Label[i]))
}


### Creating a heatmap from responders

Set <- "IDC_"

DCIS_LFC_ER <- read_csv(file.path("utilities/100_for_Flo/",Set,"/data/LFC_ER.csv"  ))
DCIS_TPM <-    read_csv(file.path("utilities/100_for_Flo/",Set,"data/TPM.csv"     )) %>% mutate(Name = str_replace_all(Name, c("C10orf10" = "DEPP1", "SEPP1" = "SELENOP")))
DCIS_Meta <-   read_csv(file.path("utilities/100_for_Flo/",Set,"data/metadata.csv"))

annot <- DCIS_Meta %>% dplyr::select(Sample_ID, ERPR) %>% column_to_rownames("Sample_ID")

Resp <- filter(Combi_MAST_Bind, Cluster == "Basal.CC", avg_logFC > 0.3, p_val_adj < 0.05)

mat <- filter(DCIS_TPM, Name %in% Resp$gene) %>% dplyr::select(2:length(DCIS_TPM)) %>% column_to_rownames("Name") %>% as.matrix()

pheatmap(mat,
         scale = "row", 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "correlation",
         color = viridis(100),
         cluster_cols = T, 
         annotation = annot, 
         show_rownames = T, 
         show_colnames = FALSE)



for (i in test) {
  
}





TM1956_2 <- subset               (TM1956, subset = percent.mt < 3 & nCount_RNA < 20000)


p1 <- VlnPlot     (CF19301_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2 <- VlnPlot     (CF3920,    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p3 <- VlnPlot     (CF19301_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p4 <- VlnPlot     (TM1956_2,  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


CombinePlots(plots = list(p1, p2, p3, p4), ncol = 2)

































