


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




col <- pal_d3("category20")(20)
options(tibble.print_max = 150, tibble.print_min = 50)
`%notin%` <- Negate(`%in%`)


# Load DataSets -----------------------------------------------------------


# load Seurat Obejct
Sobj_integrated <- readRDS("Seurat_Objects/Reboot/Final_Integration/Sobj_GoodSamples_Integrated_Final_Aug05.rds")
Sobj_integrated <- readRDS("Seurat_Objects/Reboot/Harmony/CellType/Harmony_LUM_HR-pos.rds")





# Signaling Patway n = n_top connections
Sign <- sig.read(150)

# read in Human secretome dataset

TFs <- read_csv("~/Documents/SCENIC/pySCENIC_Full/scenicdata/Resources/hs_hgnc_tfs_ANK.txt", col_names = "Gene")

Human_Secretome <- read_excel("utilities/Human_Secretome_Annotated.xlsx")
Human_Secretome <-  dplyr::select(Human_Secretome, Gene = "gene", Category = "Category") #%>% mutate(Secretion = Category) %>% mutate(Secretion = !is.na(Secretion))

# read in TM-CF response marker list (on 10 clusters)
Marker_Bind <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Output/TM-CF_Response_CellType_MAST_list.rds")
Marker_Bind <- do.call(rbind.data.frame, Marker_Bind) %>% 
  dplyr::select(Gene = "gene", 2:6, CellType = "Cluster") %>% left_join(TFs, by = "Gene")

genes <- 
  Marker_Bind %>% group_by(Gene) %>% summarise(events = n()) %>% arrange(desc(events))  %>% print (n = 150) 
Marker_Bind <- inner_join(ungroup(Marker_Bind), genes) %>% arrange(desc(events))

# read in TM-CF response marker list (on 40 clusters)
#CTPlus_Bind <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/utilities/CellType_Plus_Bind.rds") %>%
#  filter(p_val_adj <= 0.05) %>% left_join(Human_Secretome, by = "Gene") %>% left_join(TFs, by = "Gene")

genes <- 
  CTPlus_Bind %>% group_by(Gene) %>% summarise(events = n()) %>% arrange(desc(events)) %>% print (n = 150) 
CTPlus_Bind <- inner_join(ungroup(CTPlus_Bind), genes) %>% arrange(desc(events))

rm(genes)


adjacencies <- read_tsv("~/Documents/SCENIC/pySCENIC_Full/scenicdata/8kVarGenes_ANK_Added/expr_mat.adjacencies.tsv")
auc_mtx     <- read_csv("~/Documents/SCENIC/pySCENIC_Full/scenicdata/4k_Var_Genes/Results/Epithelial_auc_mtx.csv")


Average_TF_Expression <- 
  readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Output/Average_TF_Expression_CT_All.rds") %>% 
  ungroup() %>%
  filter(expression == "expressed", perc >= 33) %>% 
  dplyr::select(CellType = 2, 3) %>% distinct() %>% 
  group_by(TF) %>% 
  mutate(events = n()) %>%
  left_join(dplyr::select(Marker_Bind, Gene, CellType, avg_logFC), by = c("TF" = "Gene", "CellType"))


######################################################################
#### Step 2: Read in NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks:

ligand_target_matrix = readRDS("NicheNet/utilities/ligand_target_matrix.rds")
ligand_target_matrix[1:5,1:5]

lr_network = readRDS("NicheNet/utilities/lr_network.rds")
head(lr_network)

weighted_networks = readRDS("NicheNet/utilities/weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

weighted_networks$lr_sig <- dplyr::select(weighted_networks$lr_sig, Ligand = "from", Receptor = "to", 3)

head(weighted_networks$lr_sig)
head(weighted_networks$gr)







receptor_modules <- c("~/Documents/NicheNet/receptor_genes/")

list.files(receptor_modules)
receptors <- substr(list.files(receptor_modules), 1,nchar(list.files(receptor_modules))-4)
receptor_list <- vector("list", length = length(receptors))  

#read_tsv(paste0("utilities/receptor_genes/AR.txt")

for (i in 1:length(list.files(receptor_modules))) {
  receptor_list[[i]] <- read_tsv(paste0(receptor_modules, 
                                        receptors[i], ".txt"), 
                                 col_names = "Gene") %>% mutate(Receptor = receptors[i])
}

names(receptor_list) <- receptors
Recp_Bind <- do.call(rbind.data.frame, receptor_list) %>% dplyr::select(Receptor, Gene)

rm(receptor_list, receptors, receptor_modules)





# Gene Correlation list (Bassem)
Corr_Bind   <- readRDS("Output/Correlation_Annotated.rds")

AR_Markers  <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/utilities/AR_Markers.rds")
AR_Markers  <- rownames_to_column(AR_Markers, "Gene") %>% 
  filter(p_val_adj <= 0.05) %>% 
  dplyr::select(1,3,6) %>% 
  left_join(TFs, by = "Gene") %>% 
  arrange(-avg_logFC)


# joining TM-CF response with secretome and receptor signaling
TopTable <- 
  full_join(Marker_Bind, Corr_Bind, by = c("Gene", "CellType")) %>% 
  left_join(Sign, by = "Gene") %>% 
  left_join(Human_Secretome, by = "Gene") %>%
  mutate(Secreted = Category) %>% 
  mutate(Secreted = !is.na(Secreted)) %>% 
  dplyr::select(              "Gene", 
                              DE.evts     = "events.x", 
                              NW.evts     = "events.y", 
                              DE_FC       = "avg_logFC", 
                              DE_pval     = "p_val_adj", 
                              Corr_FC     = "Log_FC",
                              NW_Avgr     = "Mean.Expr", 
                              CellType    = "CellType",
                              #CellType_NW = "CellType.y",
                              HR_Rcptr    = "Family", 
                              Netw.ID     = "Clust.Name",
                              "Secreted", 
                              "Category") %>% 
  as_tibble()



TFs <- read_tsv("~/Documents/SCENIC/pySCENIC_Full/scenicdata/Resources/hs_hgnc_tfs.txt", col_names = "Gene") %>% mutate(TF = "YES")

Diffex_TFs <- 
CTPlus_Bind %>% dplyr::select(1,3,7,8) %>% left_join(TFs, by = "Gene") %>% filter(TF == "YES") %>% pull(Gene) %>% unique()

# determine secreted genes
secr_up <- 
  left_join(AR_Markers, Human_Secretome, by = "Gene") %>% 
  filter(avg_logFC >=0.35, p_val_adj <=0.05) %>% arrange(Category) %>% head(20) #%>% pull(Gene)

secr_do <- 
left_join(AR_Markers, Human_Secretome, by = "Gene") %>% 
  filter(avg_logFC <=-0.35, p_val_adj <=0.05) %>% arrange(Category) %>% head(15)

LumAnn <- 
rbind(secr_up, secr_do) %>% dplyr::select(-3,-4) %>% column_to_rownames("Gene")





adjacencies <- 
  adjacencies %>% dplyr::select(1, Gene = "target", 3) %>% 
  left_join(mean_logFC) %>% 
  left_join(TFs, by = "Gene") %>% dplyr::select(TF = 1, target = 2, 3:6)



query <- "CUX2"

adjacencies %>% filter(TF   == query) %>% arrange(-importance)
adjacencies %>% filter(target == query) %>% arrange(-importance)

lr_network %>% filter(to      == query) 
lr_network %>% filter(from    == query)

weighted_networks$lr_sig %>% filter(Receptor == query) %>% arrange(-weight)
weighted_networks$lr_sig %>% filter(Ligand   == query) %>% arrange(-weight)

Recp_Bind %>% filter(Gene     == query)
Recp_Bind %>% filter(Receptor == query)

CTPlus_Bind %>% filter(Gene %in% query) %>% arrange(-avg_logFC)




CTPlus_Bind %>% 
  dplyr::select(1, events) %>% 
  arrange(-events) %>% 
  distinct() %>% 
  left_join(TFs) %>% 
  left_join(mean_logFC) %>% 
  as_tibble() %>% 
  filter(events >= 10) %>% #pull(Gene)
  left_join(top5_drivers, by = "Gene") %>% arrange(Driver) %>%
  print(n = 500)
  
  
top5_drivers <- 
adjacencies %>% filter(target %in% query) %>% 
  group_by(target) %>% 
  top_n(5, importance) %>% 
  dplyr::select(Driver = "TF", Gene = "target", 3)



###  average expression across subclusters
#mean_logFC <- 
CTPlus_Bind %>% 
  group_by(Gene) %>% 
  summarise(mean_logFC = mean(avg_logFC)) %>% 
  arrange(-mean_logFC) %>% 
  left_join(dplyr::select(CTPlus_Bind, Gene, events)) %>% 
  distinct() %>% arrange(-events)




CTPlus_Bind %>% filter(Gene %in% filter(Recp_Bind, Gene == query)$Receptor) %>% 
  dplyr::select(1, 3, 7, 8, 9) %>% as_tibble() %>% arrange(Gene) %>% print(n = 50) 


CTPlus_Bind %>% filter(Gene %in% top_n(filter(weighted_networks$lr_sig, Receptor == query), 10, weight)$Ligand) %>% 
  dplyr::select(1, 3, 7:9) %>% as_tibble() %>% arrange(-events)


TopTable %>% filter(Gene %in% top_n(filter(weighted_networks$lr_sig, Receptor == query), 20, weight)$Ligand) %>% 
  dplyr::select(1:12) %>% as_tibble() %>% arrange(Gene)




direct_tfs <- 
filter(weighted_networks$lr_sig, Ligand %in% secr$Gene) %>% 
  left_join(dplyr::select(TFs, Receptor = 1, 2)) %>% 
  filter(TF == "YES") %>% #pull(Receptor) %>% unique()
  arrange(-weight) %>% dplyr::select(1, Gene = "Receptor", 3, 4)



direct_tfs %>% 
  inner_join(CTPlus_Bind, by = "Gene") %>% 
  dplyr::select(-5, -7 ,-8) %>% 
  arrange(-events) %>% print(n = 200)


x %>% group_by(Gene) %>% 
  mutate(mean = mean(avg_logFC)) %>% 
  dplyr::select(1:3, 8:10) %>% 
  distinct() %>% 
  arrange(-events) %>% 
  print(n = 200)


CTPlus_Bind %>% filter(Gene %in% direct_tfs) %>% as_tibble()



receptors_oi <- CTPlus_Bind %>% filter(avg_logFC >= 0.5 | avg_logFC <= -0.05) %>% pull(Gene) %>% unique()


weighted_networks$lr_sig %>% 
  filter(Receptor %in% receptors_oi) %>% pull(Receptor) %>% unique() %>% as.data.frame() %>% write_tsv("~/Documents/NicheNet/Receptors_To_Generate.csv")




# Add Module Score --------------------------------------------------------


genes <- Corr_Bind %>% filter(CellType == "Endothelial", str_detect(Clust.Name, "^Endothelial_2")) %>% pull(Gene)


Sobj_integrated <- AddModuleScore  (Sobj_integrated, 
                                    assay = "RNA", 
                                    #list(filter(Corr_Bind, str_detect(Clust.Name, "^KLK") | str_detect(Clust.Name, "^AR"))$Gene), 
                                    list(test$name), 
                                    nbin = 25, name = "X1_")


#saveRDS(Sobj_integrated, "Seurat_Objects/DataSubsets/Epithelial_Cells_Reduced.rds")

DefaultAssay(Sobj_integrated) <- "RNA"



Immune_Meta <- 
Sobj_integrated@meta.data %>% 
  rownames_to_column("ID") %>% 
  select("ID", "Sample", "Tissue", "Type", "Group", "Group2", "CellType", "CellType_All", "Keep") %>%
  column_to_rownames("ID")











select.cells <- CellSelector(plot = plot)

Idents(Sobj_integrated) <- "Stroma"
Idents(Immu, cells = doubles) <- "Double"

Idents(Sobj_integrated, cells = cells) <- "3"
Sobj_integrated$CellType_All <- Idents(Sobj_integrated)
Sobj_integrated$CellType_All <- Sobj_integrated$Tissue

Sobj_integrated <- RenameIdents(Sobj_integrated, c("9" = "6",
                                                   "9" = "6",
                                                   "9" = "6",
                                                   "9" = "6",
                                                   "9" = "6",
                                                   "9" = "6",))



Sobj_integrated$Tissue <- Sobj_integrated$Sample
Sobj_integrated$Sample <- recode(Sobj_integrated$Sample, "CF-249-347_2" = "CF-249-347")
Sobj_integrated$Sample <- recode(Sobj_integrated$Sample, "TM-1956_2" = "TM-1956")


# DimPlot -----------------------------------------------------------------

#image <- 
#plot1 <- 
DimPlot(Sobj_integrated, 
        #cells = WhichCells(Sobj_integrated, idents = c("TM-3937", "CF-2797")),
        reduction = "umap", 
        #group.by  = "Res.0.2_Added",
        group.by  = "CellType",
        #cells = WhichCells(Sobj_integrated, downsample = 100), 
        #group.by  = "CellType",
        #split.by = "Type",
        pt.size =1, 
        label = T,
        repel = T, 
        label.size = 8,
        #split.by  = "Sample", ncol = 6, 
        #cols = viridis::magma(n = 2, begin = 0.3, end = 0.8),
        #cols = col
        ) 

plot1 + plot2 #+ patchwork::plot_layout(nrow = 2)


save2("Immune_Markers", 1.5)

levels(Sobj_integrated) <- c("LUM_HR-pos", "Fibroblast", "Adipocyte", "Blood_EC", "Lymphatic_EC", "Vasc_Accessory", "LUM_HR-neg", "Basal", "Myeloid", "Lymphoid")
# VlnPLot -----------------------------------------------------------------

query <- c("KLF6", "ZBTB7A")#, "NR4A1", "KLF6", "FOXO3")



names <- colnames(data)[22:244]


for (i in 1:length(names)) {
  
quick <- 
VlnPlot(Sobj_integrated,
        #idents = "Adipocyte",
        #features = paste0(query,"_1"),
        features = "ADCY2",
        #group.by = "Type", 
        split.by = "Type",
        pt.size = 0, 
        ncol = 1, 
        cols = viridis::magma(n = 2, begin = 0.3, end = 0.8)
        )

save2(paste("ADCY_", names[i]), 0.5)

}
  
# Featureplot -------------------------------------------------------------

#quick2 <- 
FeaturePlot(Sobj_integrated, reduction = "umap",
            #cells = WhichCells(Sobj_integrated, expression = SCORE_1 >= 0.25),
            #features = c("ANKRD30A", "VWF", "PGR", "TP63", "nCount_RNA", "KIT", "PTPRC", "COL6A3", "RPLP1"), #cells = cells,
            features = query, 
            split.by = "Type",
            pt.size = 0.5, 
            order = T, 
            #min.cutoff = 10000,
            cols = viridis::magma(n = 100),
            #cols = c("grey", "Green", "Red")
            #ncol = 3
            )


quick1 + quick2 + plot_layout(ncol = 1)



features3 <- 
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






# Cluster Proportions Boxplot ---------------------------------------------



levels(Sobj_integrated) <- table(Sobj_integrated$CellType) %>% as.data.frame() %>% arrange(-Freq) %>% pull(Var1)


Prop <- 
  as.data.frame(prop.table(table(Sobj_integrated$Group, 
                                 Sobj_integrated$Sample), 
                           margin = 2)) %>% 
  dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq") %>% 
  mutate(Type = substring(.$Sample, 1, 2)) %>% arrange(factor(Cluster, levels = levels), desc(Freq))



Test <- 
as.data.frame(table(Sobj_integrated$RNA_snn_res.0.2, 
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
  
  
Sobj_integrated <- AddMetaData(Sobj_integrated, Test)

  filter(Cluster == 4) %>% 
  arrange(-Freq)


#image <- 
ggplot(Prop, 
       aes(reorder(Cluster, -Freq), Freq, fill = Type)) + 
  geom_boxplot(width = 0.65, outlier.alpha = 0) + 
  #geom_point(aes(color = Type)) + 
  geom_jitter(aes(color = Type), width = 0.1) +
  scale_y_sqrt() + 
  stat_compare_means(paired = F, method = "wilcox") + 
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"), 
        axis.text.x = element_text(size=15))



save("Stromal_Proportions_Boxplot", 1)




ggplot(data = Prop, aes(x = "", y = Freq, fill = Cluster )) + 
  geom_bar(stat = "identity", position = position_fill()) +
  #geom_text(aes(label = Freq), position = position_fill(vjust = 0.5)) +
  coord_polar(theta = "y") +
  facet_wrap(~ Sample, nrow = 3)  #+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) + 
  theme(legend.position='bottom') + 
  guides(fill=guide_legend(nrow=2, byrow=TRUE))



# Vulcano Plot ------------------------------------------------------------

df <- MB
df <- MB %>% filter(str_detect(CellType, "^Adipo"))
df <- mutate(df, Color = recode(Cluster, !!!palette))
key.color        <- df$Color
names(key.color) <- df$Secreted

EnhancedVolcano(filter(CTP_Secreted, Cluster == "Adipocyte"), filter(CTP_Secreted, Cluster == "Adipocyte")$Gene, 
                x = "avg_logFC",
                y = "p_val_adj", 
                #xlim = c(-1.50, -0.5),
                FCcutoff = 0.25, title = "LUM_HR-pos_DE Genes"
                )


df <- Treatment_Response %>% filter(Cluster == "LUM_HR-pos", Gene %in% AR_PGR)

EnhancedVolcano(df, df$Gene, 
                x = "avg_logFC",
                y = "p_val_adj", transcriptLabSize = 3,
                #xlim = c(-0.85, 0.85),
                FCcutoff = 0.25, title = "TFs_with_AR_motif_in_Adipocytes"
)





celltype = "LUM_HR-neg"

adjacencies <- read_tsv(paste0("Output/Reboot/SCENIC/expr_mat.adjacencies.tsv"))

percentile <- 
  adjacencies %>% 
  group_by(TF) %>% 
  mutate(percentile_rank = ntile(importance, 100))



tf = "ANKRD30A"

targets <- percentile %>% ungroup() %>% filter(TF %in% tf, percentile_rank > 90) %>% pull(target)
df <- Marker_Bind %>% ungroup() %>% filter(CellType == celltype, Gene %in% targets)

df <- percentile    %>% left_join(x, by = "target") %>% filter(TF == "ANKRD30A", percentile_rank > 99)
df <- percentile_CF %>% left_join(x, by = "target") %>% filter(TF == "ANKRD30A", percentile_rank > 99)






levels <- levels(Sobj_integrated)

for (i in 1:length(levels(Sobj_integrated))) {

celltype <- levels[i]
  
df <- Marker_Bind %>% filter(CellType == levels[i])

image <- 
EnhancedVolcano(df, df$Gene, 
                x = "avg_logFC",
                y = "p_val_adj", 
                #colCustom   = key.color,
                #xlim = c(-1.50, -0.5),
                FCcutoff = 0.25, 
                pCutoff = 0.05,
                title = paste("DE-genes in", levels[i]), 
                colAlpha = 0.8
                )

save(name = paste0("Vulcano_DE_", levels[i]), scale = 1.2)

}



display.brewer.all()

palette <- brewer.pal(n = length(unique(df$Cluster)), name = 'Set2')
shapes  <- sample(c(15, 16, 17 , 18), length(unique(df$Cluster)), replace = T)
names(palette) <- unique(df$Cluster)
names(shapes)  <- unique(df$Cluster)


palette <- brewer.pal(n = 2, name = 'Set2')
names(palette) <- c("YES", "NO")

df <- mutate(df, Color = recode(df$Secreted, !!!palette))


df <- mutate(df, 
                   Color = recode(df$TF, !!!palette), 
                   Shape = recode(df$Secreted, !!!shapes))


key.shape        <- df$Shape
key.color        <- df$Color
names(key.shape) <- df$Cluster
names(key.color) <- df$Secreted


### Drawing the Plot
EnhancedVolcano(df,
                lab = df$Gene,
                x = 'avg_logFC',
                y = 'p_val_adj',
                pointSize = 5, 
                #xlim        = c(0.2, 2),
                #ylim        = c(0, 180),
                FCcutoff    = 0.25,
                pCutoff     = 0.05, 
                title       = paste0("LUM-Pos DE Genes with CUX2 binding sites in mouse liver"), 
                #subtitle    = "Both Patients Analyzed Combined",
                #shapeCustom = key.shape,
                #colCustom   = key.color, 
                axisLabSize = 10,
                colAlpha    = 0.75)





# AR-Response Scoring -----------------------------------------------------

Sobj <- AddMetaData(Sobj, metadata = column_to_rownames(dplyr::select(X2, 1, CUX2_4k), "Cell"))

cutoff <- 0.35

FeaturePlot(LUM, cells = WhichCells(LUM, idents = "LUM_HR-pos", expression = CUX2_4k >= cutoff),
            reduction = "SCENIC", 
            features = "CUX2_4k", 
            order = T,
            cols = viridis::magma(n = 100), 
            pt.size = 2)



FeaturePlot(Sobj, cells = WhichCells(Sobj, idents = c("LUM_HR-pos", "CUX2")), #expression = `KLK3` > 0),
            reduction = "umap", 
            features = "KLK3", order = T,
            cols = viridis::magma(n = 100), 
            pt.size = 2)


Idents(LUM) <- "CellType"


Idents(LUM, cells = WhichCells(LUM, idents = "LUM_HR-pos", expression = CUX2_4k >= cutoff)) <- "CUX2"



CUX2_Marks <- FindMarkers(LUM, 
                        ident.1 = "CUX2", 
                        ident.2 = "LUM_HR-pos",
                        assay   = "RNA", 
                        slot    = "data", 
                        test.use = "MAST"
                        ) %>% rownames_to_column("Gene")

CUX2_Marks %>% filter(p_val_adj <= 0.05, str_detect(Gene, "^KLK")) %>% arrange(-avg_logFC)


CUX2_Marks <- 
CUX2_Marks %>% left_join(Human_Secretome, by = "Gene") %>% left_join(TFs, by = "Gene")




# GRNboost TF Modules -----------------------------------------------------



#Corr_TFs    <- intersect(filter(CTPlus_Bind, Cluster == "Adipocyte")$Gene, adjacencies$TF)
#Corr_target <- intersect(filter(CTPlus_Bind, Cluster == "Adipocyte")$Gene, adjacencies$target)
Corr_TFs    <- intersect(filter(CTPlus_Bind, str_detect(Cluster, "^Basal"))$Gene, adjacencies$TF)
Corr_target <- intersect(filter(CTPlus_Bind, str_detect(Cluster, "^Basal"))$Gene, adjacencies$target)


Ann <- 
  CTPlus_Bind %>% filter(Cluster == "Adipocyte") %>% 
  dplyr::select(1,3) %>% column_to_rownames("Gene")

Ann <- 
  CTPlus_Bind %>% 
    filter(str_detect(Cluster, "^Basal")) %>% 
    dplyr::select(1, 3, 7) %>% 
    pivot_wider(names_from = Cluster, values_from = avg_logFC) %>% 
    column_to_rownames("Gene")


mat <- filter(adjacencies, TF %in% Corr_TFs & target %in% Corr_target) %>% 
  pivot_wider(names_from = target, values_from = importance) %>%
  column_to_rownames("TF")


mat[is.na(mat)] <- 0


image <- 
  pheatmap(mat, scale = "column", fontsize = 10, 
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "euclidean",
           color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)),
           annotation_col = Ann,
           #annotation_colors = ann_colors, 
           treeheight_row = 0, treeheight_col = 0, 
           cutree_rows = 1, 
           cutree_cols = 10,
           fontsize_row = 9,
           fontsize_col = 6, 
           labels_col = NULL
           )




modules <- 
  rownames_to_column(as.data.frame(cutree(image$tree_col, k = 30))) %>% 
  dplyr::select(Gene = "rowname", module = 2)


for (i in 1:length(unique(modules$module))) {
#for (i in 2:2) {
    
  plot1 <- 
    pheatmap(t(mat[,filter(modules, module == i)$Gene]), 
             scale = "row", 
             fontsize = 10, 
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation",
             color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)),
             annotation_row = Ann,
             #annotation_col = Ann,
             #annotation_colors = ann_colors, 
             treeheight_row = 0, 
             treeheight_col = 0, 
             cutree_rows = 1, 
             angle_col = 90, silent = T,
             cutree_cols = 1,
             fontsize_row = 8,
             fontsize_col = 8, 
             labels_col = NULL, 
             border_color = "lightgrey"
    )
  
  
  plot2 <- reactomeBar(paste(i), filter(modules, module == i)$Gene, 10)
  

  
  plot <- plot_grid(as.grob(plot1), plot2, nrow = 2, rel_heights = c(2,1))
  print(plot)
  
}


quick <- 
  VlnPlot(Sobj_integrated,
          #idents = "Adipocyte",
          #features = paste0(query,"_1"),
          features = "Mod_19_1", 
          group.by = "CellType", 
          split.by = "Type",
          pt.size = 0, 
          ncol = 1
  )


save2("Vln_Example", 0.5)


reactomeBar("whateva", top_n(filter(adjacencies, TF == "ANKRD30A"), 30, importance)$target, 10)










### Slingshot


rd <- Sobj_integrated@reductions$umap@cell.embeddings

rd <- Sobj_integrated@reductions$harmony@cell.embeddings[,1:10]


cells <- TM


lin1 <- getLineages(rd[cells,], Idents(Sobj_integrated)[cells], start.clus = 'Progenitors')
lin1

plot(rd[cells,], col = col[Idents(Sobj_integrated)[cells]], asp = 1, pch = 16)
lines(lin1, lwd = 3, col = 'black')







levels <- c("LUM_HR-pos", "LUM_HR-neg", "Basal", "Fibroblast", "Adipocyte", "Blood_EC", "Lymphatic_EC", "Vasc_Accessory", "Myeloid", "Lymphoid")

cells = WhichCells(Sobj_integrated, downsample = 500)
Marks <- read_excel("Output/Reboot/Selected_HeatMap_Markers.xlsx") %>% pivot_longer(cols = 1:10) %>% arrange(factor(name, levels = levels), desc(value))

feats = Marks %>% pull(value)

cells <- Sobj_integrated@meta.data[cells, ] %>% dplyr::select("CellType") %>% arrange(factor(CellType, levels = levels)) %>% rownames()



mat <- as.matrix(Sobj_integrated@assays$RNA@data[feats, cells]) %>% scale(scale = T,center = T)


pheatmap((mat), 
         cluster_cols = T, 
         cluster_rows = T, 
         labels_col = F, #breaks = myBreaks,
         scale = "row", 
         color = viridis::magma(n = 100))






cells <- WhichCells(Sobj_integrated, downsample = 100)
levels <- levels(Sobj_integrated)

Marks <- read_excel("Output/Reboot/Selected_HeatMap_Markers.xlsx", sheet = 2) %>% pivot_longer(cols = 1:9) %>% filter(!is.na(value)) %>% arrange(factor(name, levels = levels), desc(value))

feats = Marks %>% pull(value)

cells <- Sobj_integrated@meta.data[cells, ] %>% dplyr::select("CNumber") %>% arrange(factor(CNumber, levels = levels)) %>% rownames()

DoHeatmap(Sobj_integrated, cells = cells, features = feats, group.colors = col)






genes <- c("ABHD5", "OSBPL8", "PPARG", "FITM1", "OSBPL11", "LPL", "CIDEA", "PNPLA2", "PPARA", "FITM2", "PLIN5", "APOC4")
