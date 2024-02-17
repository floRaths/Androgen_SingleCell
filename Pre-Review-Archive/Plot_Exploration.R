library(Seurat)
library(ggsci)
library(scales)
library(tidyverse)
library(ggpubr)
library(EnhancedVolcano)
library(pheatmap)


setwd("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/")
col <- pal_d3("category20")(20)

# load objects
Sobj_integrated <- readRDS("Seurat_Objects/Global_DataSets/Nuclei_All_Samples_Batch_Reduced.rds")
sobj <- Sobj_integrated

saveRDS(Sobj_integrated, "Seurat_Objects/DataSubsets/Epithelial_Cells_Reduced.rds")

DefaultAssay(Sobj_integrated) <- "RNA"



# function to save plot named "image"
save <- function(name, scale) {
  ggsave(plot   = image,
         path     = file.path("Figures/4.15_NicheNet_SCENIC_Refined/SVG/"),
         filename = paste0(name,".svg"),
         width    = 16, 
         height   = 9,
         scale    = scale,
         dpi      = 300)
  
  ggsave(plot   = image,
         path     = file.path("Figures/4.15_NicheNet_SCENIC_Refined/PNG/"),
         filename = paste0(name,".png"),
         width    = 16, 
         height   = 9,
         scale    = scale,
         dpi      = 300)
}



MetaData <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Output/MetaData_After_Trimming.rds")
Sobj_integrated <- AddMetaData(Sobj_integrated, dplyr::select(MetaData, Cluster = "Cluster", Group = "Group", CellType = "CellType", CellType_Plus = "CellType_Plus"))
rm(MetaData)

Idents(Sobj_integrated) <- "CellType"
levels(Sobj_integrated) <- c("LUM_HR-pos", "LUM_HR-neg", "Basal", "Endothelial", "Fibroblast", "Adipocyte", "Lymphoid", "Myeloid", "Neuronal", "Vasc_Acc")    
Sobj_integrated$CellType <- Idents(Sobj_integrated)


#c("nCount_RNA", "nFeature_RNA")
# good markers c("PTPRC", "ESR1", "VWF", "PLIN1", "ACTA2", "RELN", "CD163", "MS4A1")
# bad markers that used to be good markers  c("KRT5", "KRT8", "PECAM1", "CD8A", "ITGA6", "COL1A1", "MALAT1", "AGR2"), 

features = c("ESR1")
features = c("S.Score", "G2M.Score")
features = c("KRT5", "KRT8", "PECAM1", "CD8A", "ITGA6", "COL1A1", "MALAT1", "AGR2")

cells <- WhichCells(Sobj_integrated, idents = "Immune")



#image <- 
#plot <- 
  DimPlot(Sobj_integrated, 
        #cells = WhichCells(Sobj_integrated, cells = cells),
        reduction = "umap", 
        group.by  = "CellType_Plus",
        split.by  = "Type",
        pt.size = 0.2, 
        label = T, 
        repel = T, 
        label.size = 4,
        #ncol = 5, 
        #cols = col
        ) 



save("Strgrdhkdlomal_Cells", 1)  
  
  

  
select.cells <- CellSelector(plot = plot)
Idents(Sobj_integrated, cells = select.cells) <- "Adipocyte"
#Discard <- subset(Sobj_integrated, idents = "Discard")
Sobj_integrated <- RenameIdents(Sobj_integrated, 
                                "9" = "LUM_HR-pos_4"
                                )







  #image <- 
  VlnPlot(Sobj_integrated,  #idents = c("LUM_HR-neg"  ,     "Basal", "Endothelial",  "Fibroblast" ,  "Adipocyte"  ,  "Lymphoid" ,    "Myeloid"  ,  "Neuronal" ,   "Vasc_Acc" ),
        features = "UBC", 
        group.by = "CellType_Plus", 
        split.by = "Type",
        pt.size = 0, 
        ncol = 1)


  
#image <- 
    FeaturePlot(Sobj_integrated,
                #cells = WhichCells(Sobj_integrated, idents = "LUM_HR-pos", expression = AR_res1 > 0.5),
                features = c("PDE4B"), 
                split.by = "Type",
                pt.size = 1, 
                order = T, 
                cols = viridis::magma(n = 100), 
                #ncol = 3
                ) 

    
    
 save("Feature_ESR.ERBB4_noLUMp_split", 1)
  






left_join(rownames_to_column(AR_Marks, "Gene"), Human_Secretome, by = "Gene") %>% 
  left_join(select(Sign, Gene, Family), by = "Gene") %>% arrange(Category)

















# Cluster Proportions Boxplot ---------------------------------------------

levels(Sobj_integrated) <- c("LUM_HR-pos", "LUM_HR-neg", "Basal", "Endothelial", "Fibroblast", "Adipocyte", "Lymphoid", "Myeloid", "Vasc_Acc", "Neuronal")


Prop <- 
  as.data.frame(prop.table(table(Idents(Sobj_integrated), 
                                 Sobj_integrated$Sample), 
                           margin = 2)) %>% 
  dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq") %>% 
  mutate(Type = substring(.$Sample, 1, 2))



#image <- 
ggplot(Prop, 
       aes(Cluster, Freq, fill = Type)) + 
  geom_boxplot(width = 0.65) + 
  #geom_point(aes(color = Type)) + 
  geom_jitter(aes(color = Type), width = 0.1) +
  scale_y_sqrt()+ 
  stat_compare_means(paired = F, method = "wilcox")



save("Stromal_Proportions_Boxplot", 1)
