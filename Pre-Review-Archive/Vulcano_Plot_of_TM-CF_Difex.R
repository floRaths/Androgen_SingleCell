library(EnhancedVolcano)
library(tidyverse)
library(ggplot2)

setwd("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/")



# Reading response Output and choosing a Treatment Group (per patient and from combined DataSet) ---------

Marker_List <- readRDS("./Output/TMvsCF_Per_Cluster_Marker_List.rds")

names(Patient_A_MAST_Response)
treat  <- names(Patient_A_MAST_Response)[3] ### Choose the Treatment to plot here


### Combine the list of markers into a DataFrame (subset [[1, 2, 3, 4]] determines treatment)
Marker_Bind <- do.call(rbind.data.frame, Marker_List)

# First prepare and filter our data ----------------------

### Intersect the significant markers from A and B per cluster, and whether up or downregulated

UP   <- inner_join(filter(Pat_A_Mast_Bind, avg_logFC > 0.0, p_val_adj < 0.05),
                   filter(Pat_B_Mast_Bind, avg_logFC > 0.0, p_val_adj < 0.05),
                   by = c("gene", "Cluster"), suffix = c(".A", ".B")) %>%
  mutate(Avg_logFC = (avg_logFC.A+avg_logFC.B / 2))           %>%    # average the log_FC of A and B
  dplyr::select(gene, Cluster, avg_logFC.A, avg_logFC.B, Avg_logFC, p_val_adj.A, p_val_adj.B) %>% arrange(gene)  # selecting only the columns necessary for plotting
  
  
DOWN <- inner_join(filter(Pat_A_Mast_Bind, avg_logFC < -0.0, p_val_adj < 0.05),
                   filter(Pat_B_Mast_Bind, avg_logFC < -0.0, p_val_adj < 0.05),
                   by = c("gene", "Cluster"), suffix = c(".A", ".B")) %>%
  mutate(Avg_logFC = (avg_logFC.A+avg_logFC.B / 2))           %>%    # average the log_FC of A and B
  dplyr::select(gene, Cluster, avg_logFC.A, avg_logFC.B, Avg_logFC, p_val_adj.A, p_val_adj.B) %>% arrange(gene)  # selecting only the columns necessary for plotting


### making a DF of intersection
### This now contains all genes that were significantly UP or DOWN in both samples in the same cell types

UP   <- semi_join(filter(Combi_MAST_Bind, avg_logFC > 0, Cluster != "Low.UMI"), UP,   by = c("gene", "Cluster"))
DOWN <- semi_join(filter(Combi_MAST_Bind, avg_logFC < 0, Cluster != "Low.UMI"), DOWN, by = c("gene", "Cluster"))

TopTable   <- as_tibble(rbind (UP, DOWN))

### OR Use the DataFrames by themselves
TopTable   <- filter(Marker_Bind, p_val_adj < 0.05, 
                     #Cluster != "Low.UMI"
                     )


# Preparing to plot a combined Version -------------------------------------------------------

### assigning some key values to cluster names so that we can have custom shapes and colors per cluster
display.brewer.all()


TopTable <- filter(TopTable, Cluster == "Vascular_Accessory" | Cluster == "Neuronal_1"| Cluster == "Neuronal_2" | Cluster == "Dying?")

palette <- (pal_d3("category20")(4))
shapes  <- sample(c(15, 16, 17 , 18), length(unique(TopTable$Cluster)), replace = T)
names(palette) <- unique(TopTable$Cluster)
names(shapes)  <- unique(TopTable$Cluster)


TopTable <- mutate(TopTable, 
                   Color = recode(TopTable$Cluster, !!!palette), 
                   Shape = recode(TopTable$Cluster, !!!shapes))


key.shape        <- TopTable$Shape
key.color        <- TopTable$Color
names(key.shape) <- TopTable$Cluster
names(key.color) <- TopTable$Cluster



### Drawing the Plot
image <- 
  EnhancedVolcano(TopTable,
                lab = TopTable$gene,
                x = 'avg_logFC',
                y = 'p_val_adj',
                transcriptPointSize = 5, 
                #xlim        = c(-1.5, 0),
                #ylim        = c(0, 180),
                FCcutoff    = 0.25,
                pCutoff     = 0.05, 
                title       = "TM vs CF in Weird Cells", 
                subtitle    = "",
                shapeCustom = key.shape,
                colCustom   = key.color, 
                axisLabSize = 10,
                colAlpha    = 0.75)



save("Vulcano_WEIRD_Clusters", 1.5)

