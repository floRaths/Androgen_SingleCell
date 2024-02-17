library(EnhancedVolcano)
library(tidyverse)
library(ggplot2)

setwd("~/Box/Knott_Lab/Flo/Projects/Organoids/Analysis/")



# Reading response Output and choosing a Treatment Group (per patient and from combined DataSet) ---------

Patient_A_MAST_Response <- readRDS("./2019.02.17_Oganoid_Hormone_TimeCourse/Outputs/V3/Markers/Response.Markers/MAST_Response/Patient_A_MAST_Response.rds")
Patient_B_MAST_Response <- readRDS("./2019.02.17_Oganoid_Hormone_TimeCourse/Outputs/V3/Markers/Response.Markers/MAST_Response/Patient_B_MAST_Response.rds")
Combined_MAST_Response  <- readRDS("./2019.02.17_Oganoid_Hormone_TimeCourse/Outputs/V3/Markers/Response.Markers/MAST_Response/Combined_MAST_Response.rds")


names(Patient_A_MAST_Response)
treat  <- names(Patient_A_MAST_Response)[3] ### Choose the Treatment to plot here


### Combine the list of markers into a DataFrame (subset [[1, 2, 3, 4]] determines treatment)
Pat_A_Mast_Bind <- do.call(rbind.data.frame, Patient_A_MAST_Response[[treat]])
Pat_B_Mast_Bind <- do.call(rbind.data.frame, Patient_B_MAST_Response[[treat]])
Combi_MAST_Bind <- do.call(rbind.data.frame, Combined_MAST_Response [[treat]])

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
TopTable   <- filter(Combi_MAST_Bind, p_val_adj < 0.05, Cluster != "Low.UMI")
TopTable_A <- filter(Pat_A_Mast_Bind, p_val_adj < 0.05, Cluster != "Low.UMI")
TopTable_B <- filter(Pat_B_Mast_Bind, p_val_adj < 0.05, Cluster != "Low.UMI")


# Preparing to plot a combined Version -------------------------------------------------------

### assigning some key values to cluster names so that we can have custom shapes and colors per cluster
display.brewer.all()

palette <- brewer.pal(n = length(unique(TopTable$Cluster)), name = 'Set2')
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
EnhancedVolcano(TopTable,
                lab = TopTable$gene,
                x = 'avg_logFC',
                y = 'p_val_adj',
                transcriptPointSize = 5, 
                #xlim        = c(0.2, 2),
                #ylim        = c(0, 180),
                FCcutoff    = 0.25,
                pCutoff     = 0.05, 
                title       = paste0("Response to ", treat), 
                subtitle    = "Both Patients Analyzed Combined",
                shapeCustom = key.shape,
                colCustom   = key.color, 
                axisLabSize = 10,
                colAlpha    = 0.75)



#This actually save the plot in a image
ggsave(plot     = image,
       path     = file.path("./Figrues/10.24_Nuclei-Alignment_MAST_Responses/"),
       filename = paste0(image$labels$title, "_", image$labels$subtitle, ".svg"),
       width    = 16, 
       height   = 9,
       scale    = 1.5,
       dpi      = "screen")



# Using Combined Markers to draw the Plot ---------------------------------

TopTable   <- filter(TopTable, Cluster != "Low.UMI")

display.brewer.all()

palette <- brewer.pal(n = length(unique(TopTable$Cluster)), name = 'Set2')
shapes  <- sample(c(15, 16, 17 , 18), length(unique(TopTable$Cluster)), replace = T)
names(palette) <- unique(TopTable$Cluster)
names(shapes)  <- unique(TopTable$Cluster)

TopTable <- mutate(TopTable, Color = recode(TopTable$Cluster, !!!palette), Shape = recode(TopTable$Cluster, !!!shapes))


key.shape  <- TopTable$Shape
key.color  <- TopTable$Color
names(key.shape) <- TopTable$Cluster
names(key.color) <- TopTable$Cluster




### Draw Plots A and B

Ylim = 200
Xlim = 0

A <- EnhancedVolcano(TopTable,
                  lab = TopTable$gene,
                  x = 'avg_logFC.A',
                  y = 'p_val_adj.A',
                  #xlim     = c(-2, Xlim),
                  #ylim     = c(0, Ylim),
                  FCcutoff = 0.25,
                  pCutoff  = 0.05, 
                  transcriptPointSize = 4, 
                  title    = paste0("Response to ", treat), 
                  subtitle = "Union Response - Patient A",
                  shapeCustom = key.shape,
                  colCustom   = key.color,
                  axisLabSize = 10,
                  colAlpha = 1)


B <- EnhancedVolcano(TopTable,
                  lab = TopTable$gene,
                  x = 'avg_logFC.B',
                  y = 'p_val_adj.B',
                  #xlim     = c(-2, Xlim),
                  #ylim     = c(0, Ylim),
                  FCcutoff = 0.25,
                  pCutoff  = 0.05,
                  legendVisible = F,
                  transcriptPointSize = 4, 
                  title    = paste0("Response to ", treat), 
                  subtitle = "Union Response - Patient B",
                  shapeCustom = key.shape, 
                  colCustom   = key.color,
                  axisLabSize = 10,
                  colAlpha = 1)


cowplot::plot_grid(A, B)
ggsave(file = "./Figrues/10.24_Nuclei-Alignment_MAST_Responses/test.eps", plot = image, width=16, height=9, )

















