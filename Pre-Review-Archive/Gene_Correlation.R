library(tidyverse)
library(readr)
library(pheatmap)
library(RColorBrewer)


trace("pheatmap", edit = T) # generate_annotation_col2

generate_annotation_colours2 = function(annotation, annotation_colors, drop){
  if(is.na2(annotation_colors)){
    annotation_colors = list()
  }
  count = 0
  for(i in 1:length(annotation)){
    annotation[[i]] = annotation[[i]][!is.na(annotation[[i]])]
    if(is.character(annotation[[i]]) | is.factor(annotation[[i]])){
      if (is.factor(annotation[[i]]) & !drop){
        count = count + length(levels(annotation[[i]]))
      }
      else{
        count = count + length(unique(annotation[[i]]))
      }
    }
  }
  
  factor_colors = dscale(factor(1:count), hue_pal(l = 75))
  #factor_colors = rev(colorRampPalette(brewer.pal(11,"PRGn"))(count))
  
  
  oldseed = NULL 
  if (exists(".Random.seed")) 
    oldseed = get(".Random.seed", pos=.GlobalEnv) 
  
  set.seed(3453)
  
  cont_counter = 2
  for(i in 1:length(annotation)){
    if(!(names(annotation)[i] %in% names(annotation_colors))){
      if(is.character(annotation[[i]]) | is.factor(annotation[[i]])){
        n = length(unique(annotation[[i]]))
        if (is.factor(annotation[[i]]) & !drop){
          n = length(levels(annotation[[i]]))
        }
        ind = sample(1:length(factor_colors), n)
        annotation_colors[[names(annotation)[i]]] = factor_colors[ind]
        l = levels(as.factor(annotation[[i]]))
        l = l[l %in% unique(annotation[[i]])]
        if (is.factor(annotation[[i]]) & !drop){
          l = levels(annotation[[i]])
        }
        
        names(annotation_colors[[names(annotation)[i]]]) = l
        factor_colors = factor_colors[-ind]
      }
      else{
        annotation_colors[[names(annotation)[i]]] = brewer_pal("div", cont_counter)(4)[1:4]
        cont_counter = cont_counter + 1
      }
    }
  }
  
  if(!is.null(oldseed)){ 
    assign(".Random.seed", oldseed, pos=.GlobalEnv) 
  } 
  else{ 
    remove(.Random.seed, pos=.GlobalEnv) 
  }
  
  return(annotation_colors)
}
is.na2 = function(x){
  if(is.list(x) | length(x) > 1){
    return(FALSE)
  }
  if(length(x) == 0){
    return(TRUE)
  }
  
  return(is.na(x))
}



corr.map <- function(name, fontsize, n.clust, scale) {
  
  # read the source data
  Genes   <-  read_csv                          (paste0("transgender_breast/source/", name, "/2_Genes_TM.csv"  ), col_names = "Genes")
  Matrix  <-  as.matrix(read_csv                (paste0("transgender_breast/source/", name, "/2_Matrix_TM.csv" ), col_names = Genes$Genes))
  Avrg    <-  column_to_rownames(cbind(read_csv (paste0("transgender_breast/source/", name, "/2_AvrgGen_TM.csv"), col_names = "Gene"),
              read_csv                          (paste0("transgender_breast/source/", name, "/2_AvrgVal_TM.csv"), col_names = "Log_FC")), "Gene")
  row.names(Matrix) <- Genes$Genes
  
  
  #Ann <- 
   # left_join(dplyr::select(Genes, Gene = "Genes"), Sign, by = "Gene") %>% 
  #  arrange(Family) %>% 
   # group_by(Gene) %>% 
  #  summarise(Family = paste(Family, collapse=", ")) %>% 
  #  inner_join(rownames_to_column(Avrg, "Gene"), by = "Gene") %>%
  #  mutate(Family = na_if(Family, "NA")) %>% 
   # left_join(Human_Secretome, by = "Gene") %>% 
  #  left_join(filter(Corr_Bind, CellType == name), by = "Gene") %>%
  #  column_to_rownames("Gene") %>% dplyr::select(Log_FC = "Log_FC.x", HR_rec = "Family", Secreted = "Category", Net.Annot = "Clust.Name", n.occur = "events")
  
  
  Ann <- 
    left_join(dplyr::select(Genes, Gene = "Genes"), Sign, by = "Gene") %>% 
    arrange(Family) %>% 
    group_by(Gene) %>% 
    summarise(Family = paste(Family, collapse=", ")) %>% 
    inner_join(rownames_to_column(Avrg, "Gene"), by = "Gene") %>%
    mutate(Family = na_if(Family, "NA")) %>% 
    left_join(Human_Secretome, by = "Gene") %>% 
    left_join(filter(Corr_Bind, CellType == name), by = "Gene") %>%left_join(TFs, by = "Gene") %>%
    column_to_rownames("Gene") %>% dplyr::select(Log_FC = "Log_FC.x", HR_rec = "Family", TF = "TF", Secreted = "Category", Net.Annot = "Clust.Name", n.occur = "events")
  
  

  # create output directories
  dir.create(paste0("transgender_breast/R_Output/Figures/", name, "/SVG/"), recursive = T)
  dir.create(paste0("transgender_breast/R_Output/Figures/", name, "/PNG/"), recursive = T)
  dir.create(paste0("transgender_breast/R_Output/Results/", name), recursive = T)
  
  path =     paste0("transgender_breast/R_Output/Figures/", name) # this is used for saving the plots
  
  # draw the heatmap
  image <- 
  pheatmap(Matrix, 
           cutree_rows  = n.clust,
           cluster_rows = T, 
           breaks = seq(from = -1, to = 1, by = 2/(length(Genes$Genes))),
           cluster_cols = T, 
           clustering_distance_rows = "euclidean", 
           clustering_distance_cols = "euclidean", 
           color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(length(Genes$Genes))), 
           border_color = NA, 
           treeheight_row = 0, 
           treeheight_col = 0, 
           annotation_row = Ann, 
           annotation_legend = T, 
           main = paste(name), 
           fontsize = 12,
           fontsize_col = fontsize,
           fontsize_row = fontsize, 
           #annotation_colors = mycolors
  )
  
  
  # save the plot
  ggsave(plot     = image,
         path     = file.path(paste0(path, "/SVG")),
         filename = paste0(name,".svg"),
         width    = 16, 
         height   = 9,
         scale    = scale,
         dpi      = 300)
  
  ggsave(plot   = image,
         path     = file.path(paste0(path, "/PNG")),
         filename = paste0(name,".png"),
         width    = 16, 
         height   = 9,
         scale    = scale,
         dpi      = 300)
  
  
  # extract correaltion clusters from heatmap
  A <- inner_join(rownames_to_column(Avrg, "Gene"),  # takes gene avrg_logFC and joins it with dendrogram ID
                  rownames_to_column(as.data.frame(cutree(image$tree_row, k = n.clust)), "Gene"))
  names(A)[3] <- "Corr.Cluster" # rename column name
  A <- arrange(A, Corr.Cluster) %>% mutate(CellType = name)
  A
  
  # save the output
  write_csv(A, paste0("transgender_breast/R_Output/Results/", name, "/Corr.Results_", name, ".csv"))
  
  }


samples <- unique(Corr_Bind$CellType)

for (i in seq_along(samples)) {
  corr.map(name = samples[i], fontsize = 7, n.clust = 1, scale = 1.7)
}





