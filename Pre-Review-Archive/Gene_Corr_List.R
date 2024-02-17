
list.files("transgender_breast/R_Output/First_Parameters/Results/")
samples <- list.files("transgender_breast/R_Output/First_Parameters/Results/")

Sobj_list    <- vector("list", length = length(samples))                

### read in correaltion output
for (i in seq_along(samples)) {
  Sobj_list[[i]] <- read_csv(paste0("transgender_breast/R_Output/First_Parameters/Results/", samples[i], "/Corr.Results_", samples[i], ".csv"))
}

names(Sobj_list) <- samples # assign the correct names to the Sobj_List levels
Corr_Bind <- do.call(rbind.data.frame, Sobj_list) %>% 
  mutate(Cluster = paste(CellType, Corr.Cluster, sep = '_')) %>% 
  group_by(Cluster) %>% 
  mutate(Mean.Expr = mean(Log_FC)) 

rm(Sobj_list)


# count occurences among datasets and add it to df
genes <- 
  Corr_Bind %>% group_by(Gene) %>% summarise(events = n()) %>% arrange(desc(events)) 
Corr_Bind <- inner_join(ungroup(Corr_Bind), genes) %>% arrange(desc(events))

rm(genes)


# some early attempt to give names to clusters
Corr_Bind <- Corr_Bind %>% 
  ungroup() %>% 
  mutate(Clust.Name = Cluster) %>% 
  mutate(Clust.Name = recode(.$Clust.Name, 
                             "Adipocyte_2"     = "HER2_sig_Adp",
                             "Basal_1"         = "HER4_sig_Bsl",
                             "Basal_3"         = "Myogenesis_Bsl",
                             "Basal_7"         = "Metal_Ion_Binding_Bsl",
                             "Basal_9"         = "Heat_Stress_Bsl",
                             "Endothelial_7"   = "ANKRD_End",
                             "Fibroblast_6"    = "HER2+4_sig_Fbr",
                             "LUM_HR-neg_8"    = "HER4+RUNX2_sig_LUn",
                             "Myeloid_4"       = "HER4+Neuronal_sig_Mye",
                             "Vasc_Acc_6"      = "ANKRD_Vac",
                             "LUM_HR-pos_8"    = "AR_LUp",
                             "Fibroblast_5"    = "Heat_Stress_Fbr",
                             "Fibroblast_1"    = "ECM_interaction_Fbr",
                             "Fibroblast_3"    = "O-Glycosylation_Fbr",
                             "LUM_HR-pos_10"   = "Heat_Stress_LUp",
                             "LUM_HR-pos_7"    = "ANKRDB_LUp",
                             "LUM_HR-pos_5"    = "ESR-EGF_LUp",
                             "LUM_HR-pos_12"   = "NCAM_sig_LUp",
                             "Adipocyte_1"     = "LYPD_Adp",
                             "LUM_HR-neg_1"    = "Translation_LUn",
                             "Endothelial_4"   = "ESR-EGF_End",
                             "Neuronal_1"      = "ESR-EGF_Nrl",
                             "LUM_HR-neg_4"    = "ESR-EGF_LUn",
                             "LUM_HR-neg_7"    = "IP_Metabolism_LUn",
                             "LUM_HR-pos_9"    = "IP_Metabolism_LUp",
                             "LUM_HR-neg_6"    = "HER2+4_sig_LUn",
                             "LUM_HR-neg_2"    = "Muscle_contraction_LUn",
                             "LUM_HR-neg_9"    = "HER2_sig_LUn",
                             "LUM_HR-neg_10"   = "Heat_Stress_LUp",
                             "LUM_HR-pos_4"    = "KLKs_LUp",
                             "LUM_HR-pos_1"    = "Drug_Metabol_LUp",
                             "LUM_HR-pos_2"    = "O-Glycan_Synth_LUp",
                             "Myeloid_1"       = "Antigen-Presentation_&_(PD1)_Mye",
                             "Endothelial_5"   = "Interferon_sig_&_(PD1)_End",
                             "Myeloid_2"       = "Phagocytosis_Mye",
                             "Neuronal_4"      = "ECM_interaction_Nrl",
                             "Vasc_Acc_1"      = "Muscle_Contraction_Vac",
                             "Vasc_Acc_2"      = "ECM_interaction_Vac"
                             ))




Corr_Bind %>% 
  mutate(Regulon    = Clust.Name) %>% 
  mutate(Clust.Name = recode(.$Clust.Name, 
                             "Adipocyte_2"     = "HER2",
                             "Basal_1"         = "HER4/ESR",
                             "Basal_3"         = "Myogenesis_Bsl",
                             "Basal_7"         = "Metal_Ion_Binding_Bsl",
                             "Basal_9"         = "Heat_Stress_Bsl",
                             "Endothelial_7"   = "HER4/ESR",
                             "Fibroblast_6"    = "HER4/ESR",
                             "LUM_HR-neg_8"    = "HER4/ESR",
                             "Myeloid_4"       = "HER4/ESR",
                             "Vasc_Acc_6"      = "HER4/ESR",
                             "LUM_HR-pos_8"    = "AR_LUp",
                             "Fibroblast_5"    = "Heat_Stress_Fbr",
                             "Fibroblast_1"    = "ECM_interaction_Fbr",
                             "Fibroblast_3"    = "O-Glycosylation_Fbr",
                             "LUM_HR-pos_10"   = "Heat_Stress_LUp",
                             "LUM_HR-pos_7"    = "HER4/ESR",
                             "LUM_HR-pos_5"    = "HER4/ESR",
                             "LUM_HR-pos_12"   = "NCAM_sig_LUp",
                             "Adipocyte_1"     = "LYPD_Adp",
                             "LUM_HR-neg_1"    = "Translation_LUn",
                             "Endothelial_4"   = "HER4/ESR",
                             "Neuronal_1"      = "HER4/ESR",
                             "LUM_HR-neg_4"    = "HER4/ESR",
                             "LUM_HR-neg_7"    = "IP_Metabolism_LUn",
                             "LUM_HR-pos_9"    = "IP_Metabolism_LUp",
                             "LUM_HR-neg_6"    = "HER4/ESR",
                             "LUM_HR-neg_2"    = "Muscle_contraction_LUn",
                             "LUM_HR-neg_9"    = "HER2",
                             "LUM_HR-neg_10"   = "Heat_Stress_LUp",
                             "LUM_HR-pos_4"    = "AR_LUp",
                             "LUM_HR-pos_1"    = "Drug_Metabol_LUp",
                             "LUM_HR-pos_2"    = "O-Glycan_Synth_LUp",
                             "Myeloid_1"       = "Antigen-Presentation_&_(PD1)_Mye",
                             "Endothelial_5"   = "Interferon_sig_&_(PD1)_End",
                             "Myeloid_2"       = "Phagocytosis_Mye",
                             "Neuronal_4"      = "ECM_interaction_Nrl",
                             "Vasc_Acc_1"      = "Muscle_Contraction_Vac",
                             "Vasc_Acc_2"      = "ECM_interaction_Vac"
  ))








saveRDS(Corr_Bind, "Output/Correlation_Annotated.rds")
Corr_Bind <- readRDS( "Output/Correlation_Annotated.rds")


c("Myeloid_4", "Fibroblast_6", "Basal_1", "LUM_HR-neg_4", "Adipocyte_1", "Neuronal_1", "Endothelial_4")



filter(Corr_Bind, Clust.Name == "AR_Res?") %>% 
  arrange(Gene) %>% 
  group_by(Gene) %>% print(n = 100)
  

reactomeBar ("Lymphoid_3", filter(Corr_Bind, Cluster %in% c("Fibroblast_6"))$Gene, 15)


Clusters <- unique(Corr_Bind$Clust.Name)
dir.create(paste0("transgender_breast/R_Output/Pathways/"), recursive = T)



for (i in 1:length(Clusters)) {
  
  skip_to_next <- FALSE
  path =     paste0("transgender_breast/R_Output/Pathways/") # this is used for saving the plots
  
  tryCatch(image <- reactomeBar (name = Clusters[i], x = filter(Corr_Bind,   Clust.Name == Clusters[i])$Gene, y = 15), error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
  
  ggsave(plot   = image,
         path     = file.path(paste0(path)),
         filename = paste0(Clusters[i],".png"),
         width    = 16, 
         height   = 9,
         scale    = 1,
         dpi      = 300)
}




filter(Corr_Bind,   str_detect(Clust.Name, "^Fibroblast_9"))
filter(Corr_Bind,   CellType == "LUM_HR-pos")

filter(Corr_Bind, str_detect(Gene, "^RUN"))







Mat <- filter(Corr_Bind, events > 2) %>% 
  arrange(Gene) %>% 
  dplyr::select(1, 2, 8) %>% 
  pivot_wider(names_from = Clust.Name, values_from = Log_FC) %>%
  column_to_rownames("Gene") #%>% as.matrix()


left_join(filter(Corr_Bind, events > 2), Marker_Bind, by = c("Gene", "CellType")) %>% dplyr::select(1, 8, 10) %>% 
  pivot_wider(names_from = Clust.Name, values_from = avg_logFC) %>%
  column_to_rownames("Gene")




y <- pivot_wider(select(filter(Corr_Bind, events > 0), Gene, Log_FC, Cluster), names_from = Cluster, values_from = Log_FC)
#y <- left_join(select(filter(Corr_Bind, events > 0), Gene, CellType, Cluster), select(Marker_Bind, Gene, CellType, avg_logFC), by = c("Gene", "CellType"))
#y <- pivot_wider(select(y, Gene, Cluster, avg_logFC), names_from = Cluster, values_from = avg_logFC)

y[is.na(y)] <- 0

Mat <- column_to_rownames(y, "Gene")



pheatmap(Mat, cutree_cols = 1, cutree_rows = 1, fontsize = 7,
         cluster_rows = T, 
         cluster_cols = T, scale = "column", 
         color = viridis::viridis(n = 100),  
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation", border_color = 0)



x <- filter(Corr_Bind, events > 2) %>% 
  arrange(Log_FC) %>% 
  dplyr::select(1, 2, 8)

x <- as.data.frame(x)

ggplot(filter(Corr_Bind, events > 2), aes(Gene, Mean.Expr, color = Clust.Name)) + 
  geom_point()  + coord_flip()



