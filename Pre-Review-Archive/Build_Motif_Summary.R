
Motifs_By_TF <- readRDS("ATAC_Integration/Motifs_by_TF_database.rds")


avg_merge <- cluster.averages %>% as_tibble() %>% rename(avg_CF = "CF", avg_TM = "TM") %>% rename("CellType" = "Cluster") %>% select(-p_val_adj)


expr_merge <- LUMexpr %>% 
  filter(expression ==  "expressed")  %>% 
  select(-expression, -n, -total) %>% 
  unique() %>% 
  pivot_wider(names_from = Type, values_from = perc) %>% 
  rename(perc_CF = "CF", perc_TM = "TM")

all_expression_perc <- expr_merge %>% left_join(avg_merge, by = c("CellType", "TF" = "Gene"))


Motif_Summary_Complete <- Motifs_By_TF %>% inner_join(all_expression_perc) %>% arrange(CellType)





Immune_Combined <- 
Immune_Percent %>% ungroup() %>%
  filter(expression ==  "expressed") %>% 
  select(-expression, -n, -total) %>% 
  unique() %>%
  pivot_wider(names_from = Type, values_from = perc) %>% 
  left_join(Immune_Average_Expression, by = c("TF" = "Gene", "CellType")) %>% 
  rename(perc_CF = "CF", perc_TM = "TM") %>%
  pivot_wider(names_from = Type, values_from = value) %>%
  rename(avg_CF = "CF", avg_TM = "TM")






