AUC_Summary <- 
LUM@meta.data %>% group_by(CellType, Type) %>% 
  summarise_at(vars(colnames(LUM@meta.data[22:812])), mean) %>% 
  pivot_longer(cols = 3:793) %>% 
  pivot_wider(names_from = Type, values_from = value) %>%
  mutate(change = TM - CF) %>% 
  #filter(CellType == "Basal") %>% 
  separate(name, into = c("Gene", "input"), extra = "merge") %>% 
  arrange(change) #%>% print(n = 500)


AUC_Summary %>% group_by(CellType, Gene) %>% 
  summarise(mean = mean(change)) %>% 
  arrange(-mean)


### weeding put the weird ones tha only show up in the low input runs
Special05 <- 
AUC_Summary %>% ungroup() %>% dplyr::select(Gene, input) %>% distinct() %>%
  mutate(YES = "YES") %>% pivot_wider(names_from = input, values_from = YES) %>% 
  filter(`05` == "YES" & is.na(`1k`) & is.na(`2k`) & is.na(`4k`) & is.na(`8k`) & is.na(`8k_Ank`)) %>% 
  pull(Gene)

Specials <- 
  AUC_Summary %>% 
  filter(Gene %in% Special05) %>% 
  rbind(filter(AUC_Summary, Gene %in% Special1k))

AUC_Summary <- 
  AUC_Summary %>% filter(input %notin% c("05", "1k")) %>% rbind(Specials)






### first pick all positive and negative appearing TF modules. Then pick the max/min values of each
### then remove TFs that appear as postive AND negative... don't trust these

celltype <- "Endothelial"

AUC_Summary %>% filter(CellType == celltype) %>% ggplot(aes(x = reorder(Gene, change), y = change)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 90))


max <-   
AUC_Summary %>% filter(change > 0, CellType == celltype) %>% 
  group_by(Gene) %>% filter(change == max(change)) 

min <-   
  AUC_Summary %>% filter(change < 0, CellType == celltype) %>% 
  group_by(Gene) %>% filter(change == min(change)) 

doubles <- 
rbind(max, min) %>% group_by(Gene) %>% summarise(count = n()) %>% filter(count == 2) %>% pull(Gene)

rbind(max, min) %>% filter(Gene %notin% doubles) %>%
  ggplot(aes(x = reorder(Gene, change), y = change, color = input)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8), 
                       plot.title = element_text(size = 17, face = "bold")) + 
  ggtitle(paste("SCENIC auc_changes...", celltype)) + 
  geom_hline(yintercept = 0, color = "grey")



#### putting all results into one list
MainTFs    <- vector("list", length = length(levels(Sobj2)))

for (i in 1:length(levels(Sobj2))) {
  max <-   
    AUC_Summary %>% filter(change > 0, CellType == levels(Sobj2)[i]) %>% 
    group_by(Gene) %>% filter(change == max(change)) 
  
  min <-   
    AUC_Summary %>% filter(change < 0, CellType == levels(Sobj2)[i]) %>% 
    group_by(Gene) %>% filter(change == min(change)) 
  
  doubles <- 
    rbind(max, min) %>% group_by(Gene) %>% summarise(count = n()) %>% filter(count == 2) %>% pull(Gene)
  
  MainTFs[[i]] <- rbind(max, min) %>% filter(Gene %notin% doubles)
  
}


MainTFs_Bind <- do.call(rbind.data.frame, MainTFs)


MainTFs_Bind %>% group_by(Gene) %>% mutate(count = n()) %>% filter(count == 10, change > 0) %>% arrange(Gene) %>% print(n = 900) 

### putting back input for AUC scoring
MainTFs_Bind %>% unite(input, Gene, input) %>% filter(CellType == "LUM_HR-pos") %>% pull(input) %>% unique()
MainTFs_Bind %>% unite(input, Gene, input) %>% pull(input) %>% unique()

  