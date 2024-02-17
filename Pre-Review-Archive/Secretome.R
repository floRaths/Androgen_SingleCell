library(readxl)
setwd("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/")

options(tibble.print_max = 75, tibble.print_min = 50)

# read in nuclear receptor signaling tables
sig.read <- function(n) {
  Signalingpathways_AR <- read_excel("utilities/Signalingpathways_AR.xlsx")  %>% dplyr::select(Family = "Family", Gene = "Gene", Percentile = "Percentile", Disc.Rate = `Discovery Rate`, GMFC = "GMFC") %>% mutate(Family = "AR") %>% head(n)
  SignalingpathwaysESR <- read_excel("utilities/Signalingpathways_ESR.xlsx") %>% dplyr::select(Family = "Family", Gene = "Gene", Percentile = "Percentile", Disc.Rate = `Discovery Rate`, GMFC = "GMFC") %>% mutate(Family = "ESR") %>% head(n)
  SignalingpathwaysPGR <- read_excel("utilities/Signalingpathways_PGR.xlsx") %>% dplyr::select(Family = "Family", Gene = "Gene", Percentile = "Percentile", Disc.Rate = `Discovery Rate`, GMFC = "GMFC") %>% mutate(Family = "PGR") %>% head(n)
  Sign <- rbind(Signalingpathways_AR, SignalingpathwaysESR, SignalingpathwaysPGR)
  
  print(Sign)
  return(Sign)
  
}
Sign <- sig.read(150)

# read in Human secretome dataset
Human_Secretome <- read_excel("utilities/Human_Secretome_Annotated.xlsx")
Human_Secretome <-  dplyr::select(Human_Secretome, Gene = "gene", Category = "Category") #%>% mutate(Secretion = Category) %>% mutate(Secretion = !is.na(Secretion))

# read in TM-CF response marker list
Marker_Bind <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Output/TM-CF_Response_CellType_MAST_list.rds")
Marker_Bind <- do.call(rbind.data.frame, Marker_Bind) %>% dplyr::select(Gene = "gene", 2:6, CellType = "Cluster")

genes <- 
  Marker_Bind %>% group_by(Gene) %>% summarise(events = n()) %>% arrange(desc(events)) 
Marker_Bind <- inner_join(ungroup(Marker_Bind), genes) %>% arrange(desc(events))


CTPlus_Bind <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/utilities/CellType_Plus_Bind.rds")

genes <- 
  CTPlus_Bind %>% filter(p_val_adj <= 0.05) %>% group_by(Gene) %>% summarise(events = n()) %>% arrange(desc(events)) 
CTPlus_Bind <- inner_join(ungroup(CTPlus_Bind), genes) %>% arrange(desc(events))

rm(genes)



Corr_Bind   <- readRDS("Output/Correlation_Annotated.rds")






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
  as_tibble() #%>% 
  #filter(HR_Rcptr %in% c("AR", "ESR", "PGR"), CellType_DE %in% c("LUM_HR-pos", "Adipocyte")) %>% 
  #select(-CellType_NW) %>%
  #filter(!is.na(Category)) %>%
  #arrange(desc(DE_FC))




secr <- 
left_join(rownames_to_column(AR_Markers, "Gene"), Human_Secretome, by = "Gene") %>% 
  filter(avg_logFC >=0.5, p_val_adj <=0.05) %>% arrange(Category) %>% head(30) %>% pull(Gene)



left_join(rownames_to_column(AR_Markers, "Gene"), Human_Secretome, by = "Gene") %>% 
  filter(avg_logFC >=0.25) %>% arrange(Category) %>% head(30) %>% 
  arrange(-avg_logFC) %>% dplyr::select(-c("p_val", "pct.1", "pct.2"))


display <- 
filter(weighted_networks$lr_sig, weight >= 0.2, from %in% secr) %>% 
  pull(to) %>% unique() %>% intersect(colnames)


receptors_oi <- 
weighted_networks$lr_sig %>% 
  filter(from %in% secr) %>%
  group_by(from) %>% 
  top_n(10, weight) %>% pull(to) %>% unique() %>% as.data.frame() %>% 
  dplyr::select(Receptor = ".") %>% 
  arrange(Receptor) #%>%
  #write_csv("utilities/receptors_oi.csv")
  









left_join(Corr_Bind, Human_Secretome, by = "Gene") %>% 
  filter(CellType == "Fibroblast") %>% 
  arrange(desc(Category)) 

left_join(Marker_Bind, Human_Secretome, by = "Gene") %>% 
  filter(Cluster == "Fibroblast") %>% 
  arrange(desc(Category))













Immu_Bind <- do.call(rbind.data.frame, Marker_list)









