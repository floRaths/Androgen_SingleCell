library(tidyverse)
library(readxl)
library(pheatmap)


df <- read_excel("Documents/V1.normalized.log.tpm.bio.genes.T.M.N.xlsx", 
                 sheet = "cut.normalized.log.tpm.bio.gene")


### preparing data to sort
df4 <-   df %>%    gather(2:6166, key = "case", value = "logFC", na.rm = TRUE)  %>% 
                   mutate(case = stringr::str_replace(case, "Bladder", "bladder"),
                          case = stringr::str_replace(case, "lich", "lihc"))    %>% 
                   separate(case, c("Num", "tissue", "sample_type"), sep = "-", remove = F)

# create annotation file
annot <- df4 %>% select(case, tissue, sample_type) %>% distinct() %>% column_to_rownames("case")


### loop to produce separate matrices for different sample categories
Names    <- unique(df4$sample_type)
Types    <- vector("list", length = length(Names))
Order    <- vector("list", length = length(Names))

for (i in 1:length(Names)) { 
  Types[[i]] <- filter(df4, sample_type == Names[i]) %>% arrange(tissue)
  Order[[i]] <- unique(Types[[i]]$case)
  Types[[i]] <- Types[[i]] %>% select(-Num, -tissue, -sample_type) %>% spread(key = case, value = logFC) %>% select(Gene = Row.names, Order[[i]])
  }
  

### loop to produce separate matrices for different sample categories
Names    <- unique(df4$tissue)
Tissu    <- vector("list", length = length(Names))
Order    <- vector("list", length = length(Names))

for (i in 1:length(Names)) { 
  Tissu[[i]] <- filter(df4, tissue == Names[i]) %>% arrange(sample_type)
  Order[[i]] <- unique(Tissu[[i]]$case)
  Tissu[[i]] <- Tissu[[i]] %>% select(-Num, -tissue, -sample_type) %>% spread(key = case, value = logFC) %>% select(Gene = Row.names, Order[[i]])
}



### you could also bind the sorted dataframes into one by doing this:
Sort <- left_join(Types[[1]], Types[[2]], by =  "Gene") #%>% left_join(Norm, by =  "Gene")

Sort2 <- sample_n(Tumr, 250)




mat <- as.matrix(column_to_rownames(Tissu[[10]], "Gene"))

pheatmap(mat,
         scale = "row", 
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         color = viridis(100),
         cluster_cols = T, 
         annotation = annot, 
         show_rownames = F, 
         show_colnames = F)


B <- pheatmap(mat,
         scale = "row", 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = viridis(100),
         cluster_cols = F, 
         annotation = annot, 
         show_rownames = F, 
         show_colnames = F)

cowplot::plot_grid(A, B)


ggsave(plot     = A,
       path     = file.path("./Documents/"),
       filename = "Test.png",
       width    = 16, 
       height   = 9,
       scale    = 1.5,
       dpi      = "screen")


install.packages("grid")
install.packages("ggplotify")
library("grid")
library("ggplotify")


A <- as.ggplot(A)
B <- as.ggplot(B)
