library(nichenetr)
library(tidyverse)


## Source Files
hnscc_expression =      readRDS(url("https://zenodo.org/record/3260758/files/hnscc_expression.rds"))
ligand_target_matrix =  readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
geneset_oi =    readr::read_tsv(url("https://zenodo.org/record/3260758/files/pemt_signature.txt"), col_names = "gene") %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)]
lr_network =            readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

ligand_target_matrix[1:10,1:10]
head(geneset_oi)


### load in ligand target matrix
#ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))



#hnscc_expression = readRDS(url("https://zenodo.org/record/3260758/files/hnscc_expression.rds"))
expression       = hnscc_expression$expression  # nomalized counts (genes as columns)
sample_info      = hnscc_expression$sample_info # contains meta-information about the cells
rm(hnscc_expression)


#### just pulling sample ids of interest
tumors_remove = c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")

CAF_ids       = sample_info %>% 
  filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `non-cancer cell type` == "CAF") %>% 
  pull(cell)

malignant_ids = sample_info %>% 
  filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `classified  as cancer cell` == 1) %>% 
  pull(cell)

rm(tumors_remove)



# var genes determination?!
expressed_genes_sender   = expression[CAF_ids,]       %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
expressed_genes_receiver = expression[malignant_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()

# Check the number of expressed genes: should be a 'reasonable' number of total expressed genes in a cell type, e.g. between 5000-10000 (and not 500 or 20000)
length(expressed_genes_sender)
length(expressed_genes_receiver)

rm(CAF_ids, malignant_ids)



######################################################################
#### Step 2: Define the gene set of interest and a background of genes


# download p-EMT geneset by Puram et al.
#geneset_oi = readr::read_tsv(url("https://zenodo.org/record/3260758/files/pemt_signature.txt"), col_names = "gene") %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)] # only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.
#head(geneset_oi)

# background_expressed genes are intersection of expressed genes and expresion matrix...don't really see why they are different?!
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)




######################################################################
#### Step 3: Define a set of potential ligands

#lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

# If wanted, users can remove ligand-receptor interactions that were predicted based on protein-protein interactions and only keep ligand-receptor interactions that are described in curated databases. To do this: uncomment following line of code:
# lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")

expressed_ligands   = lr_network %>% pull(from) %>% unique() %>% intersect(expressed_genes_sender)
expressed_receptors = lr_network %>% pull(to)   %>% unique() %>% intersect(expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

rm(expressed_ligands, expressed_receptors, expressed_genes_receiver, expressed_genes_sender)

### all potantial ligands from the expressed network
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)





######################################################################
#### Step 4: Perform NicheNet’s ligand activity analysis on the gene set of interest

ligand_activities = predict_ligand_activities(geneset                    = geneset_oi,                 # what we are hoping to predict
                                              background_expressed_genes = background_expressed_genes, # background prediction
                                              ligand_target_matrix       = ligand_target_matrix, 
                                              potential_ligands          = potential_ligands
                                              )

rm(background_expressed_genes, potential_ligands)


### select the best ligands by pearson correlation
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)



# show histogram of ligand activity scores
ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="forestgreen")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()







######################################################################
#### Step 5: Infer target genes of top-ranked ligands and visualize in a heatmap

active_ligand_target_links_df = best_upstream_ligands %>% 
  lapply(get_weighted_ligand_target_links, 
         geneset              = geneset_oi, 
         ligand_target_matrix = ligand_target_matrix, 
         n = 250) %>% 
  bind_rows()

nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)



active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df     = active_ligand_target_links_df, 
                                                                 ligand_target_matrix = ligand_target_matrix, 
                                                                 cutoff               = 0.25
                                                                 )

nrow(active_ligand_target_links)
head(active_ligand_target_links)







order_ligands     = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets     = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()


p_ligand_target_network = vis_ligand_target %>% 
  make_heatmap_ggplot("Prioritized CAF-ligands","p-EMT genes in malignant cells", 
                      color = "purple",
                      legend_position = "top", 
                      x_axis_position = "top",
                      legend_title = "Regulatory potential") + 
  scale_fill_gradient2(low = "whitesmoke",  
                       high = "purple", 
                       breaks = c(0,0.005,0.01)) + 
  theme(axis.text.x = element_text(face = "italic"))

p_ligand_target_network







# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]




vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Prioritized CAF-ligands","Receptors expressed by malignant cells", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network












