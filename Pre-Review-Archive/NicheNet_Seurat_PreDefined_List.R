library(nichenetr)
library(Seurat)
library(tidyverse)


######################################################################
#### Step 1: Read the expression data of interacting cells

seuratObj = readRDS("Seurat_Objects/Global_DataSets/Nuclei_All_Samples_Batch_Integrated_MetaAdd.rds")
seuratObj[["integrated"]] <- NULL
DefaultAssay(seuratObj) <- "RNA"


Idents(seuratObj) <- "CellType"
levels(seuratObj) <- c("LUM_HR-pos", "LUM_HR-neg", "Basal", "Endothelial", "Fibroblast", "Adipocyte", "Lymphoid", "Myeloid", "Neuronal", "Vasc_Acc")    
seuratObj$CellType <- Idents(seuratObj)

seuratObj@meta.data %>% head()



seuratObj@meta.data$CellType %>% table()
DimPlot(seuratObj, reduction = "umap")

seuratObj@meta.data$Type %>% table()

DimPlot(seuratObj, reduction = "umap", group.by = "Type", pt.size = 1)


######################################################################
#### Step 2: Read in NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks:

ligand_target_matrix = readRDS("NicheNet/utilities/ligand_target_matrix.rds")
ligand_target_matrix[1:5,1:5]

lr_network = readRDS("NicheNet/utilities/lr_network.rds")
head(lr_network)

weighted_networks = readRDS("NicheNet/utilities/weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

head(weighted_networks$lr_sig)
head(weighted_networks$gr)





######################################################################
#### Step 3: Perform the NicheNet analysis

# 1. Define a “sender/niche” cell population and a “receiver/target” cell population present in your expression data and 
# determine which genes are expressed in both populations

## receiver
receiver = c("LUM_HR-neg", "Basal")
receiver = c("LUM_HR-pos")
expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.10)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender

Idents(seuratObj, cells = cells) <- "AR_Sender"
sender_celltypes = c("AR_Sender")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

expressed_genes_sender = 
  filter(rownames_to_column(AR_Markers, "Gene"), p_val_adj <= 0.05, avg_logFC >= 0.5) %>% pull(Gene)


# 2. Define a gene set of interest: these are the genes in the “receiver/target” cell population that are potentially affected 
# by ligands expressed by interacting cells (e.g. genes differentially expressed upon cell-cell interaction)


# DE_geneset_oi -----------------------------------------------------------
seurat_obj_receiver = subset(seuratObj, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["Type"]])

condition_oi        = "TM"
condition_reference = "CF" 

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, 
                                ident.1 = condition_oi, 
                                ident.2 = condition_reference, 
                                min.pct = 0.10, 
                                test.use = "MAST", 
                                only.pos = T) %>% 
  rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_logFC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]



# geneset_oi from other analysis ------------------------------------------

geneset_oi = filter(TopTable, str_detect(Netw.ID, "^ESR") | str_detect(Netw.ID, "^HER4")) %>% pull(Gene) %>% unique()
geneset_oi = filter(Marker_Bind, CellType == "LUM_HR-pos", p_val_adj <= 0.05, avg_logFC >= 0.5 | avg_logFC <= -0.5) %>% pull(Gene) %>% unique()
geneset_oi = filter(TopTable, CellType %in% c("Basal", "LUM_HR-neg")) %>% pull(Gene) %>% unique()

geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]



# 3. Define a set of potential ligands: these are ligands that are expressed by the “sender/niche” cell population and 
# bind a (putative) receptor expressed by the “receiver/target” population

expressed_ligands   = lr_network %>% pull(from) %>% unique() %>% intersect(expressed_genes_sender)
expressed_receptors = lr_network %>% pull(to)   %>% unique() %>% intersect(expressed_genes_receiver)

potential_ligands   = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()







# 4. Perform NicheNet ligand activity analysis: rank the potential ligands based on the presence of their 
# target genes in the gene set of interest (compared to the background set of genes)

ligand_activities = predict_ligand_activities(geneset                    = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix       = ligand_target_matrix, 
                                              potential_ligands          = potential_ligands
                                              )


ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities




best_upstream_ligands = 
  ligand_activities  %>% 
  top_n(25, pearson) %>% 
  arrange(-pearson)  %>% 
  pull(test_ligand)  %>% 
  unique()




DotPlot(seuratObj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()




# 5) Infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis



active_ligand_target_links_df = 
  best_upstream_ligands %>% 
  lapply(get_weighted_ligand_target_links,
         geneset              = geneset_oi, 
         ligand_target_matrix = ligand_target_matrix, 
         n = 250) %>% 
  bind_rows() %>% drop_na()


active_ligand_target_links = 
  prepare_ligand_target_visualization(ligand_target_df     = active_ligand_target_links_df, 
                                      ligand_target_matrix = ligand_target_matrix, 
                                      cutoff = 0.33
  )


order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()




p_ligand_target_network = vis_ligand_target %>% 
  make_heatmap_ggplot("Prioritized ligands","Predicted target genes", 
                      color = "purple",
                      legend_position = "top", 
                      x_axis_position = "top",
                      legend_title = "Regulatory potential")  + 
  theme(axis.text.x = element_text(face = "italic", size = 10)) + 
  scale_fill_gradient2(low = "whitesmoke",  high = "purple", 
                       #breaks = c(0,0.006,0.012)
                       )

p_ligand_target_network










#Receptors of top-ranked ligands

lr_network_top          = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()


p_ligand_receptor_network = 
  vis_ligand_receptor_network %>% 
  t() %>% 
  make_heatmap_ggplot("Ligands","Receptors", 
                      color = "mediumvioletred", 
                      x_axis_position = "top",
                      legend_title = "Prior interaction potential")
p_ligand_receptor_network



















