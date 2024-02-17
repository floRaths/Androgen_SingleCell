library(monocle3)


# Load the data
expression_matrix  <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_expression.rds"))
cell_metadata      <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_colData.rds"))
gene_annotation    <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_rowData.rds"))

# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


#Pre-Processing Data. When using PCA, you should specify the number of principal components you want Monocle to compute.
cds <- preprocess_cds(cds, method = "PCA", num_dim = 100, verbose = T)

plot_pc_variance_explained(cds)

#Reduce dimensionality and visualize the cells
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
cds <- reduce_dimension(cds, reduction_method = "tSNE", preprocess_method = "PCA")

#Plotting cells. metadata found in @colData
plot_cells(cds, cell_size = 1.5,
                group_label_size = 4,
                color_cells_by = "cao_cell_type",
                alpha = 0.5)

# You can also color your cells according to how much of a gene or set of genes they express:
plot_cells(cds, genes = c("cpna-2", "egl-21", "ram-2", "inos-1"),
                cell_size = 1.5,
                alpha = 1)

# or from which batch they were processed to spot batch effects
plot_cells(cds, color_cells_by = "plate", label_cell_groups=F, cell_size = 1.5)


# No dramatic batch effects are seen, nevertheless, we can try and remove that batch effect is by running the align_cds() function:
cds = align_cds(cds, num_dim = 100, alignment_group = "plate")

cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")

plot_cells(cds, color_cells_by = "plate", label_cell_groups=F, cell_size = 1.5)



# Monocle uses a technique called community detection to group cells. 
cds = cluster_cells(cds, resolution = 1e-5, reduction_method = "UMAP")
plot_cells(cds, cell_size = 1, group_label_size = 5)

# The cluster_cells() also divides the cells into larger, more well separated groups called partitions 
# You can visualize these partitions like this:
plot_cells(cds, color_cells_by = "partition", group_cells_by = "partition", cell_size = 1, group_label_size = 5)


# after running clutser_cells, the color_cells argument in plot_cells determines the labels automatically
plot_cells(cds, color_cells_by = "cao_tissue",    cell_size = 1, group_label_size = 5)
plot_cells(cds, color_cells_by = "cao_cell_type", cell_size = 1, group_label_size = 5)

# You can choose to label whole partitions instead of clusters by passing group_cells_by="partition". 
# You can also plot the top 2 labels per cluster by passing labels_per_group=2 to plot_cells(). 
plot_cells(cds, color_cells_by = "cao_tissue", group_cells_by="partition", cell_size = 1, group_label_size = 5)


# Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells(), 
# like this:
  
plot_cells(cds, color_cells_by="cao_cell_type", label_groups_by_cluster=FALSE, cell_size = 1, group_label_size = 5)



# Find marker genes expressed by each cluster
marker_test_res <- top_markers(cds, group_cells_by="partition", 
                               reference_cells=1000, cores=8)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))


plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=3)


top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="cluster_row_col",
                    max.size=3)



# Annotate your cells according to type
# To assign cell types based on clustering, we begin by creating a new column in colData(cds) 
# and initialize it with the values of clusters(cds):

colData(cds)$assigned_cell_type <- as.character(partitions(cds))

# Now, we can use the dplyrpackage's recode() function to remap each cluster to a different cell type:
colData(cds)$assigned_cell_type = dplyr::recode(colData(cds)$assigned_cell_type,
                                                "1"="Germline",
                                                "2"="Body wall muscle",
                                                "3"="Unclassified neurons",
                                                "4"="Vulval precursors",
                                                "5"="Failed QC",
                                                "6"="Seam cells",
                                                "7"="Pharyngeal epithelia",
                                                "8"="Coelomocytes",
                                                "9"="Am/PH sheath cells",
                                                "10"="Failed QC",
                                                "11"="Touch receptor neurons",
                                                "12"="Intestinal/rectal muscle",
                                                "13"="Pharyngeal neurons",
                                                "14"="Unknown2",
                                                "15"="flp-1(+) interneurons",
                                                "16"="Canal associated neurons",
                                                "17"="Ciliated sensory neurons",
                                                "18"="Other interneurons",
                                                "19"="Pharyngeal gland",
                                                "20"="Failed QC",
                                                "21"="Ciliated sensory neurons",
                                                "22"="Oxygen sensory neurons",
                                                "23"="Ciliated sensory neurons",
                                                "24"="Ciliated sensory neurons",
                                                "25"="Ciliated sensory neurons",
                                                "26"="Ciliated sensory neurons",
                                                "27"="Oxygen sensory neurons",
                                                "28"="Ciliated sensory neurons",
                                                "29"="Unclassified neurons",
                                                "30"="Socket cells",
                                                "31"="Failed QC",
                                                "32"="Pharyngeal gland",
                                                "33"="Ciliated sensory neurons",
                                                "34"="Ciliated sensory neurons",
                                                "35"="Ciliated sensory neurons",
                                                "36"="Failed QC",
                                                "37"="Ciliated sensory neurons",
                                                "38"="Pharyngeal muscle",
                                                "39"="Unknown")


plot_cells(cds, group_cells_by="partition", color_cells_by="assigned_cell_type", group_label_size = 5)





# Partition 7 has some substructure, and it's not obvious just from looking at the output of top_markers() 
# what cell type or types it corresponds to. So we can isolate it with the choose_cells() function for further analysis:
cds_subset <- choose_cells(cds)


# Now we have a smaller cell_data_set object that contains just the cells from the partition we'd like to drill into. 
# We can use graph_test() to identify genes that are differentially expressed in different subsets of cells from this partition:

pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))

# We can take all the genes that vary across this set of cells and group those that have similar patterns of expression into modules:
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=1e-3)

plot_cells(cds_subset, genes=gene_module_df, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE,
           cell_size = 1)




cds_subset = cluster_cells(cds_subset, resolution=1e-2)
plot_cells(cds_subset, color_cells_by="cluster", cell_size = 1)


# Based on how the patterns line up, we'll make the following assignments:

colData(cds_subset)$assigned_cell_type <- as.character(clusters(cds_subset)[colnames(cds_subset)])

colData(cds_subset)$assigned_cell_type <- dplyr::recode(colData(cds_subset)$assigned_cell_type,
                                                        "1"="Somatic gonad precursors",
                                                        "2"="Somatic gonad precursors",
                                                        "3"="Vulval precursors",
                                                        "4"="Sex myoblasts",
                                                        "5"="Sex myoblasts",
                                                        "6"="Vulval precursors",
                                                        "7"="Failed QC",
                                                        "8"="Vulval precursors",
                                                        "10"="Unclassified neurons",
                                                        "11"="Distal tip cells")

plot_cells(cds_subset, group_cells_by="cluster", color_cells_by="assigned_cell_type")



# Now we can transfer the annotations from the cds_subset object back to the full dataset. 
# We'll also filter out low-quality cells at this stage

colData(cds)[colnames(cds_subset),]$assigned_cell_type <- colData(cds_subset)$assigned_cell_type

cds <- cds[,colData(cds)$assigned_cell_type != "Failed QC" | is.na(colData(cds)$assigned_cell_type )]

plot_cells(cds, group_cells_by="partition", 
           color_cells_by="assigned_cell_type", 
           labels_per_group=5)



# To generate a Garnett file, first find the top markers that each annotated cell type expresses:

assigned_type_marker_test_res <- top_markers(cds,
                                             group_cells_by="assigned_cell_type",
                                             reference_cells=1000,
                                             cores=8)


# Next, filter these markers according to how stringent you want to be:
  
  # Require that markers have at least JS specificty score > 0.5 and
  # be significant in the logistic test for identifying their cell type:

garnett_markers <- assigned_type_marker_test_res %>%
  filter(marker_test_q_value < 0.01 & specificity >= 0.5) %>%
  group_by(cell_group) %>%
  top_n(5, marker_score)

# Exclude genes that are good markers for more than one cell type:
garnett_markers <- garnett_markers %>% 
  group_by(gene_short_name) %>%
  filter(n() == 1)



# Then call generate_garnett_marker_file:
  
generate_garnett_marker_file(garnett_markers, file = "./marker_file.txt")























