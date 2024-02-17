library(devtools)
library(velocyto.R)
library(SeuratWrappers)
library(Seurat)

### first create a loom file via velocyto.py

ldat <- ReadVelocity (file = "~/Documents/Velocyto/TM-8249/possorted_genome_bam_BP7AR.loom")
bm   <- as.Seurat      (x = ldat)

filter(rownames_to_column(Sobj_integrated@meta.data, "Cell.ID"), Sample == "TM-8249") # find out cell_id attachment

# rename velocity cell_names to proper number from integrated data
y <- separate(rownames_to_column(bm@meta.data, "Cell.ID"),Cell.ID, into = c("A", "B"), sep = ":")
y$B <- substr(y$B, 1, nchar(y$B)-1)
y <- y %>% mutate(C = "18") %>% unite(B, C, col = "Cell.ID") ### replace number with correct value for sample

# rename all cells
bm <- RenameCells(bm, new.names = y$Cell.ID)

bm <- SCTransform    (bm, assay = "spliced")
bm <- RunPCA         (bm, verbose = FALSE)

# Create Subset from whole data
Subset <- subset(Sobj_integrated, 
                 cells = intersect(rownames(bm@meta.data), 
                                   rownames(Sobj_integrated@meta.data)))

#bm <- subset(bm, cells = intersect(rownames(bm@meta.data), 
#                                   rownames(Sobj_integrated@meta.data)))


Subset <- AddMetaData(Subset, metadata = bm@meta.data) # add velocity metadata

# add velocity assays to integrated dataset
Subset@assays$spliced   <- bm@assays$spliced
Subset@assays$unspliced <- bm@assays$unspliced
Subset@assays$ambiguous <- bm@assays$ambiguous
Subset@reductions$pca   <- bm@reductions$pca



#bm <- SCTransform    (bm, assay = "spliced")
#bm <- RunPCA         (bm, verbose = FALSE)
#bm <- FindNeighbors  (bm, dims = 1:50)
#bm <- FindClusters   (bm, resolution = 0.25)
#bm <- RunUMAP        (bm, dims = 1:50)

Subset <- RunVelocity    (Subset, deltaT = 1, kCells = 25, fit.quantile = 0.02)

Idents(Subset) <- "CellType"

ident.colors <-            (scales::hue_pal())(n = length(x = levels(x = Subset)))
names(x = ident.colors) <- levels(x = Subset)
cell.colors <-             ident.colors[Idents(object = Subset)]
names(x = cell.colors) <-  colnames(x = Subset)

plot <- show.velocity.on.embedding.cor(emb = Embeddings(object = Subset, reduction = "umap"), 
                               vel = Tool(object = Subset, slot = "RunVelocity"), 
                               n = 200,
                               scale = "sqrt",
                               cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8,
                               arrow.scale = 3,
                               show.grid.flow = TRUE,
                               min.grid.cell.mass = 0.5,
                               grid.n = 40,
                               arrow.lwd = 1, 
                               do.par = FALSE,
                               cell.border.alpha = 0.1)













