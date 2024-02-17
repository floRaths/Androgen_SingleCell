library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
set.seed(1234)


setwd("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/CellRanger_Out/ATAC/")

counts <- Read10X_h5(filename = "TM-3937_ATAC/filtered_peak_bc_matrix.h5")

### explore the region size distribution
rownames(counts) %>% 
  as_tibble() %>% 
  separate("value", into = c("A", "B", "C"), convert = T) %>% 
  mutate(diff = C - B) %>% dplyr::filter(diff < 5000) %>% ### few very large fragments are filtered out
  ggplot(aes(x = diff)) + geom_histogram(bins = 250)


metadata <- read.csv(
  file = "TM-3937_ATAC/singlecell.csv",
  header = TRUE,
  row.names = 1
)

pbmc <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata
)


fragment.path <- "TM-3937_ATAC/fragments.tsv.gz"

pbmc <- SetFragments(
  object = pbmc,
  file = fragment.path
)



pbmc <- NucleosomeSignal(object = pbmc)


pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio    <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

p1 <- VlnPlot(pbmc, c('pct_reads_in_peaks', 'peak_region_fragments'), pt.size = 0.1)
p2 <- VlnPlot(pbmc, c('blacklist_ratio', 'nucleosome_signal'), pt.size = 0.1) & scale_y_log10()

p1 | p2


pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')


# create granges object with TSS positions
gene.ranges <- genes(EnsDb.Hsapiens.v75)
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

tss.ranges <- GRanges(
  seqnames = seqnames(gene.ranges),
  ranges = IRanges(start = start(gene.ranges), width = 2),
  strand = strand(gene.ranges)
)

seqlevelsStyle(tss.ranges) <- 'UCSC'
tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')

# to save time use the first 2000 TSSs
pbmc <- TSSEnrichment(object = pbmc, tss.positions = tss.ranges[1:2000])
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(pbmc, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()


pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 10 &
    TSS.enrichment > 2
)

pbmc




pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(
  object = pbmc,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)


DepthCor(pbmc)



pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3, resolution = 0.3)

DimPlot(object = pbmc, label = TRUE) + NoLegend()





# Extend coordinates upstream to include the promoter
genebodyandpromoter.coords <- Extend(x = gene.ranges, upstream = 2000, downstream = 0)

# create a gene by cell matrix
gene.activities <- FeatureMatrix(
  fragments = fragment.path,
  features = genebodyandpromoter.coords,
  cells = colnames(pbmc),
  chunk = 20
)

# convert rownames from chromsomal coordinates into gene names
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- gene.key[rownames(gene.activities)]
# add the gene activity matrix to the Seurat object as a new assay, and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)



DefaultAssay(pbmc) <- 'RNA'

FeaturePlot(object = pbmc, features = c('ESR1', 'PGR', 'CUX2', 'RYR2', 'AR', 'PTPRC'),
  pt.size = 0.8,
  max.cutoff = 'q80', order = F,
  ncol = 3,
  cols = viridis::magma(n = 100))






transfer.anchors <- FindTransferAnchors(
  reference = Sobj_TM,
  query = pbmc,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)




