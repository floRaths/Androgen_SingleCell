library(tidyverse)

geneFiltering2 <- function (exprMat, scenicOptions, minCountsPerGene = 3 * 0.01 * 
                              ncol(exprMat), minSamples = ncol(exprMat) * 0.01) 
{
  outFile_genesKept <- NULL
  dbFilePath <- NULL
  if (class(scenicOptions) == "ScenicOptions") {
    dbFilePath <- getDatabases(scenicOptions)[[1]]
    outFile_genesKept <- getIntName(scenicOptions, "genesKept")
  }
  else {
    dbFilePath <- scenicOptions[["dbFilePath"]]
    outFile_genesKept <- scenicOptions[["outFile_genesKept"]]
  }
  if (is.null(dbFilePath)) 
    stop("dbFilePath")
  if (is.data.frame(exprMat)) {
    supportedClasses <- paste(gsub("AUCell_buildRankings,", 
                                   "", methods("AUCell_buildRankings")), collapse = ", ")
    supportedClasses <- gsub("-method", "", supportedClasses)
    stop("'exprMat' should be one of the following classes: ", 
         supportedClasses, "(data.frames are not supported. Please, convert the expression matrix to one of these classes.)")
  }
  if (any(table(rownames(exprMat)) > 1)) 
    stop("The rownames (gene id/name) in the expression matrix should be unique.")
  nCountsPerGene <- Matrix::rowSums(exprMat, na.rm = T)
  nCellsPerGene <- Matrix::rowSums(exprMat > 0, na.rm = T)
  message("Maximum value in the expression matrix: ", max(exprMat, 
                                                          na.rm = T))
  message("Ratio of detected vs non-detected: ", signif(sum(exprMat > 
                                                              0, na.rm = T)/sum(exprMat == 0, na.rm = T), 2))
  message("Number of counts (in the dataset units) per gene:")
  print(summary(nCountsPerGene))
  message("Number of cells in which each gene is detected:")
  print(summary(nCellsPerGene))
  message("\nNumber of genes left after applying the following filters (sequential):")
  genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > 
                                                      minCountsPerGene)]
  message("\t", length(genesLeft_minReads), "\tgenes with counts per gene > ", 
          minCountsPerGene)
  nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
  genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > 
                                                      minSamples)]
  message("\t", length(genesLeft_minCells), "\tgenes detected in more than ", 
          minSamples, " cells")
  library(RcisTarget)
  motifRankings <- importRankings(dbFilePath)
  genesInDatabase <- colnames(getRanking(motifRankings))
  genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% 
                                                               genesInDatabase)]
  message("\t", length(genesLeft_minCells_inDatabases), "\tgenes available in RcisTarget database")
  genesKept <- genesLeft_minCells_inDatabases
  if (!is.null(outFile_genesKept)) {
    saveRDS(genesKept, file = outFile_genesKept)
    if (getSettings(scenicOptions, "verbose")) 
      message("Gene list saved in ", outFile_genesKept)
  }
  return(genesKept)
}


dir.create("~/Documents/SCENIC/SCENIC_TMCF_down", recursive = T)
     setwd("~/Documents/SCENIC/") # Or `knitr::opts_knit$set(root.dir = 'example_results/SCENIC_MouseBrain')` in the first chunk if running a notebook
     dir.create("int")
     

Idents(Sobj_integrated) <- "Type"
Subset   <- subset(Sobj_integrated, downsample = 500)
exprMat_1k  <- Subset@assays$RNA@counts
cellInfo <- select(Subset@meta.data, CellType = "CellType", nGene = "nFeature_RNA", nUMI = "nCount_RNA")
levels(cellInfo$CellType) <- droplevels(cellInfo$Class)

dim(exprMat)
head(cellInfo)

#cellInfo <- data.frame(cellInfo)

#cellTypeColumn <- "Class"
#colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "CellType"
cbind(table(cellInfo$CellType))

saveRDS(cellInfo, file="int/cellInfo.Rds")


library(ggsci)
palette        <- (pal_d3("category20")(length(levels(cellInfo$CellType))))
names(palette) <- levels(cellInfo$CellType)
colVars        <- list(CellType = palette)


# Color to assign to the variables (same format as for NMF::aheatmap)
#colVars <- list(CellType=c("LUM_HR-neg"     = "forestgreen"
                           #"Endo_2"     = "darkorange", 
                           #"Endo_3"     = "magenta4", 
                           #"Endo_4"     = "hotpink" 
 #                          ))

#colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")


plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))


library(SCENIC)
org             = "hgnc" # or hgnc, or dmel
dbDir           = "~/Documents/SCENIC/cisTarget_databases/" # RcisTarget databases location
myDatasetTitle  = "SCENIC example on TMNCF_small" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic( org = org, dbDir = dbDir, dbs = dbs, datasetTitle = myDatasetTitle, nCores = 6) 



# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars  <- "int/colVars.Rds"
# Databases:
# scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr"="mm9-tss-centered-5kb-10species.mc8nr.feather")
# scenicOptions@settings$db_mcVersion <- "v8"

# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 




# (Adjust minimum values according to your dataset)
genesKept <- geneFiltering2(exprMat, scenicOptions = scenicOptions,
                           minCountsPerGene       = 3*.01*ncol(exprMat),
                           minSamples             = ncol(exprMat)*.01)



exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)

rm(exprMat)


runCorrelation(as.matrix(exprMat_filtered), scenicOptions)


exprMat_filtered <- log2(as.matrix(exprMat_filtered)+1)
saveRDS(exprMat_filteredlog2, "int/exprMat_filteredlog2.rds")

write.csv(Matrix::as.matrix(exprMat), "~/Documents/SCENIC/example/expr_mat_tiny.csv")
saveRDS(exprMat, "int/exprMat.rds")
exprMat_filtered <- readRDS("int/exprMat_filtered.rds")


library(SCENIC)
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

# For a very quick run: 
# coexMethod=c("top5perTarget")
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
# save...

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) #** Only for toy run!!
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)






write_csv(exmat, "~/mnt/SCENIC_TMCF_local/int/expressionmatrix.csv")


build_loom("~/mnt/exmat.loom",
  
           dgem  = as.matrix(x = exprMat),
           title = "Full_TM-CF_Data",
           genome= "Human", # Just for user information, not used internally
           )



library(SCopeLoomR)


SCopeLoomR::build_loom(dgem = Matrix::as.matrix(exprmat),
                       file.name = "int/exprMat_5k.loom", 
                       title = "TMCF_Data_5k_Downsample", 
                       genome = "Human")



















