# Load required libraries
library(AUCell)      # AUCell is used to analyze gene set activity from single-cell RNA-seq data
library(GSEABase)    # GSEABase is used for handling gene sets (GMT format)
library(data.table)  # data.table is used for fast data manipulation and reading large files

txtFile <- "BrainMet_filtered.tsv"  # Expression matrix file
gmtfile <- "hallmark_msigdb.gmt"    # Gene sets file

# Read the raw counts file using fread for efficient large file handling
geoData <- fread(txtFile, sep="\t")  # Read tab-separated file into a data table

# Remove the first column (Probe_ID) as it's not needed for expression analysis
#geoData$Probe_ID <- NULL

# Extract gene names from the first column
geneNames <- unname(unlist(geoData[,1, with=FALSE]))

# Convert the remaining data (excluding gene names) into a numeric matrix for computation
exprMatrix <- as.matrix(geoData[,-1, with=FALSE])

# Remove geoData from memory to save space
rm(geoData)

# Display the dimensions of the expression matrix (genes x samples)
dim(exprMatrix)

# Assign row names to the expression matrix using extracted gene names
rownames(exprMatrix) <- geneNames

# Load the gene sets from the GMT file
geneSets <- getGmt(gmtfile)  # Read gene sets in GMT format

# Subset gene sets to only include those present in the expression matrix
geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix))

# Check the number of genes in each gene set
cbind(nGenes(geneSets))

# Rank genes for each cell/sample based on their expression levels
# nCores=1 ensures single-threaded processing (can be increased for parallel processing)
# plotStats=TRUE generates ranking statistics plots
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=TRUE)

# Calculate AUCell scores for each cell/sample based on the ranked genes and gene sets
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

# Extract the AUC matrix (activity scores of gene sets across cells)
AUCmatrix <- as.data.frame(cells_AUC@assays@data$AUC)

# Write the AUC scores to a tab-separated text file
write.table(AUCmatrix, "AUCell_scores.txt", sep="\t", row.names=TRUE, quote=FALSE)
