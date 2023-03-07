#script for getting spike in information from the klein paper to apply to the qiu paper
erccData <- read.table('ercc-info.txt', stringsAsFactors = F, sep = '\t') #got from https://github.com/jingshuw/DESCEND_manuscript_source_code which is the github page forthe wang et al paper which made the DESCEND method and the alpha-poisson model of technical noise
loadingVolume <- 1e-9 #got the loading volumes and dilutions of ERCC spike-ins for different papers, including klein et al from the supplementary table 1 of svensson et al, titled "Power Analysis of Single Cell RNA-Sequencing Experiments"
dilution <- 5000
avogadrosConstant <- 6.02214076e+23

mix1Concentrations <- erccData$concentration.in.Mix.1..attomoles.ul. * 1e-18 * 1e+6 / dilution #times concentration by 1e+6 cause its in attomoles per microlitre
names(mix1Concentrations) <- erccData$ERCC.ID
mix2Concentrations <- erccData$concentration.in.Mix.2..attomoles.ul. * 1e-18 * 1e+6 / dilution
names(mix2Concentrations) <- erccData$ERCC.ID

mix1attomolesLoaded <- mix1Concentrations * loadingVolume
mix2attomolesLoaded <- mix2Concentrations * loadingVolume

mix1expMols <- mix1attomolesLoaded * avogadrosConstant
mix2expMols <- mix2attomolesLoaded * avogadrosConstant

ERCC_IDs <- erccData[, 1]

kleinExpressionData <- read.csv('GSM1599501_K562_pure_RNA.csv', stringsAsFactors = F) #from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1599501
geneIDs <- kleinExpressionData[, 1]

cells <- ncol(kleinExpressionData) - 1
probIndices <- which(geneIDs %in% ERCC_IDs)
geneIndices <- which(! geneIDs %in% ERCC_IDs)
kleinERCCs <- t(kleinExpressionData[probIndices, 1:cells + 1])

mix1efficiencies <- sapply(1:cells, function(cell, kleinERCCs, mix1expMols) {
  return(sum(kleinERCCs[cell, ]) / sum(mix1expMols))
}, kleinERCCs, mix1expMols)

mix2efficiencies <- sapply(1:cells, function(cell, kleinERCCs, mix2expMols) {
  return(sum(kleinERCCs[cell, ]) / sum(mix2expMols))
}, kleinERCCs, mix2expMols)


alphas <- mix1efficiencies
totalProbeMoleculesPerCell <-  sum(mix1expMols)

cellSpecificUMIcounts <- sapply(1:cells + 1, function(cell, kleinExpressionData, geneIndices) {
  return(sum(kleinExpressionData[geneIndices, cell]))
}, kleinExpressionData, geneIndices)
cellTranscriptContent <- cellSpecificUMIcounts / alphas
meanContentPerCell <- mean(cellTranscriptContent)

kleinData <- list(alphas = alphas, cellContents = cellTranscriptContent)

#then calculate the alphas in the qiu data, following the processing of the qiu datasets

expectedTotalCellularMRNA <- mean(kleinData$cellContents)

qiuData <- readRDS('qiuData.Rdata') 
UMIcounts <- qiuData[["UMIcounts"]]
genes <- names(UMIcounts)
cells <- names(UMIcounts[[genes[1]]])


totalUMIsPerCell <- vector('numeric', length(cells))
names(totalUMIsPerCell) <- cells
counter <- 1
for (gene in genes) {
  totalUMIsPerCell <- totalUMIsPerCell + unlist(UMIcounts[[gene]])
  counter <- counter + 1
}

qiuAlphas <- totalUMIsPerCell / expectedTotalCellularMRNA


saveRDS(qiuAlphas, 'qiuAlphas.Rdata')



