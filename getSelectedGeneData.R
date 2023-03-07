qiuData <- readRDS('qiuData.Rdata') #we have 24208 of the total genes in the gtf here simply because any genes with 0 reads across all cells were cut out
UMIcounts <- qiuData[["UMIcounts"]]
genomicTs <- qiuData[["genomicTs"]]
conversions <- qiuData[["conversions"]]
UMIsPerGene <- sapply(UMIcounts, function(geneUMIs) {
  return(sum(unlist(geneUMIs)))
})
qiuGenes <- names(UMIcounts)

#get the T>C rates and CVs in the 4sU dataset
conversionRates <- unlist(sapply(qiuGenes, function(gene, conversions, genomicTs) { #estimated underlying conversion rate
  return(sum(unlist(conversions[[gene]])) / sum(unlist(genomicTs[[gene]])))
}, conversions, genomicTs))
conversionRateCoVs <- sapply(qiuGenes, function(gene, conversions, genomicTs, conversionRates) {
  alpha <- sum(unlist(conversions[[gene]]))
  beta <- sum(unlist(genomicTs[[gene]])) - alpha
  mu <- conversionRates[gene]
  denominator <- (((alpha + beta) ^ 2) * (alpha + beta + 1))
  variance <- (alpha / sqrt(denominator)) * (beta / sqrt(denominator))
  CoV <- sqrt(variance) / mu
  return(CoV)
}, conversions, genomicTs, conversionRates, USE.NAMES = F)
names(conversionRateCoVs) <- qiuGenes

#get the T>C rates and CVs in the control dataset
#we have gene specific backgrounds for 12746 genes simply because any genes with 0 reads across all cells in the control data set were excluded
geneSpecificBackgrounds <- read.table('geneSpecificBackgroundMutationRates.txt', header = F, stringsAsFactors = F)
geneSpecificBackgroundEstimateCOVs <- sapply(1:nrow(geneSpecificBackgrounds), function(index) {
  alpha <- geneSpecificBackgrounds[index, 2]
  beta <- geneSpecificBackgrounds[index, 3] - alpha
  mu <- geneSpecificBackgrounds[index, 4]
  denominator <- (((alpha + beta) ^ 2) * (alpha + beta + 1))
  variance <- (alpha / sqrt(denominator)) * (beta / sqrt(denominator))
  CoV <- sqrt(variance) / mu
  return(CoV)
})
names(geneSpecificBackgroundEstimateCOVs) <- geneSpecificBackgrounds[, 1]

geneSpecificBackgroundRates <- geneSpecificBackgrounds[, 4] #mean value estimates for the background rate of each gene
names(geneSpecificBackgroundRates) <- geneSpecificBackgrounds[, 1]
#look at only genes we have at least 1 read for in both datasets and for which we are confident enough in the background conversion rate
#also cut out genes with 0 observed T>C conversions in the 4sU + chemical conversion dataset
genes <- qiuGenes[which(qiuGenes %in% names(geneSpecificBackgroundRates))]
genes <- genes[which(conversionRates[genes] > 0 & geneSpecificBackgroundRates[genes] > 0)]


selectedUMIsPerGene <- UMIsPerGene[genes]
selectedUMIcountCoVs <- sapply(genes, function(gene) {
  TCs <- unlist(UMIcounts[[gene]])
  return(sd(TCs) / mean(TCs))
})

selectedUMIsPerGenePerCell <- selectedUMIsPerGene / 795


maxCoV <- 10 ^ (-0.5)
genes <- genes[which(geneSpecificBackgroundEstimateCOVs[genes] < maxCoV & conversionRateCoVs[genes] < maxCoV)]
conversionCountsPerGene <- sapply(genes, function(gene, conversions) {
  return(sum(unlist(conversions[[gene]])))
}, conversions)
genomicTCountsPerGene <- sapply(genes, function(gene, genomicTs) {
  return(sum(unlist(genomicTs[[gene]])))
}, genomicTs)
Ps <- seq(0, 0.1, 0.00001)
#can select based on only having <= 10 genes with with probability of 10 / number of genes
improbableGenes <- c()
lowestLikelihood <- c()
for (p in Ps) {
  assumeAllNewmRNAsConversionProbabilities <- pbinom(conversionCountsPerGene - 1, genomicTCountsPerGene, p + geneSpecificBackgroundRates[genes], lower.tail = F, log.p = T)
  improbableGenes <- c(improbableGenes, length(which(assumeAllNewmRNAsConversionProbabilities < log(10 / length(genes)))))
}


p <- Ps[min(which(improbableGenes < 10))] #0.07547


assumeAllNewmRNAsConversionProbabilities <- pbinom(conversionCountsPerGene - 1, genomicTCountsPerGene, p + geneSpecificBackgroundRates[genes], lower.tail = F, log.p = T)


genes <- qiuGenes[which(qiuGenes %in% names(geneSpecificBackgroundRates))]
genes <- genes[which(conversionRates[genes] > 0 & geneSpecificBackgroundRates[genes] > 0)]
maxCoV <- 10 ^ (-0.5)
#comment out line below to look at genes with apparently too low conversions for the full set of genes
genes <- genes[which(geneSpecificBackgroundEstimateCOVs[genes] < maxCoV & conversionRateCoVs[genes] < maxCoV)]
conversionCountsPerGene <- sapply(genes, function(gene, conversions) {
  return(sum(unlist(conversions[[gene]])))
}, conversions)
genomicTCountsPerGene <- sapply(genes, function(gene, genomicTs) {
  return(sum(unlist(genomicTs[[gene]])))
}, genomicTs)


genes <- qiuGenes[which(qiuGenes %in% names(geneSpecificBackgroundRates))]
genes <- genes[which(conversionRates[genes] > 0 & geneSpecificBackgroundRates[genes] > 0)]

UMIcounts <- UMIcounts[genes]
genomicTs <- genomicTs[genes]
conversions <- conversions[genes]
backgrounds <- geneSpecificBackgroundRates[genes]
qiuDataSelectedGenes <- list('UMIcounts' = UMIcounts, 'conversions' = conversions, 'genomicTs' = genomicTs, 'backgrounds' = backgrounds, 'lambdaN' = p)
saveRDS(qiuDataSelectedGenes, 'qiuDataSelectedGenes.Rdata')




