UMIcountData <- new.env()#list()
conversionData <- new.env()#list()
genomicTsData <- new.env()#list()
connection <- file('SRR11683995_data.txt', 'r')
line <- readLines(connection, n = 1)
counter <- 1
while(length(line) > 0) {
  if (counter %% 100000 == 0) {
    print(counter)
  }
  splitLine <- strsplit(line, '\t')[[1]]
  gene <- splitLine[1]
  cell <- splitLine[2]
  if (is.null(UMIcountData[[gene]])) {
    UMIcountData[[gene]] <- list()
    conversionData[[gene]] <- list()
    genomicTsData[[gene]] <- list()
  }
  UMIcount <- as.integer(splitLine[3])
  UMIcountData[[gene]][[cell]] <- UMIcount
  if (UMIcount > 0) {
    conversions <- splitLine[4:length(splitLine)]
    splitConversions <- strsplit(conversions, '/')
    mismatches <- sapply(1:length(splitConversions), function(index, splitConversions) {return(as.integer(splitConversions[[index]][1]))}, splitConversions)
    totalTs <- sapply(1:length(splitConversions), function(index, splitConversions) {return(as.integer(splitConversions[[index]][2]))}, splitConversions)
    conversionData[[gene]][[cell]] <- mismatches
    genomicTsData[[gene]][[cell]] <- totalTs
  }
  counter <- counter + 1
  line <- readLines(connection, n = 1)
}
UMIcountData <- sapply(ls(UMIcountData), function(gene, UMIcountData) {
  return(UMIcountData[[gene]])
}, UMIcountData, simplify = F)
conversionData <- sapply(ls(conversionData), function(gene, conversionData) {
  return(conversionData[[gene]])
}, conversionData, simplify = F)
genomicTsData <- sapply(ls(genomicTsData), function(gene, genomicTsData) {
  return(genomicTsData[[gene]])
}, genomicTsData, simplify = F)

qiuData <- list('UMIcounts' = UMIcountData, 'conversions' = conversionData, 'genomicTs' = genomicTsData)
saveRDS(qiuData, 'qiuData.Rdata')


