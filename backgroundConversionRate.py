import csv
import scipy
import numpy
import math
from collections import defaultdict
import matplotlib
import matplotlib.pyplot
import sys
import pickle


gtfData = list(csv.reader(open('gencode.v36.primary_assembly.annotation.gtf', 'r'), delimiter='\t')) #this is 1-indexed
gtfData = gtfData[5:]
chromosomes = list(set(line[0] for line in gtfData))

conversionsByCellByGene = pickle.load(open('SRR11683994_conversionsBySelectedCellByGene.pickle', 'rb'))
genes = list(conversionsByCellByGene.keys())
cells = list(conversionsByCellByGene[genes[0]].keys())
#pool all values together and get total proportion of possible Ts that were converted to Cs in the control sample to find the background mutation/error T>C rate
totalConversions, totalTs = 0, 0
#also get the gene-specific conversion rates
geneSpecificBackgrounds = {gene: [0, 0] for gene in genes}
for gene in genes:
    for cell in cells:
        if len(conversionsByCellByGene[gene][cell]) > 0:
            for read in conversionsByCellByGene[gene][cell]:
                totalConversions += read[0]
                totalTs += read[1]
                geneSpecificBackgrounds[gene][0] += read[0]
                geneSpecificBackgrounds[gene][1] += read[1]


backgroundMutationRate = totalConversions / totalTs
open('backgroundMutationRate.txt', 'w').write(str(backgroundMutationRate))

geneSpecificMutationRate = [[gene, geneSpecificBackgrounds[gene][0], geneSpecificBackgrounds[gene][1], geneSpecificBackgrounds[gene][0] / geneSpecificBackgrounds[gene][1]] for gene in genes if geneSpecificBackgrounds[gene][0] > 0 and geneSpecificBackgrounds[gene][1] > 0]
geneSpecificMutationRateString = '\n'.join([gene[0] + '\t' + str(gene[1]) + '\t' + str(gene[2]) + '\t' + str(gene[3]) for gene in geneSpecificMutationRate])
open('geneSpecificBackgroundMutationRates.txt', 'w').write(geneSpecificMutationRateString)
