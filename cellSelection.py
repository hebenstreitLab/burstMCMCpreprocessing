import csv
import scipy
import numpy
import math
from _collections import defaultdict
import matplotlib
import matplotlib.pyplot
import sys
import pickle

#get total reads in each cell and use this to do the final cell selection

fileName = sys.argv[1] #'SRR11683995'
targetValues = {'SRR11683994': 400, 'SRR11683995': 795}
maxValues = {'SRR11683994': 500, 'SRR11683995': 1000}
targetValue = targetValues[fileName]
maxValue = maxValues[fileName]

conversionsByCellByGene = pickle.load(open(fileName + '_conversionsByCellByGene.pickle', 'rb'))
UMIcountByCellByGene = pickle.load(open(fileName + '_UMIcountByCellByGene.pickle', 'rb'))
genes = list(conversionsByCellByGene.keys())
cells = list(conversionsByCellByGene[genes[0]].keys())
totalReadsPerCell = {cell: 0 for cell in cells}
UMIcountByCell = {cell: 0 for cell in cells}
for gene in genes:
    for cell in cells:
        totalReadsPerCell[cell] += len(conversionsByCellByGene[gene][cell])
        UMIcountByCell[cell] += UMIcountByCellByGene[gene][cell]

sortedReadsPerCell = list(reversed(sorted(zip([totalReadsPerCell[cell] for cell in cells], cells), key=lambda pair: pair[0])))
sortedUMIsPerCell = list(reversed(sorted(zip([UMIcountByCell[cell] for cell in cells], cells), key=lambda pair: pair[0])))


sortedCellsByRead = [sortedReadsPerCell[index][1] for index in range(len(sortedReadsPerCell))]
sortedCellsByUMI = [sortedUMIsPerCell[index][1] for index in range(len(sortedUMIsPerCell))]

sortedReads = [sortedReadsPerCell[index][0] for index in range(len(sortedReadsPerCell))]
sortedUMIs = [sortedUMIsPerCell[index][0] for index in range(len(sortedUMIsPerCell))]

cumsumUMIs = numpy.cumsum(sortedReads) / numpy.sum(sortedReads)
cumsumReads = numpy.cumsum(sortedUMIs) / numpy.sum(sortedUMIs)

matplotlib.pyplot.figure()
matplotlib.pyplot.plot(range(maxValue), 100 * cumsumUMIs)
matplotlib.pyplot.vlines(targetValues[fileName], 0, 100 * cumsumUMIs[-1]) #795 is no. cells specified in the dataset for the treatment (4sU) dataset, which looks to line up nicely with the point of inflection
matplotlib.pyplot.ylabel('Cumsum of transcripts (%)')
matplotlib.pyplot.xlabel('Cells')
matplotlib.pyplot.savefig('final_' + fileName + '_UMIs.pdf')
matplotlib.pyplot.figure()
matplotlib.pyplot.plot(range(maxValue), 100 * cumsumReads)
matplotlib.pyplot.vlines(targetValues[fileName], 0, 100 * cumsumReads[-1])
matplotlib.pyplot.ylabel('Cumsum of reads (%)')
matplotlib.pyplot.xlabel('Cells')
matplotlib.pyplot.savefig('final_' + fileName + '_reads.pdf')

selectedCells = sortedCellsByRead[:targetValue]
conversionsBySelectedCellByGene = {gene: {cell: conversionsByCellByGene[gene][cell] for cell in selectedCells} for gene in genes}
UMIcountBySelectedCellByGene = {gene: {cell: UMIcountByCellByGene[gene][cell] for cell in selectedCells} for gene in genes}
pickle.dump(conversionsBySelectedCellByGene, open(fileName + '_conversionsBySelectedCellByGene.pickle', 'wb'))
pickle.dump(UMIcountBySelectedCellByGene, open(fileName + '_UMIcountBySelectedCellByGene.pickle', 'wb'))
