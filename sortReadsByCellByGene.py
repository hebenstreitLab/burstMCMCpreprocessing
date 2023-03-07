import csv
import scipy
import numpy
import math
from _collections import defaultdict
import matplotlib
import matplotlib.pyplot
import sys
import pickle
import Bio
from Bio import SeqIO

#here we collapse UMIs belonging to the same gene and cell that have a hamming distance of 1 or less and record number of UMIs per cell per gene,
#and then convert read information into the info required by the algorithm, which is to record the number of T>C conversions in the read and the number of Ts in the genomic sequence

fileName = sys.argv[1] 
readInfoByCell = pickle.load(open(fileName + '_readInfoByCell.pickle', 'rb'))
readInfoByGene = pickle.load(open(fileName + '_readInfoByGene.pickle', 'rb'))
cells = list(readInfoByCell.keys())
genes = list(readInfoByGene.keys()) #try cutting out all genes that have no reads across any cells to save on RAM
readCountsByGene = [len(readInfoByGene[gene]) for gene in genes]
genes = [genes[index] for index in range(len(genes)) if readCountsByGene[index] > 0]

#sort readIDs and UMIs by cell by gene to make quicker access for later
readIDsByCellByGene = {gene: {cell: [] for cell in cells} for gene in genes}
UMIsByCellByGene = {gene: {cell: [] for cell in cells} for gene in genes}
counter = 0
for cell in cells:
    cellReadIDs = set(readInfoByCell[cell].keys())
    counter += 1
    print(counter)
    for gene in genes:
        geneReadIDs = set(readInfoByGene[gene].keys())
        # get overlap of reads found in both cell and gene
        readIDs = list(geneReadIDs.intersection(cellReadIDs)) 
        if len(readIDs) > 0:
            readIDsByCellByGene[gene][cell] = readIDs
            UMIs = list(set([readInfoByCell[cell][readID][2] for readID in readIDs]))
            UMIsByCellByGene[gene][cell] = UMIs


def hammingDistance(A, B): #A and B must be strings of the same length
    return sum(1 for a, b in zip(A, B) if a != b)


def collapseUMIs(UMIs, hammingDistance): #algorithm to optimally collapse UMIs, is to first find the row with most 1s in hamming matrix, then all UMIs 1 away from that are collapsed to it, then reform the hamming matrix and repeat until no 1s present in matrix
    hammingMatrix = numpy.array([[hammingDistance(UMI_A, UMI_B) for UMI_B in UMIs] for UMI_A in UMIs])
    rowIndexOneCounts = numpy.count_nonzero(hammingMatrix == 1, axis = 0)
    while any(rowIndexOneCounts > 0):
        bestIndex = numpy.argmax(rowIndexOneCounts)
        collapsingIndices = [index for index in range(len(hammingMatrix[bestIndex])) if hammingMatrix[bestIndex][index] == 1]
        #instead of recalculating all remaining distances and reforming matrix, just delete rows and columns as necessary
        deletionCounter = 0
        for index in collapsingIndices:
            hammingMatrix = numpy.delete(hammingMatrix, index - deletionCounter, 0)
            hammingMatrix = numpy.delete(hammingMatrix, index - deletionCounter, 1)
            del(UMIs[index - deletionCounter])
            deletionCounter += 1
        rowIndexOneCounts = numpy.count_nonzero(hammingMatrix == 1, axis=0)
    return UMIs

def getConversions(readSequence, genomicSequence): #get the number of T>C conversions in the read and the total number of Ts in the genomic sequence
    conversions = len([index for index in range(len(readSequence)) if genomicSequence[index] == 'T' and readSequence[index] == 'C'])
    genomicTs = len([index for index in range(len(genomicSequence)) if genomicSequence[index] == 'T'])
    return conversions, genomicTs


UMIcountByCellByGene = {gene: {cell: 0 for cell in cells} for gene in genes} #count of total UMIs observed in particular cell for a particular gene after collapsing those with hamming distance of 1
conversionsByCellByGene = {gene: {cell: [] for cell in cells} for gene in genes} #record the observed T>C conversions in the read and total number of Ts in the genomic sequence


counter = 0
for gene in genes:
    if counter % 1000 == 0:
        print(counter)
    counter += 1
    for cell in cells:
        # get overlap of reads found in both cell and gene
        readIDs = readIDsByCellByGene[gene][cell]
        if len(readIDs) > 0:
            UMIs = UMIsByCellByGene[gene][cell][:] # the [:] means that changes to the UMIs, when i delete some to merge based on hamming distance, will not affect the dictionary UMIsByCellByGene, because it causes UMIs to be a copy of the content rather than simply a reference to the dictionary
            collapsedUMIs = collapseUMIs(UMIs, hammingDistance)
            UMIcountByCellByGene[gene][cell] = len(collapsedUMIs)
            for readID in readIDs:
                readInfo = readInfoByCell[cell][readID]  # same as readInfoByGene[gene][readID]
                readSequence, genomicSequence = readInfo[3:]
                conversions, genomicTs = getConversions(readSequence, genomicSequence)
                conversionsByCellByGene[gene][cell] = conversionsByCellByGene[gene][cell] + [(conversions, genomicTs)]

del(readInfoByCell)
del(readInfoByGene)
del(readIDsByCellByGene)
del(UMIsByCellByGene)

pickle.dump(UMIcountByCellByGene, open(fileName + '_UMIcountByCellByGene.pickle', 'wb'))
pickle.dump(conversionsByCellByGene, open(fileName + '_conversionsByCellByGene.pickle', 'wb'))
