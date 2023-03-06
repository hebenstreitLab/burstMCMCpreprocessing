import csv
import scipy
import numpy
import math
from _collections import defaultdict
import matplotlib
import matplotlib.pyplot
import sys
import pickle

#for detecting and repairing synthesis errors in the trimmed barcodes, then selecting cell barcodes as "real cells"
#run from command line with python3 processTrimmedBarcodes.py SRR11683994 & and python3 processTrimmedBarcodes.py SRR11683995 &

fileName = sys.argv[1] #'SRR11683994'
print(fileName)
filePath = fileName + '_1.trimmed.fastq'
sequences = []
readIDs = []
i = 0
with open(filePath) as fastqFile:
    for line in fastqFile:
        if ((i + 4) / 4).is_integer():
            # start new variable for new read and add first bit of content
            readContents = [line]	
        elif ((i + 1) / 4).is_integer():
            # add last bit of read content
            readContents.append(line)
            sequences.append(readContents[1][:-1])
            readIDs.append(readContents[0].split(' ')[1].split('/')[0])
        else:
            readContents.append(line)
        if i % 1000000 == 0:
            print(i)
        i += 1
fastqFile.close()

cellIDs = []
UMIs = []
for barcode in sequences:
    cellID = barcode[:12]
    UMI = barcode[12:]
    cellIDs.append(cellID)
    UMIs.append(UMI)
del(sequences)#save on RAM

UMIsByCellID = {cellID: [] for cellID in set(cellIDs)}
UMIsByReadByByCell = {cellID: defaultdict(str) for cellID in set(cellIDs)}
for index in range(len(cellIDs)):
    cellID = cellIDs[index]
    UMI = UMIs[index]
    UMIsByCellID[cellID].append(UMI)
    readID = readIDs[index]
    UMIsByReadByByCell[cellID][readID] = UMI
del(cellIDs)
del(UMIs)
del(readIDs)

#only keep cellIDs with >=25 different UMIs
UMIsByReadByByCell = {cellID: UMIsByReadByByCell[cellID] for cellID in UMIsByCellID.keys() if len(set(UMIsByCellID[cellID])) >= 25}
UMIsByCellID = {cellID: UMIsByCellID[cellID] for cellID in UMIsByCellID.keys() if len(set(UMIsByCellID[cellID])) >= 25}

#fix the case where a cell barcode has a missing base, which is detecting by fixing of Ts at the last base of the UMI
#in these cases, add the last base of the cell ID to the start of the UMIs, and replace it in the cell barcode with an N while clipping the last base, which will be a T, from each UMI
repairedCellIDMapper = defaultdict(list) #keep track of which repaired cellIDs correspond to which original ones
counter = 0
cellIDs = list(UMIsByCellID.keys())
for cellID in cellIDs:
    UMIs = UMIsByCellID[cellID]
    if 0.8 <= numpy.sum([1 for UMI in UMIs if UMI[-1] == 'T']) / len(UMIs): #check if cell barcode has missing base
        repairedCellIDMapper[cellID] = cellID[:-1] + 'N'
        UMIs = [cellID[-1] + UMI[:-1] for UMI in UMIs]
        del(UMIsByCellID[cellID])
        if cellID[:-1] + 'N' in UMIsByCellID.keys(): #here we naturally merge cell barcodes with the same first 11bp if they had one base missing
            for UMI in UMIs:
                UMIsByCellID[cellID[:-1] + 'N'].append(UMI)
            for read in UMIsByReadByByCell[cellID].keys():
                UMIsByReadByByCell[cellID[:-1] + 'N'][read] = cellID[-1] + UMIsByReadByByCell[cellID][read][:-1]
        else:
            UMIsByCellID[cellID[:-1] + 'N'] = UMIs
            UMIsByReadByByCell[cellID[:-1] + 'N'] = {read: cellID[-1] + UMIsByReadByByCell[cellID][read][:-1] for read in UMIsByReadByByCell[cellID].keys()}
        del(UMIsByReadByByCell[cellID])
        counter += 1

#next thing is to get rid of any cell barcodes where there is a base present at a position in >=80% of the UMIs corresponding to that cell barcode
counter = 0
cellIDs = list(UMIsByCellID.keys())
nucleotides = ['A', 'T', 'C', 'G']
for cellID in cellIDs:
    UMIs = set(UMIsByCellID[cellID])
    delete = False
    for position in range(8):
        bases = [UMI[position] for UMI in UMIs]
        counts = numpy.array([len([base for base in bases if base == nucleotide]) for nucleotide in nucleotides])
        proportions = counts / len(UMIs)
        if len(numpy.argwhere(proportions >= 0.8)) == 1:
            delete = True
    if delete:
        del(UMIsByCellID[cellID])
        del(UMIsByReadByByCell[cellID])
        counter += 1
        if cellID in repairedCellIDMapper.keys():
            del(repairedCellIDMapper[cellID])


#now the required barcode processing is completed, the next thing to do is merging of cell barcodes based on either number of genes, transcripts (UMIs) or reads per cell
#where you get the cumulative sum of the reads/transcripts/genes in the cells and get all cells left of the point of inflection since those are supposed to be 'STAMPS' instead of 'empties'
cellIDs = list(UMIsByCellID.keys())
UMIsPerCell = numpy.array([len(set(UMIsByCellID[cellID])) for cellID in UMIsByCellID.keys()])
readsPerCell = numpy.array([len(UMIsByCellID[cellID]) for cellID in UMIsByCellID.keys()])
sortedUMIsPerCell, UMIsortedCellIDs = list(reversed([UMIs for UMIs, cellID in sorted(zip(UMIsPerCell, UMIsByCellID.keys()), key=lambda pair: pair[0])])), list(reversed([cellID for UMIs, cellID in sorted(zip(UMIsPerCell, UMIsByCellID.keys()), key=lambda pair: pair[0])]))
sortedReadsPerCell, readSortedCellIDs = list(reversed([reads for reads, cellID in sorted(zip(readsPerCell, UMIsByCellID.keys()), key=lambda pair: pair[0])])), list(reversed([cellID for reads, cellID in sorted(zip(readsPerCell, UMIsByCellID.keys()), key=lambda pair: pair[0])]))

targetValues = {'SRR11683994': 400, 'SRR11683995': 795}
maxValues = {'SRR11683994': 500, 'SRR11683995': 1000}

cumsumUMIs = numpy.cumsum(sortedUMIsPerCell) / numpy.sum(sortedUMIsPerCell)
cumsumReads = numpy.cumsum(sortedReadsPerCell) / numpy.sum(sortedReadsPerCell)

title = 'Control dataset' if fileName == 'SRR11683994' else '4sU dataset'
matplotlib.pyplot.figure()
matplotlib.pyplot.plot(range(2000), 100 * cumsumUMIs[:2000])
matplotlib.pyplot.vlines(targetValues[fileName], 0, 100 * cumsumUMIs[1999], linestyles='dashed') #795 is no. cells specified in the spreadsheet for the treatment (4sU) dataset, which looks to line up nicely with the point of inflection
matplotlib.pyplot.ylabel('Cumulative sum of UMI counts (%)')
matplotlib.pyplot.xlabel('Ordered cell index')
matplotlib.pyplot.title(title)
matplotlib.pyplot.savefig(fileName + '_UMIs.pdf')
matplotlib.pyplot.figure()
matplotlib.pyplot.plot(range(2000), 100 * cumsumReads[:2000])
matplotlib.pyplot.vlines(targetValues[fileName], 0, 100 * cumsumReads[1999], linestyles='dashed')
matplotlib.pyplot.ylabel('Cumulative sum of read counts (%)')
matplotlib.pyplot.xlabel('Ordered cell index')
matplotlib.pyplot.title(title)
matplotlib.pyplot.savefig(fileName + '_reads.pdf')

# probably reads per cell is the better metric at this point since some UMIs will be merged later that have edit / hamming distance <= 1
selectedCells = readSortedCellIDs[:maxValues[fileName]]
selectedReadCounts = sortedReadsPerCell[:maxValues[fileName]]
selectedUMICounts = sortedUMIsPerCell[:maxValues[fileName]]

print(len([cell for cell in selectedCells if cell in [repairedCellIDMapper[ID] for ID in repairedCellIDMapper.keys()]])) #so none of the selected cells had their barcodes repaired anyway for either datasets so no need to worry about it
selectedUMIsByCell = {cell: UMIsByCellID[cell] for cell in selectedCells}
selectedUMIsByReadByCell = {cell: UMIsByReadByByCell[cell] for cell in selectedCells}
pickle.dump(selectedUMIsByReadByCell, open(fileName + '_UMIsByReadByCell.pickle', 'wb'))

