import csv
import scipy
import numpy
import math
from _collections import defaultdict
import matplotlib
import matplotlib.pyplot
import sys
import pickle

#convert pickle files of conversions and UMI counts to txt files so can be transferred to R environment for the inference algorithms
UMIcountByCellByGene = pickle.load(open('SRR11683995_UMIcountBySelectedCellByGene.pickle', 'rb'))
conversionsByCellByGene = pickle.load(open('SRR11683995_conversionsBySelectedCellByGene.pickle', 'rb'))
genes = list(UMIcountByCellByGene.keys())
cells = list(UMIcountByCellByGene[genes[0]].keys())

with open('SRR11683995_data.txt', 'w') as txtFile:
    for gene in genes:
        for cell in cells:
            UMIcount = UMIcountByCellByGene[gene][cell]
            if UMIcount > 0:
                conversions = conversionsByCellByGene[gene][cell]
                conversions = ['/'.join([str(entry[0]), str(entry[1])]) for entry in conversions]
                conversions = '\t'.join(conversions)
                line = gene + '\t' + cell + '\t' + str(UMIcount) + '\t' + conversions + '\n'
            else:
                line = gene + '\t' + cell + '\t' + str(UMIcount) + '\n'
            txtFile.write(line)