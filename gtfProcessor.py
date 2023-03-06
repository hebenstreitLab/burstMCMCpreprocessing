import csv
import scipy
import numpy
import math
from _collections import defaultdict
import matplotlib
import matplotlib.pyplot
import sys
import pickle

#for sorting data in the gtf file into an exon by gene by chromosome structure to allow alignments to be mapped to genes
#including all genes including scaffolds, mitochondrial ones etc then can narrow down to just protein-coding genes on the main chromosomes later if desired
fileName = sys.argv[1]
gtfData = list(csv.reader(open(fileName, 'r'), delimiter='\t')) #this is 1-indexed
gtfData = gtfData[5:]
chromosomes = list(set(line[0] for line in gtfData))

#need to make a dictionary with genomic segmentation to hugely speed up the downstream mapping of reads to exons/genes
segmentSize = 1000000
strandByGene = defaultdict(str)
exonsByGeneBySegmentByChromosome = {chromosome: {index: defaultdict(list) for index in range(300)} for chromosome in chromosomes}
for line in gtfData:
    if line[2] == 'gene':
        gene = line[8].split(';')[2].split(' ')[2].split('"')[1]
        geneType = line[8].split(';')[1].split(' ')[2].split('"')[1]
        chromosome, strand = line[0], line[6]
        validGene = True #chromosome in chromosomes and geneType == 'protein_coding' # can use to get only protein coding genes
        if validGene:
            strandByGene[gene] = strand
            #find the segment/s it corresponds to
            positions = [int(line[3]), int(line[4])]
            segments = set([math.floor(position / segmentSize) for position in positions])
            for segment in segments: #in case a gene spans multiple segments
                exonsByGeneBySegmentByChromosome[chromosome][segment][gene] = []
    if line[2] == 'exon' and validGene:
        positions = [int(line[3]), int(line[4])]
        segments = set([math.floor(position / segmentSize) for position in positions])
        for segment in segments: #in case an exon spans multiple segments
            exonsByGeneBySegmentByChromosome[chromosome][segment][gene].append(positions)

pickle.dump(strandByGene, open('strandByGene.pickle', 'wb'))
pickle.dump(exonsByGeneBySegmentByChromosome, open('exonsByGeneBySegmentByChromosome.pickle', 'wb'))

