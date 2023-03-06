
import csv
import scipy
import numpy
import math
from _collections import defaultdict
import matplotlib
import matplotlib.pyplot
import sys
import pickle
from Bio import SeqIO

#go through all reads, check if read ID belongs to a valid cell, if the mapq score is high enough (10 minimum),
#then determine if it maps to a an exon of a gene, if yes record as valid read with genomic position (for checking total Ts for the TC rates)

fileName = sys.argv[1] #'SRR11683995'
filePath = fileName + '_2.sam'

targetValues = {'SRR11683994': 400, 'SRR11683995': 795}

UMIsByReadByCell = pickle.load(open(fileName + '_UMIsByReadByCell.pickle', 'rb'))

#get cells and UMIs by read ID
barcodesByRead = defaultdict(tuple)
for cell in UMIsByReadByCell.keys():
    for read in UMIsByReadByCell[cell].keys():
        barcodesByRead[int(read)] = (cell, UMIsByReadByCell[cell][read])

segmentSize = 1000000
exonsByGeneBySegmentByChromosome = pickle.load(open('exonsByGeneBySegmentByChromosome.pickle', 'rb'))
chromosomes = list(exonsByGeneBySegmentByChromosome.keys())
strandByGene = pickle.load(open('strandByGene.pickle', 'rb'))

def alignmentToGeneMapper(positions, chromosome, segments, exonsByGeneBySegmentByChromosome): #find which gene/s exon the read maps to, if any
    genes = []
    for segment in segments:
        for gene in exonsByGeneBySegmentByChromosome[chromosome][segment].keys():
            for exon in exonsByGeneBySegmentByChromosome[chromosome][segment][gene]:
                if min(positions[1], exon[1]) > max(positions[0], exon[0]):
                    genes.append(gene)
    return list(set(genes))


compliment = {'N': 'N', 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
def reverseCompliment(sequence, compliment):
    return(''.join(list(reversed([compliment[base] for base in sequence]))))



#get the fasta genome sequences so can add genome sequence straight to record
fastaSequences = defaultdict()
for fastaSequence in SeqIO.parse('GRCh38.primary_assembly.genome.fa', 'fasta'):
    fastaSequences[fastaSequence.id] = fastaSequence.seq


#make a data structure to track which reads belong to which gene readsByCellByGene = defaultdict(dict)
readInfoByGene = defaultdict(dict)
counter = 0
for chromosome in chromosomes:
    for segment in exonsByGeneBySegmentByChromosome[chromosome].keys():
        for gene in exonsByGeneBySegmentByChromosome[chromosome][segment].keys():
            readInfoByGene[gene] = defaultdict(list)

#and another to track reads by cell
readInfoByCell = {cell: defaultdict(list) for cell in UMIsByReadByCell.keys()}

#sam positions from bwa are 1-indexed, just like those in the gtf file
i = 0
with open(filePath) as samFile:
    for line in samFile:
        if i > 194: #skip the headers
            alignment = line.split('\t')
            mapqScore = int(alignment[4])
            if mapqScore >= 10:
                readID = int(alignment[0].split('.')[1])
                if len(barcodesByRead[readID]) > 0:
                    chromosome, direction = alignment[2], alignment[1]
                    if chromosome in chromosomes and direction in ['0', '16']: #potentially can cut out scaffolds and mitochondria at this stage, ensure we have a direction for the read
                        strand = '+' if direction == '0' else '-'
                        start = int(alignment[3])
                        sequence = alignment[9]
                        readLength = len(sequence)
                        end = start + readLength - 1
                        positions = [start, end]
                        segments = set([math.floor(position / segmentSize) for position in positions])
                        geneWithExonicMap = alignmentToGeneMapper(positions, chromosome, segments, exonsByGeneBySegmentByChromosome)
                        geneWithExonicMap = [gene for gene in geneWithExonicMap if strandByGene[gene] == strand] #only consider genes with the same direction as the read
                        if len(geneWithExonicMap) > 0:  # find exons the read maps to
                            cell, UMI = barcodesByRead[readID]
                            genomicSequence = fastaSequences[chromosome][(start - 1):end]
                            if strand == '-':
                                sequence, genomicSequence = reverseCompliment(sequence, compliment), reverseCompliment(genomicSequence, compliment)
                            readInfo = [geneWithExonicMap, cell, UMI, sequence, str(genomicSequence)]
                            readInfoByCell[cell][readID] = readInfo
                            for gene in geneWithExonicMap:
                                readInfoByGene[gene][readID] = readInfo
        if i % 1000000 == 0:
            print(i)
        i += 1

samFile.close()
#stop RAM becoming an issue
del(UMIsByReadByCell)
del(barcodesByRead)
del(exonsByGeneBySegmentByChromosome)
del(fastaSequences)
pickle.dump(readInfoByCell, open(fileName + '_readInfoByCell.pickle', 'wb'))
pickle.dump(readInfoByGene, open(fileName + '_readInfoByGene.pickle', 'wb'))
