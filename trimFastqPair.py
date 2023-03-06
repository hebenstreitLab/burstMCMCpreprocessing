import csv
import scipy
import numpy
import math
from _collections import defaultdict
import matplotlib
import matplotlib.pyplot
import sys
import pickle

#benefit of trimming pair of files together is that you dont need to do a follow up processing to get the overlap like you would if trimming them seperately

adaptorSequence = 'AAGCAGTGGTATCAACGCAGAGTGAATGGG' #from example in drop-seq alignment cookbook
adaptorSubsequences = [adaptorSequence[i:(i+5)] for i in range(len(adaptorSequence) - 4)]

def findAdaptorClipIndex(adaptorSubsequences, sequence): #function to find where to clip read to remove adaptor sequence
    startBaseIndex = 0
    clipIndex = 0
    for subsequence in adaptorSubsequences:
        if clipIndex > 0 and subsequence != sequence[startBaseIndex:(startBaseIndex + 5)]:
            break
        if subsequence == sequence[startBaseIndex:(startBaseIndex + 5)]:
            startBaseIndex += 1
            clipIndex = startBaseIndex + 4
    return clipIndex

lowPhredSymbols = ['!', '\"', '#', '$', '%', '&', '\'', '(', ')', '*'] #these are the ones below a quality of 10
#lowPhredSymbols = ['!', '\"', '#', '$', '%', '&', '\'', '(', ')', '*', '+', ',', '-', '.']#, '/', '0', '1', '2', '3', '4'] #these are the ones below a quality of 20

fileName = sys.argv[1]
#for the barcode file
fastqPath1 = fileName + '_1.fastq'
trimmedFastqPath1 = fileName + '_1.trimmed.fastq'
#for the DNA sequence file
fastqPath2 = fileName + '_2.fastq'
trimmedFastqPath2 = fileName + '_2.trimmed.fastq'

i = 0
with open(trimmedFastqPath1, 'a') as trimmedFastqFile1, open(fastqPath1) as fastqFile1, open(trimmedFastqPath2, 'a') as trimmedFastqFile2, open(fastqPath2) as fastqFile2:
    for line1, line2 in zip(fastqFile1, fastqFile2):
        if ((i + 4) / 4).is_integer():
            # start new variable for new read and add first bit of content
            readContents1, readContents2 = [line1], [line2]
        elif ((i + 1) / 4).is_integer():
            # add last bit of read content
            readContents1.append(line1)
            readContents2.append(line2)

            # trim read contents of adaptors and polyA tails and append to new fastq file
            # first clip polyA tails by find first occurence of six consecutive As and removing everything from that point on
            sequence1, sequence2 = readContents1[1], readContents2[1]
            quality1, quality2 = readContents1[3], readContents2[3]
            polyAindex = sequence2.find('AAAAAA')
            if polyAindex != -1:
                readContents2[1] = sequence2[:polyAindex] + '\n'
                readContents2[3] = quality2[:polyAindex] + '\n'
            adaptorClipIndex = findAdaptorClipIndex(adaptorSubsequences, sequence2)
            if adaptorClipIndex > 0:
                readContents2[1] = readContents2[1][adaptorClipIndex:]
                readContents2[3] = readContents2[3][adaptorClipIndex:]
            # now only accept reads above a certain length
            readLength = len(readContents2[1])

            tooLowQual = any([symbol in quality1 for symbol in lowPhredSymbols])
            if readLength > 10 and not tooLowQual:  # saying only to accept reads at least 10bp atm after clipping for read 2, since the \n at the end counts as a character have to say <desried min bps>+1, also only high quality barcodes (read 1)
                for trimmedLine in readContents1:
                    trimmedFastqFile1.write(trimmedLine)
                for trimmedLine in readContents2:
                    trimmedFastqFile2.write(trimmedLine)
        # print('polyA:', polyAindex, ' adaptor:', adaptorClipIndex, ' read length:', len(readContents[1]))
        # print(readContents)
        else:
            readContents1.append(line1)
            readContents2.append(line2)
        # if i > 998:
        #     break
        i += 1
fastqFile1.close()
fastqFile2.close()
trimmedFastqFile1.close()
trimmedFastqFile2.close()
