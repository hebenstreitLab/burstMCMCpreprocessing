# 4sU scRNA-seq data pre-processing pipeline

Provided here are the scripts used for processing and extracting relevant information from the 4sU scRNA-seq data sets used in the paper entitled "Synergising single-cell resolution and 4sU labelling boosts inference of transcriptional bursting". The processing of these Drop-seq data sets was carried out according to the Methods and Supplementary information sections of the paper. This document provides the commands which were used to call the scripts specifically designed for this data set. These commands assume that all of the scripts are contained within the same folder, along with the two pairs raw fasta files, with prefixes SRR11683994 and SRR11683995 representing the control (no 4sU) and 4sU data, respectively, such that we have files named SRR11683994_1.fastq, SRR11683994_2.fastq, SRR11683995_1.fastq and SRR11683995_2.fastq, with the suffixes 1 and 2 representing the pair of barcode and transcript sequence reads, respectively. This data is downloaded from https://www.ebi.ac.uk/ena/browser/view/SRR11683994 and https://www.ebi.ac.uk/ena/browser/view/SRR11683995 and was generated in the Qiu et al 2020 paper entitled "Massively parallel and time-resolved RNA sequencing in single cells with scNT-Seq". The folder must also contain a gtf and fasta file, in this case "gencode.v36.primary_assembly.annotation.gtf" and "GRCh38.primary_assembly.genome.fa", downloaded from https://www.gencodegenes.org/human/. Additionally, it must contain a sub-folder called "GRCh38.p13_genome_index" within which will be found the genome index files built from the fasta file with bwa by executing

``` {eval=F, echo=T}
./buildGenomeIndex.sh
```

One must also execute the following to generate the file "exonsByGeneBySegmentByChromosome.pickle" from the gtf to make association of read alignment positions to genes and exons more efficient.

``` {eval=F, echo=T}
python3 gtfProcessor.py gencode.v36.primary_assembly.annotation.gtf
```

The initial processing follows the "Drop-seq alignment cookbook" (https://mccarrolllab.org/wp-content/uploads/2016/03/Drop-seqAlignmentCookbookv1.2Jan2016.pdf), starting by trimmming read pairs with any base with phred quality ≤ 10 as well as clipping adaptor and polyA tail sequences with following, generating trimmed fastq files named with format "PREFFIX_1/2.trimmed.fastq".

``` {eval=F, echo=T}
./trimMaster.sh
```

Upon execution, the next chunk will repair barcodes with missing bases, discard barcodes with sequencing errors and carry out initial cell selection to discard empty droplets according to the ordered cumulative sum of total reads per cell, which will generate plots with the cut-off values shown according to Qiu et al 2020, along with output pickle data files with naming format "PREFFIX_UMIsByReadByCell.pickle" in which slightly more cells than the cut off are retained for a final cell selection further downstream.

``` {eval=F, echo=T}
python3 processTrimmedBarcodes.py SRR11683994
python3 processTrimmedBarcodes.py SRR11683995
```

The "PREFFIX_2.trimmed.fastq" files with the trimmed transcript sequence reads are then aligned with bwa using the following command to generate corresponding sam files named "PREFFIX_2.sam".

``` {eval=F, echo=T}
./bwaAligner.sh
```

The next command will associate the aligned reads with mapq score >= 10 with their corresponding genes based on their overlap with exons on the same strand, generating "PREFFIX_readInfoByCell.pickle" and "PREFFIX_readInfoByGene.pickle" output files which contains information on the associated gene, cell ID, UMI, read sequence and genomic sequence (sequence in the fasta file at the aligned position) for each read, sorted by cell ID or gene symbol.

``` {eval=F, echo=T}
./readsToGenesMaster.sh
```

The following command then collapses UMIs belonging to the same cell and gene with hamming distance of 1 in order to account for sequencing errors in the barcodes, before extracting the information as required by the inference algorithm; the UMI count for each cell and gene, the number of T>C conversions in each of the aligned reads for each cell and gene (T in the genomic sequence and C in the read sequence) and the total number of genomic Ts in each of those reads. The output files "PREFFIX_UMIcountByCellByGene.pickle" and "PREFFIX_conversionsByCellByGene.pickle" are generated.

``` {eval=F, echo=T}
./sortReadInfoMaster.sh
```

The final cell selection with based on cumulative total reads is now carried out by executing the next command, discarding any cells beyond the defined cut-off values, generating final plots and output files "PREFFIX_UMIcountBySelectedCellByGene.pickle" and "PREFFIX_conversionsBySelectedCellByGene.pickle".

``` {eval=F, echo=T}
./cellSelectionMaster.sh
```

The background T>C rates are then calculated using the control data (SRR11683994) by executing the following command. This will output a file containing the global background T>C rate across all genes and selected cells called "backgroundMutationRate.txt" as well as a file containing the cell-summed, gene-specific counts and rates which are needed for the inference algorithm called "geneSpecificBackgroundMutationRates.txt".

``` {eval=F, echo=T}
python3 backgroundConversionRate.py
```

The previously extracted information for the 4sU data set is then converted into human readable text file called "SRR11683995_data.txt" from the pickle files with the following.

``` {eval=F, echo=T}
python3 pickleToTxt.py
```

An R data structure called "qiuData.Rdata" is then generated from this using the following.

``` {eval=F, echo=T}
Rscript txtToRdata.R
```

We then infer the gene-invariant 4sU-mediated T>C rate (lambda_n) using all genes for which we have high enough confidence in their rates in both the control and 4sU data sets and finally save the data required for the inference algorithm only for selected genes, which are those with at least 1 T>C conversion observed in both the control and 4sU data set, into a file called "qiuDataSelectedGenes.Rdata" using the following command, including UMI counts, T>C counts, genomic T counts, background T>C rate and 4sU-mediated T>C rate. This file is available on Zenodo ().

``` {eval=F, echo=T}
Rscript getSelectedGeneData.R
```

Cell-specific capture efficiencies for the cells in the Qiu 4sU data set were estimated using the following command despite Qiu et al 2020 not making use of spike-in probes. This was carried out in the manner descibed in the supplementary information section of our paper by using cell-matched Drop-seq data with spike-ins from a different paper entitled "Droplet barcoding for single cell transcriptomics applied to embryonic stem cells". Our script makes use of the files "ercc-info.txt" and "GSM1599501_K562_pure_RNA.csv" which can be downloaded from https://github.com/jingshuw/DESCEND_manuscript_source_code/tree/master/spike_in_analysis and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1599501, respectively, and generates a file called "qiuAlphas.Rdata" which contains the estimated capture efficiencies. This file is available on Zenodo ().

``` {eval=F, echo=T}
Rscript ERCCanalysis.R
```