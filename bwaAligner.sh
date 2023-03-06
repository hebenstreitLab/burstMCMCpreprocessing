#!/bin/bash

bwa aln GRCh38.p13_genome_index/GRCh38.p13 SRR11683994_2.trimmed.fastq > SRR11683994_2.sai -t 14
bwa samse GRCh38.p13_genome_index/GRCh38.p13 SRR11683994_2.sai SRR11683994_2.trimmed.fastq > SRR11683994_2.sam

bwa aln GRCh38.p13_genome_index/GRCh38.p13 SRR11683995_2.trimmed.fastq > SRR11683995_2.sai -t 14
bwa samse GRCh38.p13_genome_index/GRCh38.p13 SRR11683995_2.sai SRR11683995_2.trimmed.fastq > SRR11683995_2.sam



