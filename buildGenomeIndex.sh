#!/bin/bash
#downloaded gtf and fasta primary assembly files from here: https://www.gencodegenes.org/human/
cd GRCh38.p13_genome_index/
bwa index -p GRCh38.p13 -a bwtsw ../GRCh38.primary_assembly.genome.fa
