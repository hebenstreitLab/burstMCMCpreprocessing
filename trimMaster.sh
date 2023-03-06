#!/bin/bash
> SRR11683994_1.trimmed.fastq
> SRR11683995_1.trimmed.fastq
> SRR11683994_2.trimmed.fastq
> SRR11683995_2.trimmed.fastq
python3 trimFastqPair.py SRR11683994
python3 trimFastqPair.py SRR11683995

