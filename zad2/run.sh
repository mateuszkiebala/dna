#!/bin/bash
python debruijn_algorithm.py reads/reads3.fasta contigs3.fasta
bowtie2 -a --local --mp 2,2 --rdg 10,2 --rfg 10,2 -f -x reference/reference -U contigs3.fasta > bow.out
python evaluate.py bow.out