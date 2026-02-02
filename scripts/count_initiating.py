#!/usr/bin/env python
import sys

#   Part of the RiboMethSeq analysis pipeline - count initating reads
#   USAGE:
#       python3 count_initiating.py [.genome file] [good pairs bed file]

#   how to use
if len(sys.argv) < 3:
    exit("\nPart of the RiboMethSeq analysis pipeline - count initating reads:\n\t"
         "python3 count_initiating.py [.genome file] [good pairs bed file]\n")

genome_dict = {}  # [chrom] -> length

# genome file
with open(sys.argv[1], 'r') as genome_file:
    for line in genome_file:
        chrom, length = line.strip().split('\t')
        genome_dict[chrom] = [0] * (int(length))

# bed file
with open(sys.argv[2], 'r') as bed_file:
    for line in bed_file:
        chrom, start, end, location, quality, strand = line.strip().split('\t')
        if int(start) < len(genome_dict[chrom]):
            genome_dict[chrom][int(start)] += 1

for keys, array in genome_dict.items():
    count = 0
    for index, value in enumerate(array):
        print(keys, count, value, sep='\t')
        count += 1
