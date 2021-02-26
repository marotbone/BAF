#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 14:46:31 2020

@author: maro
"""
"""
Takes as sys.argv[1] a reference FASTA file and as sys.argv[2] a corresponding .bam file.
Counts the total reads aligned per contig sequence and returns json file.
sys.argv[3]=output filename

"""
import sys
import json
import pysam
from Bio import SeqIO
ref_names=[]
counts={}

fasta_path=sys.argv[1]
bam_path=sys.argv[2]
target_filename=sys.argv[3]

for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
    ref_names.append(seq_record.id)

samfile = pysam.AlignmentFile(sys.argv[2], "rb")
for i in ref_names:
    if samfile.count(i) > 0:
        print(samfile.count(i),i)
        counts[i]=samfile.count(i)
#print(counts)
sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)

with open('{}'.format(target_filename), 'w') as outfile:
    json.dump(sorted_counts,outfile)
samfile.close()
