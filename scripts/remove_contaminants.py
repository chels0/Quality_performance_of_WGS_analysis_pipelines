#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 13:21:38 2022

@author: chelsea
"""
import sys

fasta_name = sys.argv[1]
contig_size_to_remove = int(sys.argv[2])


with open(fasta_name+'.fasta', "r") as f:
    seq = f.readlines()

every_other = seq[1::2]
to_change = seq.copy()

test = []
    
for i in range(len(every_other)):
    change = i
    change += change + 1
    test.append(change)
    if len(every_other[i]) < contig_size_to_remove:
        to_change[change] = 0
        to_change[change-1] = 0
    
no_int = [x for x in to_change if not isinstance(x,int)]

filtr = str(contig_size_to_remove)

with open(fasta_name+'_filtered'+filtr+'.fasta', 'w') as file2:
    for line in no_int:
        file2.write("%s" % line)
            
    
            
    
