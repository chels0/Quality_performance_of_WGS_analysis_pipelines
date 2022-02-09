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
hej = len(every_other[0])
to_change = seq.copy()
    
for i in range(len(every_other)):
    if len(every_other[i]) < contig_size_to_remove:
        to_change[i+2] = 0
        to_change[i+1] = 0
    
no_int = [x for x in to_change if not isinstance(x,int)]

with open(fasta_name+'_polished.fasta', 'w') as file2:
    for line in no_int:
        file2.write("%s" % line)
            
    
            
    
