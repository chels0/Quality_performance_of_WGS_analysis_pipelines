#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 11:27:15 2022

@author: chelsea
"""

import re
import sys

fasta_name = sys.argv[1]
contig_size_to_remove = int(sys.argv[2])

with open(fasta_name+'.fasta', "r") as f:
    seq = f.readlines()

lista = seq.copy()

indices = []
bajs = []
contigs = []
empty = []


r = re.compile(">")
newlist = list(filter(r.match, lista))

for element in newlist:
    if element in lista:
       bajs.append(lista.index(element))

for element in newlist:
    contigs.append(element.split('_'))
              
       
       
for i in range(len(contigs)):
    contigs[i].append(bajs[i])
    if int(contigs[i][3]) < contig_size_to_remove:       
        empty.append(int(contigs[i][6]))

del lista[empty[0]:]

with open(fasta_name+'_filtered.fasta', 'w') as file2:
    for line in lista:
        file2.write("%s" % line)
               
        
    
