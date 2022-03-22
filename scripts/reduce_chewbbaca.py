#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 08:38:03 2022

@author: chelsea
"""

import pandas as pd
import numpy as np
import sys
import os
import itertools
import pathlib

#name of directory
directory = sys.argv[1]

list_of_files = []

#append filenames into list_of_files
for filename in os.listdir(directory+'/chewbbaca_quast_tables/'):
    list_of_files.append(filename.split('.')[0]) #split on . to separate filename from file extension

combos = [] #list of combinations of filenames

#loop for creating all pairwise combinations of filenames
for L in range(2, 3):
    for subset in itertools.combinations(list_of_files, L):
        combos.append(subset)

#relevant columns
columns = ['Sample','# contigs', 'Largest contig', 'Total length', 'Reference length', 
               'Genome fraction (%)', 'GC (%)', 'Reference GC (%)', 'N50', 'NG50', '# misassemblies',
               '# mismatches per 100 kbp']


#Do pairwise comparison of two dataframes based on filename combinations        
for files in combos:
    filename1 = files[0] #name of first combination of files
    filename2 = files[1] # name of second combination of files
    
    #create dataframes
    df = pd.read_csv(directory+'/chewbbaca_quast_tables/'+filename1+'.tsv', sep='\t')
    df2 = pd.read_csv(directory+'/chewbbaca_quast_tables/'+filename2+'.tsv', sep= '\t')
   
    
    indices = []
    indices2 = []
    
    #add sample IDs of dataframe to list
    for i in df['Sample']:
        new_index = i.split('_') 
        indices2.append(new_index)
    
    #add sample names to list
    for i in range(len(indices2)):
        #a length of 1 means that this is the reference which get apended
        if len(indices2[i]) == 1: 
            indices.append(indices2[i][0])
        #otherwise only add first 3 indices of indices2 list
        else:
            indices.append(indices2[i][0]+'_'+indices2[i][1])
    
    #set elements of indices as index
    df.index = indices 
    df2.index = indices
    
    #extract chewbbaca columns to new dataframes
    df_chewbbaca = df.drop(columns, axis=1)
    df2_chewbbaca = df2.drop(columns, axis=1)
    
    #extract differences between two chewbbaca dataframes
    test2 = df_chewbbaca.compare(df2_chewbbaca, align_axis=1).rename(columns={'self': filename1, 'other': filename2}, level=-1)
    test2 = test2.fillna('')
    
    #extract quast columns to new dataframe
    df_quast = df[columns]
    df2_quast = df2[columns]
    
    #drop sample column
    df_quast = df_quast.drop('Sample', axis=1)
    df2_quast = df2_quast.drop('Sample', axis=1)
    
    #define filenames for csv table
    filename3 = filename1.split('_results')
    filename4 = filename2.split('_results')
 
    #save chewbbaca comparisons
    #test2.to_csv(filename3[0]+'_vs_'+filename4[1]+'chewbbaca.tsv', sep='\t', encoding='utf-8')
    
    #extract differences between two quast dataframes
    test3 = df_quast.compare(df2_quast, align_axis=1).rename(columns={'self': filename, 'other': filename2}, level=-1)
    
    #save quast dataframes
    #test3.to_csv(filename3[0]+'_vs_'+filename4[0]+'_quast.tsv', sep='\t', encoding='utf-8')
    
    #merge quast and chewbbaca results
    #merged = pd.merge(test3,test2,left_index= True, right_index=True, how="outer")
    #merged = merged.fillna('')
    #merged.sort_index(inplace=True, ascending=False)
    
    comp = df.compare(df2, align_axis=1).rename(columns={'self': filename1, 'other': filename2}, level=-1)
    comp = comp.fillna('')
    
    
    #save to csv
    pathlib.Path('Results/Comparisons').mkdir(parents=True, exist_ok=True)
    comp.to_csv('Results/Comparisons/'+ filename3[0]+'_vs_'+filename4[0]+'.tsv', sep='\t', encoding='utf-8')

# #test3 = pd.concat([df_quast, df2_quast], keys=[filename, filename2], axis=1, sort=False)

# cols = list(test3.columns.values)
# length_of_cols = len(cols)

# half_cols = cols[0:int(len(cols)/2)]
# other_half = cols[int(len(cols)/2):len(cols)]

# columns = []

# for i in range(len(half_cols)):
#     columns.append(half_cols[i])
#     columns.append(other_half[i])

# test4 = test3[columns]

#test4.to_csv(filename3[0]+'_vs_'+filename4[0]+'_quast.tsv', sep='\t', encoding='utf-8')

#merged = pd.merge(test4,test2,left_index= True, right_index=True, how="outer")

#merged.to_csv(filename3[0]+'_vs_'+filename4[0]+'_results.tsv', sep='\t', encoding='utf-8')


