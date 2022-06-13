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
import copy


list_of_files = []
#append filenames into list_of_files
for filename in os.listdir('Results/chewbbaca_quast_tables/placeholder/'):
    list_of_files.append(filename.split('_')[0]) #split on . to separate filename from file extension

spades = []
skesa = []
string_lengths = []
first_letters = []

for strings in list_of_files:
    first_letters.append(strings[0]) #Append first character of filenames
    length = len(strings) #Length of each filename
    string_lengths.append(length) #Append lengths
    if 'Ske' in strings:
        skesa.append(strings) #append SKESA assemblies into list skesa
    elif 'Sp' in strings:
        spades.append(strings) #append SPAdes assemblies into list skesa

#skesa_char = skesa[0][1:4] #character denotation for SKESA
spades_char = spades[0][1:3] #character denotation for SPAdes
spades_set = spades[0][3:4] #SPAdes setting denotation

unique_letters = list(set(first_letters)) #remove duplicate first characters
length_of_letters = len(unique_letters) #length of list
unique_lengths = list(set(string_lengths))
length_of_lengths = len(unique_lengths)

#Create empty lists based on how many values are in lists
spades_lists = [[] for _ in range(length_of_lengths)]
skesa_lists = [[] for _ in range(length_of_lengths)]
spades_same = [[] for _ in range(length_of_letters)]
skesa_same = [[] for _ in range(length_of_letters)]

count = 0 #counter

#Append filenames for SKESA and SPAdes in lists
for length in unique_lengths:
    if count == 0:
        for strings in spades:
            if len(strings) == length:
                spades_lists[count].append(strings)
        for strings in skesa:
            if len(strings) == length:
                skesa_lists[count].append(strings)
        count = count + 1
    else:
        for strings in spades:
            if len(strings) == length:
                spades_lists[count].append(strings)
        for strings in skesa:
            if len(strings) == length:
                skesa_lists[count].append(strings)
        count = count+1

count = 0
for letter in unique_letters:
    if count == 0:
        for strings in spades:
            if strings[0] == letter:
                spades_same[count].append(strings)
        for strings in skesa:
            if strings[0] == letter:
                skesa_same[count].append(strings)
        count = count + 1
    else:
        for strings in spades:
            if strings[0] == letter:
                spades_same[count].append(strings)
        for strings in skesa:
            if strings[0] == letter:
                skesa_same[count].append(strings)
        count = count + 1 

skesa_highest = copy.deepcopy(skesa_same)
skesa_lowest = copy.deepcopy(skesa_same)
spades_highest = copy.deepcopy(spades_same)
spades_lowest = copy.deepcopy(spades_same)

len_max = max(string_lengths)
len_min = min(string_lengths)

for i in range(len(skesa_same)):
    shortest_string = [s for s in skesa_highest[i] if len(s) == len_min]
    skesa_highest[i].remove(shortest_string[0])
    longest_string = [s for s in skesa_lowest[i] if len(s) == len_max]
    skesa_lowest[i].remove(longest_string[0])
    
for i in range(len(spades_same)):
    shortest_string = [s for s in spades_highest[i] if len(s) == len_min]
    spades_highest[i].remove(shortest_string[0])
    longest_string = [s for s in spades_lowest[i] if len(s) == len_max]
    spades_lowest[i].remove(longest_string[0]) 

skesa_and_spades = []
for element in spades:
    spades_element = element
    skesa_element = element.replace('Sp', 'Ske')

for element in skesa:
    skesa_element = element
    spades_element = element.replace(skesa_char, spades_char+spades_set)
    duplicate = [skesa_element, spades_element]
    skesa_and_spades.append(duplicate)

all_lists = spades_lists+spades_lowest+spades_highest+skesa_lists+skesa_highest+skesa_lowest+skesa_and_spades

combos = [] #list of combinations of filenames


#loop for creating all pairwise combinations of filenames
for lists in all_lists:
    for L in range(2, 3):
        for subset in itertools.combinations(lists, L):
            combos.append(subset)
to_remove = []

combos2 = list(set(combos))

for tuples in combos2:
    if 'f' in tuples[0] and 'fP' not in tuples[0] and 'P' in tuples[1] and 'fP' not in tuples[1]:
        to_remove.append(tuples)
    if 'f' in tuples[1] and 'fP' not in tuples[1] and 'P' in tuples[0] and 'fP' not in tuples[0]:
        to_remove.append(tuples)



for tuples in to_remove:
    combos2.remove((tuples[0],tuples[1]))    
#relevant columns
columns = ['Sample','# contigs', 'Largest contig', 'Total length', 'Reference length', 
               'Genome fraction (%)', 'GC (%)', 'Reference GC (%)', 'N50', 'NG50', '# misassemblies',
               '# mismatches per 100 kbp']

#Do pairwise comparison of two dataframes based on filename combinations        
for files in combos2:
    filename1 = files[0] #name of first combination of files
    filename2 = files[1] # name of second combination of files
    
    #create dataframes
    df = pd.read_csv('Results/chewbbaca_quast_tables/placeholder/'+filename1+'_results.tsv', sep='\t')
    df2 = pd.read_csv('Results/chewbbaca_quast_tables/placeholder/'+filename2+'_results.tsv', sep='\t')
   
    
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
    df_chewbbaca = df_chewbbaca.applymap(str) #turn dataframe to string
    df2_chewbbaca = df2.drop(columns, axis=1)
    df2_chewbbaca = df2_chewbbaca.applymap(str) #turn dataframe to string

    #extract differences between two chewbbaca dataframes
    test2 = df_chewbbaca.compare(df2_chewbbaca, align_axis=1).rename(columns={'self': filename1, 'other': filename2}, level=-1)
    test2 = test2.fillna('')

    
    #extract quast columns to new dataframe
    df_quast = df[columns]
    df2_quast = df2[columns]
    
    #drop sample column
    df_quast = df_quast.drop('Sample', axis=1)
    df2_quast = df2_quast.drop('Sample', axis=1)
    df_quast = df_quast.applymap(str) #turn dataframe to string
    df2_quast = df2_quast.applymap(str) #turn dataframe to string

    #extract differences between two quast dataframes
    test3 = df_quast.compare(df2_quast, align_axis=1).rename(columns={'self': filename1, 'other': filename2}, level=-1)
    
    #define filenames for csv table
    filename3 = filename1.split('_results')
    filename4 = filename2.split('_results')
 
   
    
    #merge quast and chewbbaca results
    #merged = pd.merge(test3,test2,left_index= True, right_index=True, how="outer")
    #merged = merged.fillna('')
    #merged.sort_index(inplace=True, ascending=False)
    
    df = df.applymap(str) #turn dataframe to string
    df2 = df2.applymap(str) #turn dataframe to string

    comp = df.compare(df2, align_axis=1).rename(columns={'self': filename1, 'other': filename2}, level=-1)
    comp = comp.fillna('')
        
    #save to csv
    pathlib.Path('Results/Comparisons/Chewbbaca_comparisons/').mkdir(parents=True, exist_ok=True)
    #save chewbbaca comparisons
    test2.to_csv('Results/Comparisons/Chewbbaca_comparisons/'+filename3[0]+'_vs_'+filename4[0]+'_chewbbaca.tsv', sep='\t', encoding='utf-8')
    
    
    
    #save quast dataframes
    pathlib.Path('Results/Comparisons/Quast_comparisons/').mkdir(parents=True, exist_ok=True)
    test3.to_csv('Results/Comparisons/Quast_comparisons/'+filename3[0]+'_vs_'+filename4[0]+'_quast.tsv', sep='\t', encoding='utf-8')
    pathlib.Path('Results/Comparisons/Quast+Chewbbaca_comparisons/').mkdir(parents=True, exist_ok=True)
    comp.to_csv('Results/Comparisons/Quast+Chewbbaca_comparisons/'+ filename3[0]+'_vs_'+filename4[0]+'.tsv', sep='\t', encoding='utf-8')



