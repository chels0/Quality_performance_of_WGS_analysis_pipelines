#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 18:42:31 2022

@author: chelsea
"""
import pandas as pd
import numpy as np
import sys
import os
import itertools
import pathlib
import copy

#name of directory
directory = sys.argv[1]
#directory = '/mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Results/chewbbaca_quast_tables/placeholder/'
#directory = 'Results/chewbbaca_quast_tables/placeholder/'

list_of_files = []

#append filenames into list_of_files
for filename in os.listdir(directory):
    list_of_files.append(filename.split('_')[0]) #split on . to separate filename from file extension

#__________________________________________________________________
#                  FINDING PAIRWISE COMPARISONS

assembler_2 = [] #list in which filenames assembled with first assembler are contained
assembler_1 = [] #list in which filesnames assembled with second assembler are contained
string_lengths = [] #list in which length of filenames are contained
first_letters = [] #list of the first characters of filenames are contained, aka trimming options

assembler_1_char = sys.argv[2]
assembler_2_char = sys.argv[3]

#assembler_1_char = 'SpC'
#assembler_2_char = 'SpI'

for strings in list_of_files:
    first_letters.append(strings[0]) #Append first characters of filename 
    length = len(strings)
    string_lengths.append(length) #Append length of filenames
    if assembler_1_char in strings:
        assembler_1.append(strings) #Append first assembly option
    elif assembler_2_char in strings: #append second assembly option
        assembler_2.append(strings)

# assembler_1_char = assembler_1[0][1:4]
# assembler_2_char = assembler_2[0][1:3]
# assembler_2_set = assembler_2[0][3:4]

unique_letters = list(set(first_letters)) #remove duplicates
length_of_letters = len(unique_letters) #amount of trimming options

unique_lengths = list(set(string_lengths)) #remove duplicates
length_of_lengths = len(unique_lengths) #amount of parameters in each filename

assembler_2_trimming_sep = [[] for _ in range(length_of_lengths)] #Create list of empty lists, one for every trimming option
assembler_1_trimming_sep = [[] for _ in range(length_of_lengths)] 
assembler_2_sep_post_improv = [[] for _ in range(length_of_letters)] #Create list of empty lists, one for every post improvement option
assembler_1_sep_post_improv = [[] for _ in range(length_of_letters)]

count = 0 #counter

#Append list of lists with filenames so each list contains different trimming options but the same post improv options
for length in unique_lengths:
    if count == 0:
        for strings in assembler_2:
            if len(strings) == length: #if 
                assembler_2_trimming_sep[count].append(strings)
        for strings in assembler_1:
            if len(strings) == length:
                assembler_1_trimming_sep[count].append(strings)
        count = count + 1
    else:
        for strings in assembler_2:
            if len(strings) == length:
                assembler_2_trimming_sep[count].append(strings)
        for strings in assembler_1:
            if len(strings) == length:
                assembler_1_trimming_sep[count].append(strings)
        count = count+1

count = 0
#Append list of lists with filenames so each list contains different post improv options but the same trimming option
for letter in unique_letters:
    if count == 0:
        for strings in assembler_2:
            if strings[0] == letter:
                assembler_2_sep_post_improv[count].append(strings)
        for strings in assembler_1:
            if strings[0] == letter:
                assembler_1_sep_post_improv[count].append(strings)
        count = count + 1
    else:
        for strings in assembler_2:
            if strings[0] == letter:
                assembler_2_sep_post_improv[count].append(strings)
        for strings in assembler_1:
            if strings[0] == letter:
                assembler_1_sep_post_improv[count].append(strings)
        count = count + 1 

#create deep copies of all list of lists which will contain the longest and shortest filenames respeticvely
assembler_1_highest = copy.deepcopy(assembler_1_sep_post_improv)
assembler_1_lowest = copy.deepcopy(assembler_1_sep_post_improv)
assembler_2_highest = copy.deepcopy(assembler_2_sep_post_improv)
assembler_2_lowest = copy.deepcopy(assembler_2_sep_post_improv)

len_max = max(string_lengths) #length of the filenames with the most software combinations
len_min = min(string_lengths) #length of the filnames with the least software combinations

#Remove filenames which has no post assembly improvement
for i in range(len(assembler_1_sep_post_improv)):
    shortest_string = [s for s in assembler_1_highest[i] if len(s) == len_min] #find shortest filename
    assembler_1_highest[i].remove(shortest_string[0]) #remove list containing shortest filename from highest list
    longest_string = [s for s in assembler_1_lowest[i] if len(s) == len_max] #find longest filename
    assembler_1_lowest[i].remove(longest_string[0]) #remove list containging longest filename from the smallest list
    
for i in range(len(assembler_2_sep_post_improv)):
    shortest_string = [s for s in assembler_2_highest[i] if len(s) == len_min]
    assembler_2_highest[i].remove(shortest_string[0])
    longest_string = [s for s in assembler_2_lowest[i] if len(s) == len_max]
    assembler_2_lowest[i].remove(longest_string[0]) 

assembler_1_and_assembler_2 = []

# for element in assembler_2:
#     assembler_2_element = element
#     assembler_1_element = element.replace(assembler_2_char, assembler_1_char)

#produce a list containing all combinations for each assembly 
for element in assembler_1:
    assembler_1_element = element
    assembler_2_element = element.replace(assembler_1_char, assembler_2_char)
    duplicate = [assembler_1_element, assembler_2_element] #create list of assembler1 and 2
    assembler_1_and_assembler_2.append(duplicate) #append duplicate to list

#combine all lists containig comparisons that should be done in one list
all_lists = assembler_2_trimming_sep+assembler_2_lowest+assembler_2_highest+assembler_1_trimming_sep+assembler_1_highest+assembler_1_lowest+assembler_1_and_assembler_2

combos = [] #list of combinations of filenames

#loop for creating all pairwise combinations of filenames
for lists in all_lists:
    for L in range(2, 3):
        for subset in itertools.combinations(lists, L):
            combos.append(subset)
to_remove = []

combos2 = list(set(combos)) #remove duplicates

#Find pairwise comparisons between Pilon and Filtering and add to a to_remove list
for tuples in combos2:
    if 'f' in tuples[0] and 'fP' not in tuples[0] and 'P' in tuples[1] and 'fP' not in tuples[1]:
        to_remove.append(tuples)
    if 'f' in tuples[1] and 'fP' not in tuples[1] and 'P' in tuples[0] and 'fP' not in tuples[0]:
        to_remove.append(tuples)

#Remove Pilon and Filtering tuples from cominations
for tuples in to_remove:
    combos2.remove((tuples[0],tuples[1]))    
    
#_______________________________________________________________________
#                   PERFORMING PAIRWISE COMPARISONS

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
