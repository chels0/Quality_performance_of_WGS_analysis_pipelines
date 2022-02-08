#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 11:05:08 2022

@author: chelsea
"""

#Script for calculating which quality threshold trimmomatic should have for its trailing command line option. The script uses a cut off of 28 as a good quality score and everything below will be cut.

#Import pandas and numpy
import pandas as pd
import numpy as np

#Read csv file containing quality scores
df_files= pd.read_csv('xx03.csv', sep=' ')

#splits dataset into two parts
splt = np.array_split(df_files, 2)

#create list with the lower values from the split
small_end = splt[1]
small_end = small_end[small_end['Mean'] < 28] #keep only values cells with smaller mean quality score than 28

to_cut_end = small_end['#Base'].tolist() #add the bases with low quality score to list

list_of_values= [] #empty list to add values to

#Iterate through list of bases with low quality score and split bases on - and append values to list
for value in to_cut_end:
    to_cut_end = value.split('-')
    list_of_values.append(to_cut_end[0])
    list_of_values.append(to_cut_end[1])    

list_end = [] #empty list to add values to

#Iterate through list of bases and append those values to empty list    
for i in range(len(list_of_values)):
    list_end.append(int(list_of_values[i]))

#Assign the trailing parameter
if not list_end: #If there are no values smaller than 28, trailing is 1
    trail = 1
else:
    min_end = min(list_end)
    max_end = max(list_end)
    trail = max_end - min_end 

#If the amount of bases needed to be cut is bigger than 20, allow quality to be moderate
if trail > 20:
    trail = 20
else:
    trail = 28 #otherwise get the best quality
    
print(trail) #define output for nextflow

    
