#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 16:41:35 2022

@author: chelsea
"""
import pandas as pd
import numpy as np
import os
import itertools
import pathlib


directory = 'Results/Comparisons/Chewbbaca_comparisons'
#directory = '/mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Results/Comparisons/Chewbbaca_comparisons/'

pathlib.Path('Results/Conclusions/').mkdir(parents=True, exist_ok=True) #Create path if non existant

#Column names
diff = '# Differences'
corr = '# Corrections'
wrong = '# Errors'
change = '# Changes'


list_of_files = []
list_of_files1 = []
list_of_files2 = []

#append filenames into list_of_files
for filename in os.listdir(directory):
    split_files = filename.split('_') #split filename on _
    tupls = (split_files, filename) #add split filename and filename in tuple
    list_of_files.append(tupls) #append tupl to list
    list_of_files1.append(split_files[0]) #append first software combination in pairwise comparison in list
    list_of_files2.append(split_files[2]) #append second
    
all_lists = list_of_files1 + list_of_files2 #concatenate lists

files = []
word_count = []
word_count2 = []
df_coverages = []

characters = ['P','500f','500fP', '200f', '200fP'] #post improvement characters

#for each post improvement, count the number of differences observed between using the post assembly improvement and not using it, and output in a table for each software combination
for char in characters:
    if 'final_result' in locals():
        del final_result #remove the local variable final_result if it exists to not append results
    
    files = []
    #append filenames containing post assembly improvement software
    for filename in all_lists:
        if char in filename:
            files.append(filename)

    no_dup = list(set(files)) #remove duplicates
    
    not_char = []
    #Remove the post assembly improvement character denotation from software combination
    #and append to list
    for value in no_dup:
        combination_wo_char = value.replace(char, '') #replace character
        tupls = (value, combination_wo_char) #tuple original software combination with replace software combination to compare the two with each other later
        not_char.append(tupls) #append to list
    
    
    #create template dataframes
    template = {diff: [], corr: [], wrong: [], change: []} #create template dictionary
    template_columns = [diff, corr, wrong, change] #columns for dataframe
    x20x = pd.MultiIndex.from_product([['20x'], template_columns], names=["Coverage", ""]) #create a 20x template multicolumn
    x50x = pd.MultiIndex.from_product([['50x'], template_columns], names=["Coverage", ""])#create a 50x template multicolumn
    x100x = pd.MultiIndex.from_product([['100x'], template_columns], names=["Coverage", ""])#create a 100x template multicolumn
    df_20x_template = pd.DataFrame(template, columns = x20x) #template df
    df_50x_template = pd.DataFrame(template, columns = x50x)
    df_100x_template = pd.DataFrame(template, columns = x100x)
    
    list_of_columns = []
    list_of_columns1 = []
    list_of_columns2 = []
    word_list = []
    word_list2 = []
    
    #fill dataframe with differences
    for tuple_ in not_char:
        for value in list_of_files:
            if ((tuple_[0] == value[0][0] and tuple_[1] == value[0][2]) 
                or (tuple_[0] == value[0][2] and tuple_[1] == value[0][0])): 
                df = pd.read_csv(directory+'/'+value[1], index_col = 0, sep='\t') #Using if statement to load dataframe of pairwise comparison between software combination containing character and not containing it
                df.index = df.index.fillna('No label') #set header as no label
                #separate assemblies on coverage into several dataframes and append to list
                df_20x = df[df.index.str.contains('20x|No label', regex = True)] #keep 20x assemblies
                df_50x = df[df.index.str.contains('50x|No label', regex = True)] #keep 50x assemblies
                df_100x = df[df.index.str.contains('100x|No label', regex = True)] # keep 100x assemblies
                dfs = []
                dfs.append((df_20x, '20x')) #append dataframe along with the coverage
                dfs.append((df_50x, '50x' ))
                dfs.append((df_100x, '100x'))                
                
                #Loop through all dataframes with different coverages
                for element in dfs:
                    dataframe = element[0] #the dataframe
                    coverage  = element[1] #the coverage of the dataframe
                    list_of_columns = []
                    list_of_columns1 = []
                    list_of_columns2 = []
                    for col in dataframe.columns:
                        list_of_columns.append(col) #append dataframe columns to dataframe
                    
                    #Separate duplicated loci into two lists
                    for i in range(len(list_of_columns)):
                        if i % 2 == 0:
                            list_of_columns1.append(list_of_columns[i])
                        else:
                            list_of_columns2.append(list_of_columns[i])
                    
                    #counter to count differences
                    count_right = 0
                    count_wrong = 0
                    count_change = 0
                    count_diffs = 0
                    
                    #clear lists
                    word_list = []
                    word_list2 = []
                    word_count = []
                    word_count2 = []
                    
                    if not dataframe.empty:
                        if char in dataframe.iloc[0,0]: #if post assembly character is on the left side of header continue
                    
                            #loop through columns and count differences. A zero is the same as reference. Non-zeroes are errors from the reference
                             for i in range(len(list_of_columns1)):
                                word_list = []
                                df_copy = dataframe.copy() #copy dataframe
                                df_copy[wrong]  = np.where((dataframe[list_of_columns1[i]] != '0') & 
                                                           (dataframe[list_of_columns2[i]] == '0'), 
                                                           dataframe[list_of_columns1[i]], np.nan) #find assemblies where post assembly improvement created an error from not using the improvement            
                                count_wrong = count_wrong + len(df_copy[wrong].dropna()) #add the amount of errors to the error counter
                                
                                no_top_row = df_copy.iloc[1:, :] #remove header row
                                mask = (no_top_row[wrong].str.extract(('((?!^\d+$)+(?!^\W)^.+$)'), expand=False)) #find all string errors in error column
                                new_mask = mask.dropna(axis=0) #extract string errors
                                if len(new_mask) != 0:
                                    word_list.append(new_mask[0]) #if there are string errors append them to list
                                uniques = list(set(word_list)) #remove duplicates
                                
                                for words in uniques:
                                    counter = word_list.count(words) #count how many of the same error exists in list
                                    word_count.append((words, counter)) #append the error and how many times it shows
                                    
                                df_copy = dataframe.copy() #copy dataframe
                                df_copy[diff]  = np.where((dataframe[list_of_columns1[i]] != dataframe[list_of_columns2[i]]), 
                                                          dataframe[list_of_columns1[i]], np.nan) #count total amount of differences           
                                count_diffs = count_diffs + len(df_copy.iloc[1:][diff].dropna()) #add amount of differences to counter, not counting header row
                                
                                df_copy = dataframe.copy()
                                df_copy[corr]  = np.where((dataframe[list_of_columns1[i]] == '0') & 
                                                          (dataframe[list_of_columns2[i]] != '0'), 
                                                          dataframe[list_of_columns1[i]], np.nan) #count corrections           
                                count_right = count_right + len(df_copy[corr].dropna()) #add amount of corrections to counter
                                df_copy = dataframe.copy()
                                df_copy[change]  = np.where((dataframe[list_of_columns1[i]] != '0') & (dataframe[list_of_columns2[i]] != '0') & (dataframe[list_of_columns1[i]] != dataframe[list_of_columns2[i]]), dataframe[list_of_columns1[i]], np.nan)            
                                count_change = count_change + len(df_copy.iloc[1:][change].dropna())
                                
                        
                        if char in dataframe.iloc[0,1]: #if post assembly character is on the right side of header
                            for i in range(len(list_of_columns1)):
                                word_list2 = []
                                df_copy = dataframe.copy()
                                df_copy[wrong]  = np.where((dataframe[list_of_columns2[i]] != '0') & (dataframe[list_of_columns1[i]] == '0'), dataframe[list_of_columns2[i]], np.nan)            
                                count_wrong = count_wrong + len(df_copy[wrong].dropna())
                                no_top_row = df_copy.iloc[1:, :]
                                
                                mask = (no_top_row[wrong].str.extract(('((?!^\d+$)+(?!^\W)^.+$)'), expand=False))
                                new_mask = mask.dropna(axis=0)
                                if len(new_mask) != 0:
                                    word_list2.append(new_mask[0])
                                
                                uniques2 = list(set(word_list2))
                                for words in uniques2:
                                    counter = word_list2.count(words)
                                    word_count2.append((words, counter))
                                
                                df_copy = dataframe.copy()
                                df_copy[diff]  = np.where((dataframe[list_of_columns2[i]] != dataframe[list_of_columns1[i]]), dataframe[list_of_columns2[i]], np.nan)            
                                count_diffs = count_diffs + len(df_copy.iloc[1:][diff].dropna())
                                
                                df_copy = dataframe.copy()
                                df_copy[corr]  = np.where((dataframe[list_of_columns2[i]] == '0') & (dataframe[list_of_columns1[i]] != '0'), dataframe[list_of_columns2[i]], np.nan)            
                                count_right = count_right + len(df_copy[corr].dropna())
                                
                                df_copy = dataframe.copy()
                                df_copy[change]  = np.where((dataframe[list_of_columns2[i]] != '0') & (dataframe[list_of_columns1[i]] != '0') & (dataframe[list_of_columns1[i]] != dataframe[list_of_columns2[i]]), dataframe[list_of_columns2[i]], np.nan)            
                                count_change = count_change + len(df_copy.iloc[1:][change].dropna())
                                   
                                
                        post_assembly_software = tuple_[0]
                        not_post_assembly_software = tuple_[1]
                        correct_result = {(coverage,diff) : count_diffs, (coverage,corr) : count_right, 
                                          (coverage, wrong) : count_wrong, (coverage, change) : count_change}  #dictionary with results                       
                        df_list_of_columns = pd.DataFrame(correct_result, 
                                                          index = [post_assembly_software+ ' vs ' + not_post_assembly_software]) #dataframe of correct results
                        
                        uniques = list(set(word_count)) #remove duplicates
                        uniques2 = list(set(word_count2))#remove duplicates
                        uniques3 = list(set(uniques + uniques2)) #concat lists and remove duplicates
                        
                        #Series of if statements to calculate the amount of string errors 
                        if word_count: #if there are string errors do the following
                            for unique_value in uniques:
                                df_list_of_columns[coverage, '% ' + unique_value[0]] = 0 #set string error as column and give it value 0 in dataframe
                                
                            for words in word_count:
                                df_list_of_columns[coverage, '% ' + words[0]] = (words[1] + (df_list_of_columns[coverage, '% ' + words[0]][0])) #add the number of string error to the existing number of string errors in the dataframe
        
                        if word_count2 and not word_count:
                            for unique_value in uniques2:
                                df_list_of_columns[coverage, '% ' + unique_value[0]] = 0
                            for words2 in word_count2:
                                df_list_of_columns[coverage, '% ' + words2[0]] = (words2[1] + (df_list_of_columns[coverage, '% ' + words2[0]][0]))
                            
                        if word_count2 and word_count:
                            for words in word_count2:
                                df_list_of_columns[coverage, '% ' + words[0]] = ((words[1] + df_list_of_columns[coverage, words[0]][0]))
                        
                        #calculate the amount of wrong called alleles, aka errors which are not string errors if there are no string errors at all
                        for i in range(len(df_list_of_columns)):
                            if df_list_of_columns[coverage, wrong][i] != 0: #if there are errors in the dataframe
                                if not word_count2 and not word_count: #and if there are no string errors
                                    df_list_of_columns[coverage, '% Wrong-called'] = count_wrong #wrong called is the amount of errors counted
                        
                        column = int(df_list_of_columns[coverage, wrong][0]) #convert error column values to int
                        to_subtract_list = []
                        
                        #Calculate the amount of wrong called alleles, aka errors that are not string errors
                        for value in uniques3:
                            to_subtract = int(df_list_of_columns[coverage, '% ' + value[0]][0]) #how many string errors
                            to_subtract_list.append(to_subtract) #append to list
                            sum_ = sum(to_subtract_list) #sum the amount of string errors
                            wrong_called = column - sum_ #subtract the amount of string errors from the total amount of errors to get wrong called errors
                            
                            df_list_of_columns[coverage, '% Wrong-called'] = (wrong_called) #add wrong called
                        
                        #Convert what the percentages of errors are made up of
                        for col in df_list_of_columns:
                            if '%' in col[1]:
                                wrong_col = df_list_of_columns[coverage, wrong][0] #error column value
                                df_list_of_columns[col] = df_list_of_columns[col][0]*(100/wrong_col) #calculate percentage and add to dataframe
                                
                        #concatenate template with dataframe
                        if coverage == '20x':
                            df_20x_template= pd.concat([df_20x_template, df_list_of_columns], sort = False)
                        elif coverage == '50x':
                            df_50x_template = pd.concat([df_50x_template, df_list_of_columns], sort = False)
                        else:
                            df_100x_template = pd.concat([df_100x_template, df_list_of_columns], sort = False)
                                        
        #merge together the different coverage dataframes into one dataframe
        final_result = df_20x_template.merge(df_50x_template, left_index = True, right_index = True)
        final_result = final_result.merge(df_100x_template, left_index = True, right_index = True)
       
        #add post assembly improvement software combinations which did not differ from combinations without the improvement as an empty row
        if len(final_result) < len(no_dup):
             empty = {('20x', diff): [0], ('50x', diff): [0], ('100x', diff): [0]} #empty dictionary if final_result table is shorter than the amount of software combinations with post assembly improvements
             list_of_not_in = []
             for tuple_ in not_char:
                 if (tuple_[0]+' vs '+tuple_[1]) not in final_result.index:
                     df_empty = pd.DataFrame(data=empty, index=[tuple_[0]+' vs '+tuple_[1]]) #create dataframe of empty dictionary and add comparison as index
                     list_of_not_in.append(df_empty) #append to list
             
             df_all_empty = pd.concat(list_of_not_in) #concatenate all dataframes in list to one dataframe
             
             final_result = pd.concat([final_result, df_all_empty]) #concat empty dataframes with final_result
        final_result = final_result.replace('nan', 0) #replace nan with 0
        final_result = final_result.fillna(0) #replace nan with 0
        final_result = final_result.astype(int) #convert to int
        final_result = final_result.sort_index(axis=0) #sort index
        final_result.to_csv('Results/Conclusions/'+char+'results.tsv', sep='\t', encoding='utf-8')

