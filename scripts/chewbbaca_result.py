#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import sys
import os
import pathlib


#Script for comparing the reference genome's allele calling results to the assemblies' allele calling results.
#The script subtracts the reference's type for each loci from all the assemblies' loci. 
#Thus the reference will have a row of zeroes representing a correct result
#Assemblies with a type of 0 at loci is a correct allele calling

#_______________________________________________________________________________________________________________

directory = sys.argv[1]


#Define relevant columns to keep in dataframe
columns = ['# contigs', 'Largest contig', 'Total length', 'Reference length', 
               'Genome fraction (%)', 'GC (%)', 'Reference GC (%)', 'N50', 'NG50', '# misassemblies',
               '# mismatches per 100 kbp']

#Loop through all files in cgMLST result
for direct in os.listdir(directory):
    filename_chew = directory + '/' + direct + '/chewBBACA/cgMLST_results_jejuni/results_alleles.tsv'
   
    df = pd.read_csv(filename_chew, sep='\t', index_col=0) #create tab separated dataframe
    df = df.applymap(str) #change dataframe to string
    
    df.index.name = 'Sample' #set Sample as index name of dataframe
    df.index = df.index.str.replace('.fasta','', regex=True) #remove .fasta from sample names
    
    reference = df.iloc[[0],:] #get reference row, the first row of dataframe with all its columns
    reference = reference[[]] #empty the reference row
    string_errors = []
    
    #Find all letter characters in each column and add column to list containing either character or nan
    mask = ([df[col].str.extract(('((?!^\d+$)^.+$)'), expand=False) for col in df])
    #Append letter characters to list
    for element in mask:
        mask2 = element.dropna().drop_duplicates() #drop all nans and duplicates
        for col in mask2:
            string_errors.append(col) #append letters from each column
  
    
    string_errors_no_dup = set(string_errors) #remove duplicate letters
    string_errors_no_dup = list(string_errors_no_dup) #create list of letters instead of set
    count = list(range(10000, (len(string_errors_no_dup)+1)*10000, 10000)) #create list of ints ranging from 10000 to the length of letter list times 10000
    
    string_error_int_tuple = list(zip(string_errors_no_dup, count)) #Give each string error an int ID from count 
    
    #Replace all string errors in dataframe with their IDs in order to be able to perform subtractions on entire dataframe
    for values in string_error_int_tuple:
        df = df.replace(values[0], values[1], regex=True)
    
    df = df.applymap(int) #convert dataframe to int
    
    to_compare = df.iloc[[0],:] #reference row that all other rows will be compared to
    
    empty = []
    pathlib.Path("Results/chewbbaca_quast_tables").mkdir(parents=True, exist_ok=True)
    
    #Compare each row in dataframe to reference row and calculate difference. Also count how many loci differ from reference
    for i in range(len(df)):
        nextt = df.iloc[[0,i],:] #row to be compared to reference along with reference
	#Count how many loci differ from reference at each row
        ref = nextt.iloc[[0]].reset_index() #reference without index
        ref.drop('Sample', axis=1, inplace=True) #remove sample column from reference
        to_change = nextt.iloc[[1]] #row to be compared to reference
        ref.insert(0, 'Sample', to_change.index) #insert new sample column with index of row to be compared to reference
        ref.set_index('Sample', inplace = True) # set sample as index
        comp = ref.compare(to_change, align_axis=1).rename(columns={'self': 'Reference', 'other': to_change.index[0]}, level=-1) #compare row and reference 
        colu = comp.columns.get_level_values(0) #extract how many differences at each loci (column) for that row
        colu = colu[~colu.duplicated()] #remove duplicates
        
        #print how many different columns there are compared to reference in txt file
        with open('Results/chewbbaca_quast_tables/'+direct+'_samples.txt', 'a') as file2:
            file2.write("%s" % to_change.index[0] + '\t' + str(len(colu)) + '\t')
            file2.writelines(colu + ' ')
            file2.write('\n')
        
        #subtract reference value from each row one at the time
        nextt = nextt.diff() #calculate difference between rows
        to_change = nextt.iloc[[1],:] #extract row compared to reference
        subtracted_row = to_change #set variable as row
        to_compare = pd.concat([to_compare, subtracted_row]) #append row to to_compare
    
    #Convert ints back to their letter characters
    for col in to_compare:
        to_compare[col] =to_compare[col].astype(int) #had to set as int because diff function outputs floats
        for values in string_error_int_tuple:
            cond = (df[col]<= values[1]) & (df[col]> (values[1]-2000)) #set condition if value in dataframe is in range of a string errors ID
            to_compare.loc[cond, col] = values[0] #replace int with proper string error
        
    to_compare = to_compare[~to_compare.index.duplicated(keep='last')] #remove duplicate reference row
    
    #COLUMNS
    columns_with_no_diff = to_compare.apply(pd.value_counts).dropna(thresh=10) #columns where a lot of samples have the same values
    
    #QUAST
    filename_quast = directory + '/' + direct + '/MultiQC/multiqc_quast.tsv'
    #filename_quast = sys.argv[2] #filename for quast result
    df2 = pd.read_csv(filename_quast, sep='\t', index_col=0) #create tab separated dataframe
    df4 = df2[columns] #keep only relevant columns
    df3 = pd.concat([reference, df4], sort=False) #add empty reference row to dataframe
    df3 = df3.reset_index() #remove sample names as index

    df3.index = df.index #set index of quast dataframe as the index of chewbbaca dataframe
    df3 = df3.drop('Sample', axis=1) #remove sample column which is not index
    
    df3 = df3.applymap(str) #turn dataframe to string
    
    #MERGER
    df_quast_chew = pd.merge(df3, to_compare, how ="right", left_index= True, right_index = True) #merge two dataframes 
    df_quast_chew.fillna('', inplace=True) #fill nans with empty string
    df_quast_chew.replace('nan', '', inplace=True) #replace string nan with empty string
    name_of_run = df.index[1] #sample name
    name = name_of_run.split('x_') #split sample name
   
    
    df_quast_chew.to_csv('Results/chewbbaca_quast_tables/'+ direct+'_results.tsv', sep='\t', encoding='utf-8') #save to csv with name of run
    columns_with_no_diff.to_csv('Results/chewbbaca_quast_tables/'+ direct+'_common_columns.tsv', sep='\t', encoding='utf-8') #save to csv with name of run

    
