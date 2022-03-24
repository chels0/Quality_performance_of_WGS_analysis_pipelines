#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import sys
import os
import pathlib

list_of_directories = []

directory = sys.argv[1]

#Define relevant columns to keep in dataframe
columns = ['# contigs', 'Largest contig', 'Total length', 'Reference length', 
               'Genome fraction (%)', 'GC (%)', 'Reference GC (%)', 'N50', 'NG50', '# misassemblies',
               '# mismatches per 100 kbp']



for direct in os.listdir(directory):
    filename_chew = directory + '/' + direct + '/chewBBACA/cgMLST_results_jejuni/results_alleles.tsv'

    #filename_chew = direct + '/chewBBACA/cgMLST_results_jejuni/results_alleles.tsv'
    

    #filename_chew = '/mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Results/alleles/results_alleles1.tsv'
    
    #filename_chew = sys.argv[1] #filename for chewbbaca result
    df = pd.read_csv(filename_chew, sep='\t', index_col=0) #create tab separated dataframe
    df = df.applymap(str) #change dataframe to string
    
    df.index.name = 'Sample' #set Sample as index name of dataframe
    df.index = df.index.str.replace('.fasta','', regex=True) #remove .fasta from sample names
    
    # indices = df.index.str.split('_') #split sample names on underscores
    
    # ids = [] #empty list for sample IDs
    # #Add sample IDs to list
    # for values in indices:
    #     ids.(values[0]) #append sample ID
    
    reference = df.iloc[[0],:] #get reference row, the first row of dataframe with all its columns
    reference = reference[[]] #empty the reference row
    tom = []
    hej = []
    
    #Find all letter characters in each column and add column to list containing either character or nan
    mask = ([df[col].str.extract(('((?!^\d+$)^.+$)'), expand=False) for col in df])
    
    #Append letter characters to list
    for element in mask:
        mask2 = element.dropna().drop_duplicates() #drop all nans and duplicates
        for col in mask2:
            hej.append(col) #append letters from each column
    
    test = set(hej) #remove duplicate letters
    test = list(test) #create list of letters instead of set
    
    count = list(range(10000, (len(test)+1)*10000, 10000)) #create list of ints ranging from 10000 to the length of letter list times 10000
    
    ko = list(zip(test, count)) #tuple letter character with an int from count
    
    #Replace all letter characters in dataframe with ints
    for values in ko:
        df = df.replace(values[0], values[1], regex=True)
    
    df = df.applymap(int) #convert dataframe to int
    
    to_compare = df.iloc[[0],:] #reference row that all other rows will be compared to
    
    empty = []
    
    #Compare each row in dataframe to reference row and calculate difference
    for i in range(len(df)):
        nextt = df.iloc[[0,i],:] #row to be compared to reference along with reference
        nextt = nextt.diff() #calculate difference between rows
        to_change = nextt.iloc[[1],:] #extract row compared to reference
        bajs = to_change #set variable as row
        to_compare = pd.concat([to_compare, bajs]) #append row to to_compare
    
    #Convert ints back to their letter characters
    for col in to_compare:
        to_compare[col] =to_compare[col].astype(int) #had to set as int because diff function outputs floats
        for values in ko:
            cond = (df[col]<= values[1]) & (df[col]> (values[1]-2000)) #set condition if value in dataframe is same int as letter blabla
            to_compare.loc[cond, col] = values[0] #replace int with proper letter character
    
        
    to_compare = to_compare[~to_compare.index.duplicated(keep='last')] #remove duplicate reference row
    
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
    kuk = pd.merge(df3, to_compare, how ="right", left_index= True, right_index = True) #merge two dataframes 
    kuk.fillna('', inplace=True) #fill nans with empty string
    kuk.replace('nan', '', inplace=True) #replace string nan with empty string
    name_of_run = df.index[1] #sample name
    name = name_of_run.split('x_') #split sample name
    #name_of_run = name_of_run.split('_contigs')
    
    #if len(name_of_run) == 1:
    #    name_of_run = name_of_run[0].split('_scaffolds')
   
    pathlib.Path("Results/chewbbaca_quast_tables").mkdir(parents=True, exist_ok=True)
    
    kuk.to_csv('Results/chewbbaca_quast_tables/'+ direct+'_results.tsv', sep='\t', encoding='utf-8') #save to csv with name of run


# for index, row in df.iterrows():
#     hej.append(row[0])
#     print(row[0])


# for col in df:
#     hejda = df[col].value_counts()
#     maximum = hejda.max()
#     test2 = hejda[hejda < maximum]
#     new = test2.index
#     hej.append(test2)
    
# bajs = []
# piss = []
# for element in hej:
#     grej = element.reset_index()
#     for kuk in element:
#         bajs.append(kuk)
#     bajs = []
#     piss.append(bajs)
    
