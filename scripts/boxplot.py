#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 14:20:12 2022

@author: chelsea
"""

import pandas as pd
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns


#Script for plotting QUAST metrics, one boxplot for each metric at different coverages
#run with command: python3 path_to_script/boxplot.py flag1 flag2
#where flag1 is the assembler compared to SPADes --isolate. Put Ske for SKESa or SpC for --careful
#flag2 is the trimming option compared to Trimmomatic, put no_trim to plot no trimming compared to Trimmomatic or Fastp to compare Fastp to Trimmomatic
#_____________________________________________________________________________________


directory = 'Results/chewbbaca_quast_tables/placeholder/'

assembler_comp = sys.argv[1]
trimming_software = sys.argv[2]



#Define which assemblers have been used
assembler_1 = assembler_comp
assembler_2 = 'SpI'
#Define what name should be plotted
if assembler_1 == 'Ske':
    assembler_1_full = 'SKESA'
if assembler_1 == 'SpC':
    assembler_1_full = 'SPAdes --careful'
assembler_2_full = 'SPAdes --isolate'

only_trimming = []
list_of_trimming_column_vales = []
list_of_software_combinations = []
list_of_files = []
list_of_column_values_lists = []

#columns to plot
columns = ['# contigs', 'Largest contig', 'Total length', 
               'Genome fraction (%)', 'GC (%)', 'N50', 'NG50', '# misassemblies']


columns = ['# contigs', 
               'Genome fraction (%)', 'N50', '# misassemblies']

#coverages
coverages = ['20x', '50x', '100x']
empty_df = pd.DataFrame({}) #create empty dataframe

#add filenames to list
for filename in os.listdir(directory):
    list_of_files.append(filename)

list_of_files.sort() #sort list


count=0 #counter
#plot all quast metrics (columns) for each coverage
for col in columns:
    df = pd.read_csv(directory+filename, sep='\t', index_col=0) #create tab separated dataframe
    df = df[col]
    largest_value = []
    smallest_value = []
    for cov in coverages:
        only_cov_df = df[df.index.str.contains(cov)]
        only_cov_df.astype(int)
        max_value = only_cov_df.to_numpy().max()
        min_value = only_cov_df.to_numpy().min()
        largest_value.append(max_value)
        smallest_value.append(min_value)
    
    y_axis_max = max(largest_value)  
    y_axis_min = min(smallest_value)
    if y_axis_min == 0:
        y_axis_min = -1
    
    
    for cov in coverages:
        count = count + 1 #increase counter
        #empty lists
        only_trimming = []
        list_of_trimming_column_vales = []
        list_of_software_combinations = []
        list_of_column_values_lists = []
        
        #extract data from files
        for filename in list_of_files:
            df = pd.read_csv(directory+filename, sep='\t', index_col=0) #create tab separated dataframe
            df = df[columns] #keep only relevant columns
            df = df.iloc[1: , :] #remove reference row
            df = df[df.index.str.contains(cov)] #keep only samples of certain coverage
            df_assembler_1olumn_values_list = df[col].to_list() #turn dataframe column into list
            list_of_software_combinations.append(filename.split('_')[0]) #append software combination in list
            list_of_column_values_lists.append(df_assembler_1olumn_values_list) #append dataframe in list form to list
            
            if 'P' not in filename and 'f' not in filename:
                only_trimming.append(filename.split('_')[0]) #append software combinations without post assembly improvements
                list_of_trimming_column_vales.append(df_assembler_1olumn_values_list) #append column values of files with no post assembly improvements
                
            
        #transform dataframe for plotting on the bassis of assembler   
        df_2 = pd.DataFrame(list_of_trimming_column_vales) #create dataframe
        df_2 = df_2.transpose() #transpose
        df_2.columns = only_trimming #set software combinations without psot assembly improvements as columns 
        df_assembler_1 = df_2[['N'+assembler_1, 'T'+assembler_1]] #create dataframe with only values from first assembler
        df_assembler_2 = df_2[['N'+assembler_2, 'T'+assembler_2 ]] # create dataframe with only values from second assembler
        
        df_assembler_1_melt = df_assembler_1.melt().assign(assembler=assembler_1_full) #unpivot dataset for plotting
        df_assembler_2_melt = df_assembler_2.melt().assign(assembler=assembler_2_full) #unpivot dataset for plotting
        
        df_assembler_1_f = df_2[['F'+assembler_1, 'T'+assembler_1]] #create dataframe with only values from first assembler
        df_assembler_2_f = df_2[['F'+assembler_2, 'T'+assembler_2 ]] # create dataframe with only values from second assembler
        
        df_assembler_1_melt_f = df_assembler_1_f.melt().assign(assembler=assembler_1_full) #unpivot dataset for plotting
        df_assembler_2_melt_f = df_assembler_2_f.melt().assign(assembler=assembler_2_full) #unpivot dataset for plotting
        
        
        #concated = pd.concat([df_assembler_1_melt, df_assembler_2_melt]) #combine dataframes
        if trimming_software == 'Fastp' or trimming_software == 'fastp' or trimming_software == 'F' :
            concated = pd.concat([df_assembler_1_melt_f, df_assembler_2_melt_f]) #combine dataframes
        else:
            concated = pd.concat([df_assembler_1_melt, df_assembler_2_melt])
        
        #Replace denotations to full names for plotting
        concated['variable'] = concated['variable'].str.replace('SpI','')
        concated['variable'] = concated['variable'].str.replace('SpC','')
        concated['variable'] = concated['variable'].str.replace('Ske','')
        concated['variable'] = concated['variable'].str.replace('N','No trim '+cov)
        concated['variable'] = concated['variable'].str.replace('F','Fastp '+cov)
        concated['variable'] = concated['variable'].str.replace('T','Trimmomatic '+cov)
        concated = concated.rename(columns = {'variable': 'Trimming software', 'value': col })
        
        #set figure labels        
        if col == '# contigs':
            fig_nr = 'A'
        if col == 'Largest contig':
            fig_nr = 'B'
        if col == 'N50':
            fig_nr = 'B'
        if col == 'Genome fraction (%)':
            fig_nr = 'C'
        if col == '# misassemblies':
            fig_nr = 'D'
            
        if count == 9: #if statement to print legend on only the 100x figure
            plt.figure()
            sns.set_theme()
            sns.set(rc={"figure.figsize":(7, 5), "figure.dpi":300, 'savefig.dpi':300})
            ax = sns.boxplot(data=concated, x='assembler', y=col, hue='Trimming software', palette="vlag")
            ax = sns.swarmplot(data=concated, x='assembler', y=col, hue='Trimming software', dodge = True, s=2.5)
            ax.set_title(fig_nr+cov+': '+col+' at coverage '+cov, fontsize = 15)
            plt.tight_layout()
            plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0.)
        
        else:
            plt.figure()
            sns.set_theme()
            sns.set(rc={"figure.figsize":(7, 5), "figure.dpi":300, 'savefig.dpi':300})
            ax = sns.boxplot(data=concated, x='assembler', y=col, hue='Trimming software', palette="vlag").set(title=cov)
            plt.legend([],[], frameon=False)
            ax = sns.swarmplot(data=concated, x='assembler', y=col, hue='Trimming software', dodge = True, s=2.5)
            ax.set_title(fig_nr+cov+': '+col+' at coverage '+cov, fontsize = 15)
            plt.legend([],[], frameon=False)
            plt.tight_layout()

        plt.savefig('Results/Conclusions/'+trimming_software+'_'+col+'_'+cov+'_.png', bbox_inches="tight")
        
        
