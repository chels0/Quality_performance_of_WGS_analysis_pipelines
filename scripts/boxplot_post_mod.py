#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 11:30:57 2022

@author: chelsea
"""

import pandas as pd
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns


#Script for plotting QUAST metrics, one boxplot for each metric at different coverages
#_____________________________________________________________________________________

directory = 'Results/chewbbaca_quast_tables/placeholder/'

assembler_comp = sys.argv[1]
#assembler_comp = 'Ske'
trimming_software = sys.argv[2]
#trimming_software = 'both'

if trimming_software == 'Pilon':
    trimming_software = 'f'
    full_name_software = 'Pilon'

if trimming_software == 'filtering':
    trimming_software = 'P'
    full_name_software = 'filtering'

elif trimming_software == 'both':
    full_name_software = 'filtering_+_pilon'

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
    for cov in coverages:
        count = count + 1 #increase counter
        #empty lists
        only_trimming = []
        list_of_trimming_column_vales = []
        list_of_software_combinations = []
        list_of_column_values_lists = []
        
        #extract data from files
        for filename in list_of_files:
            split = filename.split('_')[0]
            if trimming_software == 'both':
                if 'P' in split or 'f' in split:
                    df = pd.read_csv(directory+filename, sep='\t', index_col=0) #create tab separated dataframe
                    df = df[columns] #keep only relevant columns
                    df = df.iloc[1: , :] #remove reference row
                    df = df[df.index.str.contains(cov)] #keep only samples of certain coverage
                    df_assembler_1olumn_values_list = df[col].to_list() #turn dataframe column into list
                    list_of_software_combinations.append(filename.split('_')[0]) #append software combination in list
                    list_of_column_values_lists.append(df_assembler_1olumn_values_list) #append dataframe in list form to list
                    
                    only_trimming.append(filename.split('_')[0]) #append software combinations without post assembly improvements
                    list_of_trimming_column_vales.append(df_assembler_1olumn_values_list) #append column values of files with no post assembly improvements
                    
            elif trimming_software not in filename and len(split)<9:
                df = pd.read_csv(directory+filename, sep='\t', index_col=0) #create tab separated dataframe
                df = df[columns] #keep only relevant columns
                df = df.iloc[1: , :] #remove reference row
                df = df[df.index.str.contains(cov)] #keep only samples of certain coverage
                df_assembler_1olumn_values_list = df[col].to_list() #turn dataframe column into list
                list_of_software_combinations.append(filename.split('_')[0]) #append software combination in list
                list_of_column_values_lists.append(df_assembler_1olumn_values_list) #append dataframe in list form to list
                
                only_trimming.append(filename.split('_')[0]) #append software combinations without post assembly improvements
                list_of_trimming_column_vales.append(df_assembler_1olumn_values_list) #append column values of files with no post assembly improvements
                    
            
        #transform dataframe for plotting on the bassis of assembler   
        df_2 = pd.DataFrame(list_of_trimming_column_vales) #create dataframe
        df_2 = df_2.transpose() #transpose
        df_2.columns = only_trimming #set software combinations without psot assembly improvements as columns 
        
        
        test= df_2.filter(like='Ske', axis=1)
        
        df_assembler_1 = df_2.filter(like=assembler_1, axis=1) #create dataframe with only values from first assembler
        df_assembler_2 = df_2.filter(like=assembler_2, axis=1) # create dataframe with only values from second assembler
        
        df_assembler_1_melt = df_assembler_1.melt().assign(assembler=assembler_1_full) #unpivot dataset for plotting
        df_assembler_2_melt = df_assembler_2.melt().assign(assembler=assembler_2_full) #unpivot dataset for plotting
        
        
        concated = pd.concat([df_assembler_1_melt, df_assembler_2_melt]) #combine dataframes
        
        #Replace denotations to full names for plotting
        concated['variable'] = concated['variable'].str.replace('SpI','')
        concated['variable'] = concated['variable'].str.replace('SpC','')
        concated['variable'] = concated['variable'].str.replace('Ske','')
        concated['variable'] = concated['variable'].str.replace('F','Fastp ')
        concated['variable'] = concated['variable'].str.replace('N','No trimming ')
        concated['variable'] = concated['variable'].str.replace('T','Trimmomatic ')
        concated['variable'] = concated['variable'].str.replace('P','Pilon ')
        concated['variable'] = concated['variable'].str.replace('f','filtering ')
        concated['variable'] = concated['variable'].str.replace('fP','filtering Pilon')
        concated = concated.rename(columns = {'variable': 'Trimming software', 'value': col })
        
        #set figure labels        
        if col == '# contigs':
            fig_nr = 'A'
            #range_ =
        if col == 'Largest contig':
            fig_nr = 'B'
        if col == 'N50':
            fig_nr = 'B'
        if col == 'Genome fraction (%)':
            fig_nr = 'C'
        if col == '# misassemblies':
            fig_nr = 'D'
            
        if count == 6: #if statement to print legend on only the 100x figure
            plt.figure()
            sns.set_theme()
            sns.set(rc={"figure.figsize":(7, 5), "figure.dpi":300, 'savefig.dpi':300})
            ax = sns.boxplot(data=concated, x='assembler', y=col, hue='Trimming software', palette="Paired")
            ax.set_title(fig_nr+cov+': '+col+' at coverage '+cov, fontsize = 15)
            plt.tight_layout()
            plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0.)
        
        else:
            plt.figure()
            sns.set_theme()
            sns.set(rc={"figure.figsize":(7, 5), "figure.dpi":300, 'savefig.dpi':300})
            ax = sns.boxplot(data=concated, x='assembler', y=col, hue='Trimming software', palette="Paired")
            ax.set_title(fig_nr+cov+': '+col+' at coverage '+cov, fontsize = 15)
            plt.legend([],[], frameon=False)
            plt.tight_layout()

        plt.savefig('Results/Conclusions/'+full_name_software+'_'+col+'_'+cov+'_.png', bbox_inches="tight")
        
        
        
