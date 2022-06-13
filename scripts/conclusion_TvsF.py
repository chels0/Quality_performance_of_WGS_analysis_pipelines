#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 12:24:15 2022

@author: chelsea
"""
import pandas as pd
import matplotlib.pyplot as plt

directory = 'Results/Conclusions'

df_T = pd.read_csv(directory+'/N+T_results.tsv' , header = [0,1],index_col = 0, sep='\t')
df_F = pd.read_csv(directory+'/N+F_results.tsv' , header = [0,1], index_col = 0, sep='\t')


diff = '# Differences'
corr = '# Corrections'
wrong = '# Errors'
change = '# Changes'
coverages = ['20x', '50x', '100x']

df_T_index = []
df_F_index = []

new_index_F = []
new_index_T = []

new_index_tot = []

for index in df_T.index:
    to_use = index.split(' ')
    df_T_index.append(to_use[len(to_use)-1])

for index in df_F.index:
    to_use = index.split(' ')
    df_F_index.append(to_use[len(to_use)-1])
    
for i in range(len(df_T_index)):
    first_char = df_T_index[i][0]
    word = df_T_index[i].replace(first_char, '')
    index_in_F = df_F_index.index(df_F_index[0][0]+word)
    new_index_T.append((df_T_index[i], i))
    new_index_F.append((df_F_index[0][0]+word, index_in_F))
    #new_index_tot.append(df_T_index[i] + ' vs ' + df_F_index[0][0]+word)


# sorted_ = [None]*len(new_index_F)

# for i in range(len(new_index_F)):
#     index = new_index_F[i][1]
#     element = new_index_F[i][0]
#     sorted_.insert(index, element)

# sorted_ = list(filter(None, sorted_))

temp = df_F.iloc[[1]]
temp = temp.iloc[1: ,:]                

#change order of dataframe F to math T
for value in new_index_F:
    f_row = df_F.iloc[[value[1]]]
    bajs = f_row
    temp = pd.concat([temp, bajs])

temp2 = temp.iloc[[value[1]]]
temp2 = temp2.iloc[1:, :]
temp3 = temp.iloc[[value[1]]]
temp3 = temp3.iloc[1:, :]

for i in range(len(df_T)):
    t_row = df_T.iloc[[i]]      #row to be compared to reference along with reference
    f_row = temp.iloc[[i]]
    bajs = pd.concat([t_row, f_row])
    temp3 = pd.concat([bajs, temp3])
    bajs = bajs.diff()
    test = bajs.iloc[[1]]
    temp2 = pd.concat([temp2, test])
    
temp2.index = temp2.index.str.replace('N', 'T')
counter3 = 0
for coverage in ['20x', '50x', '100x']:
    counter3 = counter3 + 1
    if counter3 < 3:
        ax = temp2[[(coverage, diff), (coverage, corr), (coverage, wrong), (coverage, change)]].plot(kind='bar', title=coverage, legend = False)
        ax.set_ylim(-200, 200)

    else:
        ax = temp2[[(coverage, diff), (coverage, corr), (coverage, wrong), (coverage, change)]].plot(kind='bar', title=coverage)
        ax.set_ylim(-200, 200)
        ax.legend(labels = [diff, corr, wrong, change], title='Labels', loc='center left',bbox_to_anchor=(1, 0.5))
    plt.savefig('Results/Conclusions/F+T_'+coverage+'_plot.png', bbox_inches='tight')
temp2.to_csv('Results/Conclusions/F+T_results.tsv', sep='\t', encoding='utf-8')

to_plot5 = temp3[[('20x', diff), ('50x', diff), ('100x', diff)]]
to_plot6 = temp3[[('20x', corr), ('50x', corr), ('100x', corr)]]
to_plot7 = temp3[[('20x', wrong), ('50x', wrong), ('100x', wrong)]]
to_plot8 = temp3[[('20x', change), ('50x', change), ('100x', change)]]

col_length = len(to_plot5.columns)
for j in range(col_length):
    for i in range(len(to_plot5)):
        diffe = to_plot5.iloc[i,j]
        to_plot6.iloc[i,j] = to_plot6.iloc[i,j]*(100/diffe)
        to_plot7.iloc[i,j] = to_plot7.iloc[i,j]*(100/diffe)
        to_plot8.iloc[i,j] = to_plot8.iloc[i,j]*(100/diffe)

fig, axes = plt.subplots(nrows=2, ncols=2, figsize = (20,20))
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=1)
to_plot5.plot(kind='bar', title=' % Differences', ax = axes[0,0])
to_plot6.plot(kind='bar', title='% Corrections', ax = axes[0,1])
to_plot7.plot(kind='bar', title='% Errors', ax = axes[1,0])
to_plot8.plot(kind='bar', title='% Change', ax = axes[1,1])
plt.savefig('Results/Conclusions/TF_plot_coverage.png')
