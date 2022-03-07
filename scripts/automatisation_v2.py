#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 12:52:40 2022

@author: chelsea
"""
import itertools

filtering = ['\tfilter_contigs']
improv = ['\tassembly_improvement']
no_trim = ['\tno_trim']
fastp_trim_qc = ['\tfastp_trim_qc']
trimmomatic = ['\ttrimmomatic']
bools = [' = false\n', ' = true\n']
trues = [' = true\n']
falses = [' = false\n']

no_trim_set = list(itertools.product(no_trim, trues, filtering, bools, improv, bools))
fastp_trim_set = list(itertools.product(fastp_trim_qc, trues, filtering, bools, improv, bools))
trimmomatic_set = list(itertools.product(trimmomatic, trues, filtering, bools, improv, bools))

no_trim_settings = []
fastp_settings = []
trimmomatic_settings = []

for i in range(len(no_trim_set)):
    no_trim_set[i] = list(no_trim_set[i])
    fastp_trim_set[i] = list(fastp_trim_set[i])
    trimmomatic_set[i] = list(trimmomatic_set[i])
    
    no_fastp_trimmomatic = fastp_trim_qc + falses + trimmomatic + falses
    no_no_trim_fastp = no_trim + falses + fastp_trim_qc + falses
    no_no_trim_trimmomatic = no_trim + falses + trimmomatic + falses
    no_trim_set_new = no_trim_set[i]+no_fastp_trimmomatic
    no_trim_settings.append(no_trim_set_new)
    trimmomatic_set_new = trimmomatic_set[i] + no_no_trim_fastp
    trimmomatic_settings.append(trimmomatic_set_new)
    fastp_trim_set_new = fastp_trim_set[i] + no_no_trim_trimmomatic
    fastp_settings.append(fastp_trim_set_new)

count = 0

test_list = []

length = len(no_trim_settings[0])

for i in range(len(no_trim_settings)):
    no_trim_settings[i][0:length] = [''.join(no_trim_settings[i][0:length])]
    fastp_settings[i][0:length] = [''.join(fastp_settings[i][0:length])]
    trimmomatic_settings[i][0:length] = [''.join(trimmomatic_settings[i][0:length])]
    
for element in no_trim_settings:
    count = count + 1
    if count > 9:
    	with open('config_files/Run'+ str(count), 'w') as file2:
        	for line in element:
            		file2.write("%s" % line)
    else:
    	with open('config_files/Run'+'0'+str(count), 'w') as file2:
        	for line in element:
            		file2.write("%s" % line)
    
for element in fastp_settings:
    count = count + 1
    if count > 9:
    	with open('config_files/Run'+ str(count), 'w') as file2:
        	for line in element:
            		file2.write("%s" % line)
    else:
    	with open('config_files/Run'+'0'+str(count), 'w') as file2:
        	for line in element:
            		file2.write("%s" % line)

for element in trimmomatic_settings:                
    count = count + 1
    if count > 9:
    	with open('config_files/Run'+ str(count), 'w') as file2:
        	for line in element:
            		file2.write("%s" % line)
    else:
    	with open('config_files/Run'+'0'+str(count), 'w') as file2:
        	for line in element:
            		file2.write("%s" % line)
    


# for element in no_trim_settings:
#     txt_variables.append(element[0]+' = ' + element[1] + '\n' + element[2]+' = ' + element[3] + '\n' + element[4]+' = ' + element[5] + '\n'+ element[6]+' = ' + element[7] + '\n' + element[8]+' = ' + element[9] + '\n')                       
#     count = count + 1
#     with open('Run'+ str(count), 'w') as file2:
#         for line in txt_variables:
#             file2.write("%s" % line)
#     txt_variables = []
    
# for element in fastp_settings:
#     txt_variables.append(element[0]+' = ' + element[1] + '\n' + element[2]+' = ' + element[3] + '\n' + element[4]+' = ' + element[5] + '\n'+ element[6]+' = ' + element[7] + '\n' + element[8]+' = ' + element[9] + '\n')                                              
#     count = count + 1
#     with open('Run'+ str(count), 'w') as file2:
#         for line in txt_variables:
#             file2.write("%s" % line)
#     txt_variables = []
    
# for element in trimmomatic_settings:
#     txt_variables.append(element[0]+' = ' + element[1] + '\n' + element[2]+' = ' + element[3] + '\n' + element[4]+' = ' + element[5] + '\n'+ element[6]+' = ' + element[7] + '\n' + element[8]+' = ' + element[9] + '\n')                              
#     count = count + 1
#     with open('Run'+ str(count), 'w') as file2:
#         for line in txt_variables:
#             file2.write("%s" % line)
#     txt_variables = []
    
#length = len(fastp_settings[0])
#placeholder = []
#fastp_settings2 = []
#chunk_size = 2

#for i in range(len(fastp_settings)):
#    for j in range(0, length, chunk_size):
#        chunk = fastp_settings[i][j:j+chunk_size]
#        placeholder.append(chunk)
#    fastp_settings2.append(placeholder[:len(fastp_settings)])
    
#for i in range(len(fastp_settings2)):
#    for j in range(len(fastp_settings2[i])):
#        thing = fastp_settings2[i][j][0]+ " = " + fastp_settings2[i][j][1] + "\n"
#        txt_variables.append(thing)
#    count = count + 1
#    with open('Run'+ str(count), 'w') as file2:
#        for line in txt_variables:
#            file2.write("%s" % line)
#    txt_variables = []
