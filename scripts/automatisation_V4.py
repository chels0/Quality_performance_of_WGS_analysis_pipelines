#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 16:43:29 2022

@author: chelsea
"""
#Import modules
import itertools
import re

#Define parameters to be used in lists
filtering = ['\tfilter_contigs']
improv = ['\tassembly_improvement']
no_trim = ['\tno_trim']
fastp_trim_qc = ['\tfastp_trim_qc']
trimmomatic = ['\ttrimmomatic']
assembler = ['\tassembler']
assembly_programs = [" = 'skesa'\n", "= 'spades'\n"]
bools = [' = false\n', ' = true\n']
trues = [' = true\n']
falses = [' = false\n']

#Put your own dependant variables here
filter_set = ['\tfilter_set']
filter_variables = [' = 200\n', ' = 500\n' ]
#filter_variables = [' = 200\n']
spades_set = ['\tspades_set']
spades_variables = [' = --careful\n', ' = --isolate\n']
#spades_variables = [' = --careful\n']
#Put them in a list here with the negative dependancy
dependant_variables = [(filtering[0]+falses[0], filter_set[0]) , (assembler[0]+assembly_programs[0], spades_set[0])]

#Create all possible combinations when no trimming
testtt = list(itertools.product(no_trim, trues, fastp_trim_qc, falses, trimmomatic, falses, filtering, bools, improv, bools,assembler, assembly_programs))
testtt2 = list(itertools.product(filter_set, filter_variables, spades_set, spades_variables))
no_trim_set = list(itertools.product(testtt, testtt2))


#Create all combinations of parameters when fastp trimming is chosen
testtt = list(itertools.product(fastp_trim_qc, trues, no_trim, falses, trimmomatic, falses, filtering, bools, improv, bools,assembler, assembly_programs))
fastp_set = list(itertools.product(testtt, testtt2))

#Create all combinations of parameters when trimmomatic is chosen
testtt = list(itertools.product(trimmomatic, trues, no_trim, falses, fastp_trim_qc, falses, filtering, bools, improv, bools,assembler, assembly_programs))
trimmomatic_set = list(itertools.product(testtt, testtt2))

#List with all the software
all_diff = [no_trim_set, fastp_set, trimmomatic_set]

#Counter for text creator later
count2 = 0

#For each software in software list, create txt files containing all the combinations
#of possible software. 
for software in all_diff:
    no_trim_list = []
    test = []
    #For each start software combination, concat each parameter along with its settings
    #and add to no_trim_list
    for val in software:
        val = list(val) #turn tuple to list
        no_trim_list = [] #clear lists
        bajs = [] 
        
        for i in range(len(val[0])): #for each non dependant value in software list
            val[0] = list(val[0]) #non-dependant variables, turn tuple to list
            val[1] = list(val[1]) #dependant variables, turn tuple to list
            if i % 2: #find which indices are not dividable by 2 to be able to concat values with indices 0 and 2, 1 and 3 and so on
                start = i-1 #start index is before index not divided by 2
                other = i #oter index is not divided by 2
                valll = val[0][start] + val[0][other] #concat non-dependable parameters with its settings
                no_trim_list.append(valll) #add concated parameter to list
        length = len(no_trim_list) #length of list
        no_trim_list[0:length] = [''.join(no_trim_list[0:length])] #merge all elements in list to one element
        
        #same as above but for dependable variables
        for i in range(len(val[1])):
            if i % 2:
                start = i-1
                other = i
                vall2 = val[1][start] + val[1][other]
                bajs.append(vall2)
        length2 = len(bajs)
        bajs[0:length2] = [''.join(bajs[0:length2])]
        
        test.append(no_trim_list[0]+bajs[0]) #join together lists of dependable and non-dependable variables
    
    kuken = []
    kuken2 = []
    test3 = []
    
    #Remove dependant variables from element in list where the parameter the variable is dependant on does not exist
    for element in test:
        count = 0 #counter
        for variables in dependant_variables: 
            if variables[0] in element: #if parameter variables are dependant on is not in element of list
                count = count + 1 #counter goes up if if fulfilled
                element = re.sub(variables[1]+'.*\n', '',element) #remove dependant variable 
                changed = element #new element
                
        if count == 0: #if counter is zero, all parameters the dependant variables are dependant on are in this element and the element is not changed
            kuken2.append(element) #append element
        
        kuken.append(changed) #append changed element
        
    no_dup = list(set(kuken)) #remove duplicates
    no_dup_2 = list(set(kuken2)) #remove dupliaces
    total = no_dup+no_dup_2 #the full list of variables
    
    #Create txt files for each combination
    for element in total:
        count2 = count2 + 1 #counter for id of config file
    
        #If count is less than 9 add a 0 before count
        if count2 > 9:
            with open('config_files/Run'+ str(count2), 'w') as file2:
                for line in element:
                    file2.write("%s" % line)
        else:
             with open('config_files/Run'+'0'+str(count2), 'w') as file2:
                 for line in element:
                    file2.write("%s" % line)

