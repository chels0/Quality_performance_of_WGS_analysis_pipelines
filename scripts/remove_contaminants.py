#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 13:21:38 2022

@author: chelsea
"""
import sys

with open(sys.argv[1], "r") as f:

    seq = f.readlines()
    every_other = seq[1::2]
    hej = len(every_other[0])
    to_change = seq.copy()
    
    for i in range(len(every_other)):
        if len(every_other[i]) < 1000:
            to_change[i+2] = 0
            to_change[i+1] = 0
    
    no_int = [x for x in to_change if not isinstance(x,int)]


            
    
            
    
