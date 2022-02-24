#!/usr/bin/env bash

#create path to cgMLST

cat Cjejuni_cgMLST_678_listGenes.txt | sed -e "s|^|/data/Raw_data/Campylobacter_jejuni/Schema/Cjejuni_wgMLST_2795_schema/schema_seed_campy_roary_V5/|g" > fullpath_cgMLSTschema.txt

#Allele call

chewBBACA.py AlleleCall -i /home/chelsea/Documents/Quality_performance_of_WGS_analysis_pipelines/Results/all_assemblies/ -g fullpath_cgMLSTschema.txt -o cgMLST --cpu 6
