#!/usr/bin/env bash

#create path to cgMLST

path=/mnt/bigdisk

mkdir Results/chewbbaca/Jejuni

cat Raw_data/Campylobacter_jejuni/Schema/Cjejuni_cgMLST_678_listGenes.txt | sed -e "s|^|${path}/Quality_performance_of_WGS_analysis_pipelines/Raw_data/Campylobacter_jejuni/Schema/Cjejuni_wgMLST_2795_schema/schema_seed_campy_roary_V5/|g" > Results/chewbbaca/Jejuni/fullpath_cgMLSTschema.txt

#Allele call

chewBBACA.py AlleleCall -i ${path}/Quality_performance_of_WGS_analysis_pipelines/Results/all_assemblies/ -g Results/chewbbaca/Jejuni/fullpath_cgMLSTschema.txt -o cgMLST --cpu 2
