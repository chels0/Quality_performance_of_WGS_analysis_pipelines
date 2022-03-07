#!/usr/bin/env bash

#create path to cgMLST

path=/home/chelsea/Documents

mkdir Results/chewbbaca/Jejuni

cat Raw_data/Campylobacter_jejuni/Schema/Cjejuni_cgMLST_678_listGenes.txt | sed -e "s|^|${path}/Quality_performance_of_WGS_analysis_pipelines/Raw_data/Campylobacter_jejuni/Schema/Cjejuni_wgMLST_2795_schema/schema_seed_campy_roary_V5/|g" > Results/chewbbaca/Jejuni/fullpath_cgMLSTschema.txt

#Allele call

ls Results/ > directory_list.txt

for dir in $(cat directory_list.txt)
do
	ls Results/${dir}/PT* > Results/${dir}/samples_list.txt
	sed -i '/Fast/d' samples_list.txt #remove fast things
	mkdir Results/${dir}/ChewBBACA
	
	for dir2 in $(cat samples_list.txt)
	do	
		chewBBACA.py AlleleCall -i Results/${dir}/Assemblies/ -g fullpath_cgMLSTschema.txt --cpu 8 -o Results/${dir}/ChewBBACA
	done	
done


