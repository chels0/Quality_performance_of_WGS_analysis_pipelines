#!/usr/bin/env bash

#för att kunna skriva config till det man vill

python3 scripts/automatisation_v2.py

pwd > path1.txt
sed 's./.\\/.g' path1.txt > path.txt

path=$(cat path.txt) 

sed "s/PATHS/$path/g" parameter_settings.txt > parameter_path.txt

> log.txt

source ../../miniconda3/bin/activate env_version1 
if [[ -z $(ls -A Results) ]]
then
	for file in config_files/*;
	do
		head -15 parameter_path.txt > nextflow.config
		cat $file >> nextflow.config
		grep fastqc_set parameter_path.txt -A 10 >> nextflow.config
	
		date
		echo "Starting pipeline with the following paramaters"
		echo $file
		echo | cat $file
	
		date >> log.txt
		echo "Starting pipeline with the following paramaters" >> log.txt 
        	echo $file >> log.txt
		echo | cat $file >> log.txt
		
		nextflow pipeline6.nf -profile conda -resume
		#multiqc Results/. -o Results/MultiQC
	
		rm -r work
		>nextflow.config 
	done
	bash scripts/sort_assemblies_v3.sh
	
	#ls Results > directory_list.txt
	path2=$(pwd)
	cat Raw_data/Campylobacter_jejuni/Schema/Cjejuni_cgMLST_678_listGenes.txt | sed -e "s|^|${path2}/Raw_data/Campylobacter_jejuni/Schema/Cjejuni_wgMLST_2795_schema/schema_seed_campy_roary_V5/|g" > fullpath_cgMLSTschema.txt
	for dir in Results/*
	do
		chewBBACA.py AlleleCall -i ${dir}/Assemblies/ -g fullpath_cgMLSTschema.txt --cpu 8 -o ${dir}/chewBBACA/cgMLST_results_jejuni
	done        
else
    	#ls Results > directory_list.txt
	path2=$(pwd)
	cat Raw_data/Campylobacter_jejuni/Schema/Cjejuni_cgMLST_678_listGenes.txt | sed -e "s|^|${path2}/Raw_data/Campylobacter_jejuni/Schema/Cjejuni_wgMLST_2795_schema/schema_seed_campy_roary_V5/|g" > fullpath_cgMLSTschema.txt
	for dir in Results/*
	do
		chewBBACA.py AlleleCall -i ${dir}/Assemblies/ -g fullpath_cgMLSTschema.txt --cpu 8 -o ${dir}/chewBBACA/cgMLST_results_jejuni
	done 
fi

rm path1.txt
rm path.txt
rm directory_list.txt




