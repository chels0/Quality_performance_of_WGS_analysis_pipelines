#!/usr/bin/env bash

#för att kunna skriva config till det man vill

python3 scripts/automatisation_v2.py

pwd > path1.txt
sed 's./.\\/.g' path1.txt > path.txt

path=$(cat path.txt) 

sed "s/PATHS/$path/g" parameter_settings.txt > parameter_path.txt

> log.txt

for file in config_files/*;
do
	head -15 parameter_path.txt > nextflow.config
	cat $file >> nextflow.config
	grep fastqc_set parameter_path.txt -A 10 >> nextflow.config
	
	echo "Starting pipeline with the following paramaters"
	echo $file
	echo | cat $file
	
	echo "Starting pipeline with the following paramaters" >> log.txt 
        echo $file >> log.txt
	echo | cat $file >> log.txt


	source ../../miniconda3/bin/activate env_version1 
	nextflow pipeline6.nf -profile conda -resume
	#multiqc Results/. -o Results/MultiQC
	
	rm -r work
	>nextflow.config 
done
