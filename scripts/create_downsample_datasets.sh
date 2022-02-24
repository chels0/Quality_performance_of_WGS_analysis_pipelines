#!/usr/bin/bash

# Script for downsampling reads

#Define coverage values to downsample to
coverages=( 20 50 100 )

#Create directories for each coverage value
mkdir Raw_data/${coverages[0]}x_coverage
mkdir Raw_data/${coverages[1]}x_coverage
mkdir Raw_data/${coverages[2]}x_coverage

#Loop through each of the samples and coverage values to calculate the coverages using Rasusa.
for sample in {1..4} #loops through samples
do
	ls Raw_data/PT28-${sample}* > Raw_data/path.txt #create file containing all samples
	sed 's/Raw_data\///g' Raw_data/path.txt > Raw_data/sample.txt #remove Raw_data from file names 
	sed 's/_.*//g' Raw_data/sample.txt > Raw_data/prefix_dupl.txt #remove everything after sample prefix
	uniq Raw_data/prefix_dupl.txt Raw_data/prefix.txt #keep only single values and not duplicates

	for prefix in $(cat Raw_data/prefix.txt) #loop through file with prefixes 
	do
		for coverage_value in ${coverages[@]} #loop through coverages
		do
		rasusa -i Raw_data/${prefix}_R1.fastq.gz -i Raw_data/${prefix}_R2.fastq.gz --coverage ${coverage_value} --genome-size 1.8m -O g -o Raw_data/${coverage_value}x_coverage/${prefix}_${coverage_value}x_R1.fastq.gz -o Raw_data/${coverage_value}x_coverage/${prefix}_${coverage_value}x_R2.fastq.gz #calculate coverages
		done
	done
done

