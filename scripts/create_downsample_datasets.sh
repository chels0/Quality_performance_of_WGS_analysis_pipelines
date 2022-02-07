#!/usr/bin/bash

coverages=( 20 50 100 )

mkdir /data/Raw_data/${coverages[0]}x_coverage
mkdir /data/Raw_data/${coverages[1]}x_coverage
mkdir /data/Raw_data/${coverages[2]}x_coverage

for sample in {1..4}
do
	ls Raw_data/PT28-${sample}* > Raw_data/path.txt
	sed 's/Raw_data\///g' Raw_data/path.txt > Raw_data/sample.txt
	sed 's/_.*//g' Raw_data/sample.txt > Raw_data/prefix_dupl.txt
	uniq Raw_data/prefix_dupl.txt Raw_data/prefix.txt

	for prefix in $(cat Raw_data/prefix.txt)
	do
		for coverage_value in ${coverages[@]}
		do
		rasusa -i Raw_data/${prefix}_R1.fastq.gz -i Raw_data/${prefix}_R2.fastq.gz --coverage ${coverage_value} --genome-size 1.8m -O g -o /data/Raw_data/${coverage_value}x_coverage/${prefix}_${coverage_value}xR1.fastq.gz -o /data/Raw_data/${coverage_value}x_coverage/${prefix}_${coverage_value}x_R2.fastq.gz
		done
	done
done



