#!/usr/bin/env bash

#script for sorting assemblies into chewbbaca for every combination

ls Results/ > directory_list.txt

for dir in $(cat directory_list.txt)
do
	rm -r Results/${dir}/Assemblies
	ls Results/${dir} > Results/${dir}/samples_list.txt
	sed -i '/Fast/d' Results/${dir}/samples_list.txt #remove fast things
	sed -i '/samples_list.txt/d' Results/${dir}/samples_list.txt
	mkdir Results/${dir}/Assemblies
	
	for dir2 in $(cat Results/${dir}/samples_list.txt)
	do
		ls Results/${dir}/${dir2} > Results/${dir}/${dir2}/results_list.txt
		sed -i '/results_list.txt/d' Results/${dir}/${dir2}/results_list.txt
		if grep -Fxq "Pilon" Results/${dir}/${dir2}/results_list.txt
		then
			cp Results/${dir}/${dir2}/Pilon/SPAdes/* Results/${dir}/Assemblies/.
			cp Results/${dir}/${dir2}/Pilon/SKESA/* Results/${dir}/Assemblies/.
		else
			cp Results/${dir}/${dir2}/SPAdes/${dir2}* Results/${dir}/Assemblies/.
			cp Results/${dir}/${dir2}/SKESA/* Results/${dir}/Assemblies/.
		fi
		
	done
	rm Results/${dir}/samples_list.txt	
done

rm directory_list.txt
