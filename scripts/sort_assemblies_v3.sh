#!/usr/bin/env bash

#script for sorting assemblies into chewbbaca for every combination

outdir=$1

ls $outdir > directory_list.txt

for dir in $(cat directory_list.txt)
do
	rm -r $outdir/${dir}/Assemblies
	ls $outdir/${dir} > $outdir/${dir}/samples_list.txt
	sed -i '/Fast/d' $outdir/${dir}/samples_list.txt #remove fast things
	sed -i '/samples_list.txt/d' $outdir/${dir}/samples_list.txt
	mkdir $outdir/${dir}/Assemblies
	
	for dir2 in $(cat $outdir/${dir}/samples_list.txt)
	do
		ls $outdir/${dir}/${dir2} > $outdir/${dir}/${dir2}/results_list.txt
		sed -i '/results_list.txt/d' $outdir/${dir}/${dir2}/results_list.txt
		if grep -Fxq "Pilon" $outdir/${dir}/${dir2}/results_list.txt
		then
			cp $outdir/${dir}/${dir2}/Pilon/SPAdes/* $outdir/${dir}/Assemblies/.
			cp $outdir/${dir}/${dir2}/Pilon/SKESA/* $outdir/${dir}/Assemblies/.
		else
			cp $outdir/${dir}/${dir2}/SPAdes/${dir2}* $outdir/${dir}/Assemblies/.
			cp $outdir/${dir}/${dir2}/SKESA/* $outdir/${dir}/Assemblies/.
		fi
		
	done
	rm $outdir/${dir}/samples_list.txt	
done

rm directory_list.txt

