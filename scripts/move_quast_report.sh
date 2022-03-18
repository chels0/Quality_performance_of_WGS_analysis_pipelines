#!/usr/bin/env bash

outdir=/data/Results_all

ls $outdir > directory_list.txt

for dir in $(cat directory_list.txt)
do
	ls $outdir/${dir} > $outdir/${dir}/samples_list.txt
	sed -i '/Fast/d' $outdir/${dir}/samples_list.txt #remove fast things
	sed -i '/samples_list.txt/d' $outdir/${dir}/samples_list.txt
	sed -i '/Assemblies/d' $outdir/${dir}/samples_list.txt
	sed -i '/chewBBACA/d' $outdir/${dir}/samples_list.txt
	
	for dir2 in $(cat $outdir/${dir}/samples_list.txt)
	do
		cp ${outdir}/${dir}/${dir2}/Quast/Reports/${dir2}_SPAdes* ${outdir}/${dir}/${dir2}/Quast/SPAdes/report.tsv
		cp ${outdir}/${dir}/${dir2}/Quast/Reports/${dir2}_SKESA* ${outdir}/${dir}/${dir2}/Quast/SKESA/report.tsv
	done
	rm $outdir/${dir}/samples_list.txt	
done

rm directory_list.txt
