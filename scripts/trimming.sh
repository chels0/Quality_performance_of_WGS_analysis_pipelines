# unzip fastqc results
		unzip '*.zip' 
		cd ${sampleID}_R1_fastqc
		# split the data file after each >>
		csplit fastqc_data.txt '/>>/' '{*}' 
		# save first and second column to csv file
		awk '{ print \$1,\$2 }' xx03 > xx03_relevant
		# remove first row
		tail -n +2 xx03_relevant > xx03.csv 
		
		# determine leading for trimmomatic
		echo "\$(python3.6 $script_folder/fastqc_trim_beginning.py)" > lead_R1.txt
		# determine trailing for trimmomatic
		echo "\$(python3.6 $script_folder/fastqc_trim.py)" > trail_R1.txt 
		
		mv lead_R1.txt ..
		mv trail_R1.txt ..
		
		cd ../${sampleID}_R2_fastqc
		csplit fastqc_data.txt '/>>/' '{*}'
		awk '{ print \$1,\$2 }' xx03 > xx03_relevant
		tail -n +2 xx03_relevant > xx03.csv
		
		echo "\$(python3.6 $script_folder/fastqc_trim_beginning.py)" > lead_R2.txt
		echo "\$(python3.6 $script_folder/fastqc_trim.py)" > trail_R2.txt
		
		mv lead_R2.txt ..
		mv trail_R2.txt ..
		
		cd ../
		
		# define parameters for python
		l_R1="\$(cat lead_R1.txt)" 
		l_R2="\$(cat lead_R2.txt)"
		t_R1="\$(cat trail_R1.txt)"
		t_R2="\$(cat trail_R2.txt)"
		# save highest leading value in txt file
		echo "\$(python3.6 $script_folder/check_trim.py \$l_R1 \$l_R2)" > result.txt 
		# save highest trailing value in txt file
		echo "\$(python3.6 $script_folder/check_trim_trail.py \$t_R1 \$t_R2)" >> result.txt 
