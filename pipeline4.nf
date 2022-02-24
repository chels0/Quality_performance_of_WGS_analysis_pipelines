#!/usr/bin/env nextflow


genomes = Channel.fromFilePairs(params.path_to_reads)
genomes.into { raw_data_for_Fastp; raw_data_for_trimmomatic; raw_data_for_spades; raw_data_for_bowtie2_spades_no_trim ; raw_data_for_bowtie2_skesa_no_trim ; raw_data_for_bowtie2_spades_fastp ; raw_data_for_bowtie2_skesa_fastp ; raw_data_for_bowtie2_spades_trimmomatic ; raw_data_for_bowtie2_skesa_trimmomatic ; raw_data_for_fastqc ; raw_data_for_skesa}

reference = Channel.value(params.path_to_reference)
reference.into { ref_for_quast_no_trim ; ref_for_quast_fastp_trim ; ref_for_quast_trimmomatic_trim }

//reference_index = Channel.value(params.path_to_bowtie2_index)

if ( params.fastp_trim_qc == false )
	folder_name = 'No_trimming'
	
else if ( params.fastp_trim_qc  == true )
	folder_name = 'Trimmed_w_fastp'


folder = Channel.value(folder_name)

scripts_folder = Channel.value(params.path_to_scripts)
scripts_folder.into { script_fastqc ; script_spades_no_trim ; script_skesa_no_trim ; script_spades_fastp ; script_skesa_fastp ; script_spades_trimmomatic ; script_skesa_trimmomatic}

trimmomatic_setting = Channel.value(params.trimmomatic_set)
spades_setting = Channel.value(params.spades_set)
spades_setting.into { no_trim_settings ; trim_trimmomatic_setting ; trim_fastp_setting }

process fastqc_raw {
	
	publishDir "./Results/${sampleID}/FastQC_Raw_reads", mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from raw_data_for_fastqc
	val script_folder from script_fastqc

	output:
	file "*" into fastqc_raw_output
	file "result.txt" into trimmomatic_input
	
	script:
	if ( params.trimmomatic == true )
	
		"""
		# fastqc
		fastqc ${reads[0]} ${reads[1]} --outdir . 
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
		
		"""
	else
		"""
		fastqc ${reads[0]} ${reads[1]} --outdir .
		echo "no trimming" > result.txt
		"""
}
	
process fastp{
	
	publishDir "./Results/${sampleID}/"+folder_name+"/Fastp/", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from raw_data_for_Fastp
	val folders from folder.collect()
	
	output:
	tuple sampleID, file("${sampleID}*") into fastp_output
	file "Reports/*" into fastp_reports

	
	script:
	if ( params.fastp_trim_qc  == true )

		"""
		mkdir Reports
		fastp -i ${reads[0]} -I ${reads[1]} -o ${sampleID}_R1_trimmed_fastp.fq.gz -O ${sampleID}_R2_trimmed_fastp.fq.gz
		mv fastp* Reports/
		"""
	
	else if ( params.fastp_trim_qc == false )
	
		"""
		mkdir Reports
		fastp -A -L -Q -G -i ${reads[0]} -I ${reads[1]} -o ${sampleID}_R1.fq.gz -O ${sampleID}_R2.fq.gz
		mv fastp* ./Reports
		""" 
	
}

fastp_output.into { fastp_output_fastqc ; fastp_output_SPAdes ; fastp_output_skesa }

process fastqc_post_fastp_trim {
	
	publishDir "./Results/${sampleID}/"+folder_name+"/FastQC_post_trim", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from fastp_output_fastqc
	
	output:
	file "*" into fastqc_trim_fastp
	
	when:
	params.fastp_trim_qc == true
	
	"""
	fastqc ${reads[0]} ${reads[1]} --outdir .

	"""
}


process trimmomatic{

	publishDir "./Results/$sampleID/Trimmed_w_Trimmomatic/Trimmomatic/", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from raw_data_for_trimmomatic 
	file f_r from trimmomatic_input
	val settings from trimmomatic_setting
	
	output:
	tuple sampleID, file("*_paired.fq.gz") into trimmomatic_output 
	
	when:
	params.trimmomatic == true
	
	"""
	lead="\$(head -n 1 ${f_r})"
	trail="\$(tail -n 1 ${f_r})"
	
	trimmomatic PE ${reads[0]} ${reads[1]} ${sampleID}_R1_paired.fq.gz ${sampleID}_R1_unpaired.fq.gz ${sampleID}_R2_paired.fq.gz ${sampleID}_R2_unpaired.fq.gz LEADING:\$lead TRAILING:\$trail ${settings}
	"""

}

trimmomatic_output.into { trimmomatic_output_for_SPAdes ; trimmomatic_output_for_FastQC ; trimmomatic_output_for_skesa}

process fastqc_post_trimmomatic {
	
	publishDir "./Results/${sampleID}/Trimmed_w_Trimmomatic/FastQC_post_trim", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from trimmomatic_output_for_FastQC

	output:
	file "*" into fastqc_raw
	
	
	"""
	fastqc ${reads[0]} ${reads[1]} --outdir .

	"""
}

process spades_no_trim{
	
	publishDir "./Results/${sampleID}/No_trimming/SPAdes", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from raw_data_for_spades
	val settings from no_trim_settings
	val script_folder from script_spades_no_trim
	
	output:
	tuple sampleID, file("${sampleID}*") into spades_output
	file "*" into spades_all
	
	when:
	params.no_trim == true

	"""
	spades.py -1 ${reads[0]} -2 ${reads[1]} -o . ${settings}
	mv scaffolds.fasta ./${sampleID}_spades_scaffolds.fasta
	python3 ${script_folder}/remove_contaminants_spades.py ${sampleID}_spades_scaffolds 300
	"""

}


process spades_after_fastp{
	
	publishDir "./Results/${sampleID}/"+folder_name+"/SPAdes", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from fastp_output_SPAdes
	val script_folder from script_spades_fastp
	val settings from trim_fastp_setting

	output:
	tuple sampleID, file("${sampleID}*") into spades_output_fastp
	file "*" into spades_all_fastp

 	
 	when:
 	params.fastp_trim_qc == true
 	
 	"""
	spades.py -1 ${reads[0]} -2 ${reads[1]} -o . ${settings}
	mv scaffolds.fasta ${sampleID}_trimmed_fastp_spades_scaffolds.fasta
	python3 ${script_folder}/remove_contaminants_spades.py ${sampleID}_trimmed_fastp_spades_scaffolds 300
	"""

}

process spades_after_trimmomatic{
	
	publishDir "./Results/${sampleID}/Trimmed_w_Trimmomatic/SPAdes", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from trimmomatic_output_for_SPAdes
	val script_folder from script_spades_trimmomatic
	val settings from trim_trimmomatic_setting

	output:
	tuple sampleID, file("${sampleID}*") into spades_output_trimmomatic
	file "*" into spades_all_trimmomatic

 	when:
 	params.trimmomatic==true

	"""
	spades.py -1 ${reads[0]} -2 ${reads[1]} -o . ${settings}
	mv scaffolds.fasta ${sampleID}_trimmed_trimmomatic_spades_scaffolds.fasta
	python3 ${script_folder}/remove_contaminants_spades.py ${sampleID}_trimmed_trimmomatic_spades_scaffolds 300
	"""


}

spades_output.into { spades_output_for_quast_no_trim ; spades_output_for_pilon_no_trim ; spades_no_trim_sorting }
spades_output_fastp.into {spades_output_for_quast_fastp ; spades_output_fastp_for_pilon_trim ; spades_fastp_sorting }
spades_output_trimmomatic.into { spades_output_quast_trimmomatic ; spades_output_trimmomatic_for_pilon_trim ; spades_trimmomatic_sorting}


process skesa_no_trim{
	
	publishDir "./Results/${sampleID}/No_trimming/SKESA", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from raw_data_for_skesa
	val settings from no_trim_settings
	val script_folder from script_skesa_no_trim
	
	output:
	tuple sampleID, file("${sampleID}*") into skesa_output_no_trim
	
	when:
	params.no_trim == true

	"""
	skesa --reads ${reads[0]},${reads[1]} --use_paired_ends > ${sampleID}_skesa_contigs.fasta
	python3 ${script_folder}/remove_contaminants.py ${sampleID}_skesa_contigs 300
	"""

}

process skesa_after_fastp{
	
	publishDir "./Results/${sampleID}"+folder_name+"/SKESA", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from fastp_output_skesa
	val script_folder from script_skesa_fastp
	
	output:
	tuple sampleID, file("${sampleID}*") into skesa_output_fastp
	
	when:
	params.fastp_trim_qc == true

	"""
	skesa --reads ${reads[0]},${reads[1]} --use_paired_ends > ${sampleID}_trimmed_fastp_skesa_contigs.fasta
	python3 ${script_folder}/remove_contaminants.py ${sampleID}_trimmed_fastp_skesa_contigs 300
	"""
}

process skesa_after_trimmomatic{
	
	publishDir "./Results/${sampleID}/Trimmed_w_Trimmomatic/SKESA", mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from trimmomatic_output_for_skesa
	val script_folder from script_skesa_trimmomatic
	
	output:
	tuple sampleID, file("${sampleID}*") into skesa_output_trimmomatic
	
	when:
	params.trimmomatic == true

	"""
	skesa --reads ${reads[0]},${reads[1]} --use_paired_ends > ${sampleID}_trimmed_trimmomatic_skesa_contigs.fasta
	python3 ${script_folder}/remove_contaminants.py ${sampleID}_trimmed_trimmomatic_skesa_contigs 300
	"""
}

skesa_output_no_trim.into { skesa_output_quast_no_trim ; skesa_no_trim_sorting ; skesa_no_trim_pilon }
skesa_output_fastp.into { skesa_output_quast_fastp ; skesa_fastp_sorting ; skesa_fastp_pilon }
skesa_output_trimmomatic.into { skesa_output_quast_trimmomatic ; skesa_trimmomatic_sorting ; skesa_trimmomatic_pilon }

//join the raw data for bowtie with the spades results by using sampleID as key
spades_raw_data_join = spades_output_for_pilon_no_trim.join(raw_data_for_bowtie2_spades_no_trim)
spades_raw_trimmomatic_join = spades_output_fastp_for_pilon_trim.join(raw_data_for_bowtie2_spades_trimmomatic)
spades_raw_fastp_join = spades_output_trimmomatic_for_pilon_trim.join(raw_data_for_bowtie2_spades_fastp)

//join the raw data for bowtie with skesa results by using sampleID as key
skesa_raw_data_join = skesa_no_trim_pilon.join(raw_data_for_bowtie2_skesa_no_trim)
skesa_raw_trimmomatic_join = skesa_trimmomatic_pilon.join(raw_data_for_bowtie2_skesa_trimmomatic)
skesa_raw_fastp_join = skesa_fastp_pilon.join(raw_data_for_bowtie2_skesa_fastp)

process pilon_post_spades_no_trim {

	publishDir "./Results/${sampleID}/No_trimming/Pilon/SPAdes", mode: 'copy'
	
	input:
	tuple sampleID, file(scaffold), file(reads) from spades_raw_data_join
	
	output:
	tuple sampleID, file("${sampleID}*") into pilon_output_spades_no_trim
	
	when:
	params.assembly_improvement==true && params.no_trim==true
	
	script:
	index_base = scaffold[0].toString() - ~/.fasta/
	index_base_polish = scaffold[1].toString() - ~/.fasta/
	"""
	mkdir index
	bowtie2-build ${scaffold[0]} index/$index_base
	bowtie2 -x index/${index_base} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_sorted_alignment.bam
	samtools index ${sampleID}_sorted_alignment.bam
	
	pilon --genome ${scaffold[0]} --frags ${sampleID}_sorted_alignment.bam --output ${index_base}_improved
	
	rm *.bam
	rm *.bai
	
	mkdir index_polish
	bowtie2-build ${scaffold[1]} index_polish/$index_base_polish
	bowtie2 -x index_polish/${index_base_polish} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_polished_sorted_alignment.bam
	samtools index ${sampleID}_polished_sorted_alignment.bam
	
	pilon --genome ${scaffold[1]} --frags ${sampleID}_polished_sorted_alignment.bam --output ${index_base_polish}_improved
	
	rm *.bam
	rm *.bai
	"""
	

}

process pilon_post_spades_fastp {

	publishDir "./Results/${sampleID}/"+folder_name+"/Pilon/SPAdes", mode: 'copy'
	
	input:
	tuple sampleID, file(scaffold), file(reads) from spades_raw_fastp_join
	
	output:
	tuple sampleID, file("${sampleID}*") into pilon_output_spades_fastp
	
	when:
	params.assembly_improvement==true && params.fastp_trim_qc==true
	
	script:
	index_base = scaffold[0].toString() - ~/.fasta/
	index_base_polish = scaffold[1].toString() - ~/.fasta/
	
	"""
	mkdir index
	bowtie2-build ${scaffold[0]} index/$index_base
	bowtie2 -x index/${index_base} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_sorted_alignment.bam
	samtools index ${sampleID}_sorted_alignment.bam
	
	pilon --genome ${scaffold[0]} --frags ${sampleID}_sorted_alignment.bam --output ${index_base}_improved
	
	rm *.bam
	rm *.bai
	
	mkdir index_polish
	bowtie2-build ${scaffold[1]} index_polish/$index_base_polish
	bowtie2 -x index_polish/${index_base_polish} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_polished_sorted_alignment.bam
	samtools index ${sampleID}_polished_sorted_alignment.bam
	
	pilon --genome ${scaffold[1]} --frags ${sampleID}_polished_sorted_alignment.bam --output ${index_base_polish}_improved
	
	rm *.bam
	rm *.bai
	"""
	

}

process pilon_post_spades_trimmomatic {

	publishDir "./Results/${sampleID}/Trimmed_w_Trimmomatic/Pilon/SPAdes", mode: 'copy'
	
	input:
	tuple sampleID, file(scaffold), file(reads) from spades_raw_trimmomatic_join
	
	output:
	tuple sampleID, file("${sampleID}*") into pilon_output_spades_trimmomatic
	
	when:
	params.assembly_improvement==true && params.trimmomatic==true
	
	script:
	index_base = scaffold[0].toString() - ~/.fasta/
	index_base_polish = scaffold[1].toString() - ~/.fasta/
	"""
	mkdir index
	bowtie2-build ${scaffold[0]} index/$index_base
	bowtie2 -x index/${index_base} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_sorted_alignment.bam
	samtools index ${sampleID}_sorted_alignment.bam
	
	pilon --genome ${scaffold[0]} --frags ${sampleID}_sorted_alignment.bam --output ${index_base}_improved
	
	rm *.bam
	rm *.bai
	
	mkdir index_polish
	bowtie2-build ${scaffold[1]} index_polish/$index_base_polish
	bowtie2 -x index_polish/${index_base_polish} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_polished_sorted_alignment.bam
	samtools index ${sampleID}_polished_sorted_alignment.bam
	
	pilon --genome ${scaffold[1]} --frags ${sampleID}_polished_sorted_alignment.bam --output ${index_base_polish}_improved
	
	rm *.bam
	rm *.bai
	"""
	

}

process pilon_post_skesa_no_trim {

	publishDir "./Results/${sampleID}/No_trimming/Pilon/SKESA", mode: 'copy'
	
	input:
	tuple sampleID, file(scaffold), file(reads) from skesa_raw_data_join
	
	output:
	tuple sampleID, file("${sampleID}*") into pilon_output_skesa_no_trim
	
	when:
	params.assembly_improvement==true && params.no_trim==true
	
	script:
	index_base = scaffold[0].toString() - ~/.fasta/
	index_base_polish = scaffold[1].toString() - ~/.fasta/
	"""
	mkdir index
	bowtie2-build ${scaffold[0]} index/$index_base
	bowtie2 -x index/${index_base} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_sorted_alignment.bam
	samtools index ${sampleID}_sorted_alignment.bam
	
	pilon --genome ${scaffold[0]} --frags ${sampleID}_sorted_alignment.bam --output ${index_base}_improved
	
	rm *.bam
	rm *.bai
	
	mkdir index_polish
	bowtie2-build ${scaffold[1]} index_polish/$index_base_polish
	bowtie2 -x index_polish/${index_base_polish} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_polished_sorted_alignment.bam
	samtools index ${sampleID}_polished_sorted_alignment.bam
	
	pilon --genome ${scaffold[1]} --frags ${sampleID}_polished_sorted_alignment.bam --output ${index_base_polish}_improved
	
	rm *.bam
	rm *.bai
	"""
	

}

process pilon_post_skesa_fastp {

	publishDir "./Results/${sampleID}/"+folder_name+"/Pilon/SKESA", mode: 'copy'
	
	input:
	tuple sampleID, file(scaffold), file(reads) from skesa_raw_fastp_join
	
	output:
	tuple sampleID, file("${sampleID}*") into pilon_output_skesa_fastp
	
	params.assembly_improvement==true && params.fastp_trim_qc==true
	
	script:
	index_base = scaffold[0].toString() - ~/.fasta/
	index_base_polish = scaffold[1].toString() - ~/.fasta/
	"""
	mkdir index
	bowtie2-build ${scaffold[0]} index/$index_base
	bowtie2 -x index/${index_base} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_sorted_alignment.bam
	samtools index ${sampleID}_sorted_alignment.bam
	
	pilon --genome ${scaffold[0]} --frags ${sampleID}_sorted_alignment.bam --output ${index_base}_improved
	
	rm *.bam
	rm *.bai
	
	mkdir index_polish
	bowtie2-build ${scaffold[1]} index_polish/$index_base_polish
	bowtie2 -x index_polish/${index_base_polish} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_polished_sorted_alignment.bam
	samtools index ${sampleID}_polished_sorted_alignment.bam
	
	pilon --genome ${scaffold[1]} --frags ${sampleID}_polished_sorted_alignment.bam --output ${index_base_polish}_improved
	
	rm *.bam
	rm *.bai
	"""
	

}

process pilon_post_skesa_trimmomatic {

	publishDir "./Results/${sampleID}/Trimmed_w_Trimmomatic/Pilon/SKESA", mode: 'copy'
	
	input:
	tuple sampleID, file(scaffold), file(reads) from skesa_raw_trimmomatic_join
	
	output:
	tuple sampleID, file("${sampleID}*") into pilon_output_skesa_trimmomatic
	
	when:
	params.assembly_improvement==true && params.trimmomatic==true
	
	script:
	index_base = scaffold[0].toString() - ~/.fasta/
	index_base_polish = scaffold[1].toString() - ~/.fasta/
	"""
	mkdir index
	bowtie2-build ${scaffold[0]} index/$index_base
	bowtie2 -x index/${index_base} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_sorted_alignment.bam
	samtools index ${sampleID}_sorted_alignment.bam
	
	pilon --genome ${scaffold[0]} --frags ${sampleID}_sorted_alignment.bam --output ${index_base}_improved
	
	rm *.bam
	rm *.bai
	
	mkdir index_polish
	bowtie2-build ${scaffold[1]} index_polish/$index_base_polish
	bowtie2 -x index_polish/${index_base_polish} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_polished_sorted_alignment.bam
	samtools index ${sampleID}_polished_sorted_alignment.bam
	
	pilon --genome ${scaffold[1]} --frags ${sampleID}_polished_sorted_alignment.bam --output ${index_base_polish}_improved
	
	rm *.bam
	rm *.bai
	"""
	

}

process quast_no_trimming{
	
	publishDir "./Results/${sampleID}/No_trimming/Quast/No_improvement", mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from spades_output_for_quast_no_trim
	tuple sampleID_skesa, file(assembly_skesa) from skesa_output_quast_no_trim
	val reference from ref_for_quast_no_trim
	
	output:
	file "SPAdes/*" into quast_output_spades_no_trim
	file "SKESA/*" into quast_output_skesa_no_trim

	
	"""
	quast.py ${assembly[0]} -o SPAdes/Not_polished -r $reference
	quast.py ${assembly[1]} -o SPAdes/Polished -r $reference
	quast.py ${assembly_skesa[0]} -o SKESA/Not_polished -r $reference
	quast.py ${assembly_skesa[1]} -o SKESA/Polished -r $reference
	"""
}

process quast_after_fastp{
	
	publishDir "./Results/${sampleID}/"+folder_name+"/Quast/No_improvement", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from spades_output_for_quast_fastp
	tuple sampleID_skesa, file(assembly_skesa) from skesa_output_quast_fastp
	val reference from ref_for_quast_no_trim
	
	output:
	file "SPAdes/*" into quast_output_spades_fastp
	file "SKESA/*" into quast_output_skesa_fastp

	
	"""
	quast.py ${assembly[0]} -o SPAdes/Not_polished -r $reference
	quast.py ${assembly[1]} -o SPAdes/Polished -r $reference
	quast.py ${assembly_skesa[0]} -o SKESA/Not_polished -r $reference
	quast.py ${assembly_skesa[1]} -o SKESA/Polished -r $reference
	"""
}

process quast_after_trimmomatic{
	
	publishDir "./Results/${sampleID}/Trimmed_w_Trimmomatic/Quast/No_improvement", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from spades_output_quast_trimmomatic
	tuple sampleID_skesa, file(assembly_skesa) from skesa_output_quast_trimmomatic
	val reference from ref_for_quast_no_trim
	
	output:
	file "SPAdes/*" into quast_output_spades_trimmomatic
	file "SKESA/*" into quast_output_skesa_trimmomatic

	
	"""
	quast.py ${assembly[0]} -o SPAdes/Not_polished -r $reference
	quast.py ${assembly[1]} -o SPAdes/Polished -r $reference
	quast.py ${assembly_skesa[0]} -o SKESA/Not_polished -r $reference
	quast.py ${assembly_skesa[1]} -o SKESA/Polished -r $reference
	"""
}

process quast_no_trimming_improved{
	
	publishDir "./Results/${sampleID}/No_trimming/Quast/Improved_assembly", mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from pilon_output_spades_no_trim
	tuple sampleID_skesa, file(assembly_skesa) from pilon_output_skesa_no_trim
	val reference from ref_for_quast_no_trim
	
	output:
	file "SPAdes/*" into quast_output_spades_improved_no_trim
	file "SKESA/*" into quast_output_skesa_improved_no_trim

	
	"""
	quast.py ${assembly[0]} -o SPAdes/Not_polished -r $reference
	quast.py ${assembly[1]} -o SPAdes/Polished -r $reference
	quast.py ${assembly_skesa[0]} -o SKESA/Not_polished -r $reference
	quast.py ${assembly_skesa[1]} -o SKESA/Polished -r $reference
	"""
}

process quast_after_fastp_improved{
	
	publishDir "./Results/${sampleID}/No_trimming/Quast/Improved_assembly", mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from pilon_output_spades_fastp
	tuple sampleID_skesa, file(assembly_skesa) from pilon_output_skesa_fastp
	val reference from ref_for_quast_no_trim
	
	output:
	file "SPAdes/*" into quast_output_spades_improved_fastp
	file "SKESA/*" into quast_output_skesa_improved_fastp

	
	"""
	quast.py ${assembly[0]} -o SPAdes/Not_polished -r $reference
	quast.py ${assembly[1]} -o SPAdes/Polished -r $reference
	quast.py ${assembly_skesa[0]} -o SKESA/Not_polished -r $reference
	quast.py ${assembly_skesa[1]} -o SKESA/Polished -r $reference
	"""
}

process quast_after_trimmomatic_improved{
	
	publishDir "./Results/${sampleID}/No_trimming/Quast/Improved_assembly", mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from pilon_output_spades_trimmomatic
	tuple sampleID_skesa, file(assembly_skesa) from pilon_output_skesa_trimmomatic
	val reference from ref_for_quast_no_trim
	
	output:
	file "SPAdes/*" into quast_output_spades_improved_trimmomatic
	file "SKESA/*" into quast_output_skesa_improved_trimmomatic

	
	"""
	quast.py ${assembly[0]} -o SPAdes/Not_polished -r $reference
	quast.py ${assembly[1]} -o SPAdes/Polished -r $reference
	quast.py ${assembly_skesa[0]} -o SKESA/Not_polished -r $reference
	quast.py ${assembly_skesa[1]} -o SKESA/Polished -r $reference
	"""
}



//result.view { it }
