#!/usr/bin/env nextflow


genomes = Channel.fromFilePairs(params.path_to_reads)
genomes.into { raw_data_for_Fastp; raw_data_for_trimmomatic; raw_data_for_spades; raw_data_for_bowtie2 ; raw_data_for_fastqc ; raw_data_for_skesa}

reference = Channel.value(params.path_to_reference)
reference.into { ref_for_quast_no_trim ; ref_for_quast_fastp_trim ; ref_for_quast_trimmomatic_trim }

reference_index = Channel.value(params.path_to_bowtie2_index)

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
	
	publishDir './Results/FastQC_Raw_reads', mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from raw_data_for_fastqc
	val script_folder from script_fastqc

	output:
	file "${sampleID}/*" into fastqc_raw_output
	file "${sampleID}/result.txt" into trimmomatic_input
	
	script:
	if ( params.trimmomatic == true )
	
		"""
		mkdir ${sampleID}
		# fastqc
		fastqc ${reads[0]} ${reads[1]} --outdir ${sampleID} 
		cd ${sampleID}
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
		mkdir ${sampleID}
		fastqc ${reads[0]} ${reads[1]} --outdir ${sampleID}
		echo "no trimming" > ${sampleID}/result.txt
		"""
}
	
process fastp{
	
	publishDir './Results/'+folder_name+'/Fastp/', mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from raw_data_for_Fastp
	val folders from folder.collect()
	
	output:
	file "${sampleID}/*.gz" into fastp_output
	file "${sampleID}" into fastp_reports
	val "${sampleID}" into fastp_id_name
	
	script:
	if ( params.fastp_trim_qc  == true )

		"""
		mkdir ${sampleID}
		mkdir ${sampleID}/Reports
		cd ${sampleID}
		fastp -i ../${reads[0]} -I ../${reads[1]} -o ${sampleID}_R1_trimmed_fastp.fq.gz -O ${sampleID}_R2_trimmed_fastp.fq.gz
		mv fastp* ./Reports
		"""
	
	else if ( params.fastp_trim_qc == false )
	
		"""
		mkdir ${sampleID}
		mkdir ${sampleID}/Reports
		cd ${sampleID}
		fastp -A -L -Q -G -i ../${reads[0]} -I ../${reads[1]} -o ${reads[0]} -O ${reads[1]}
		mv fastp* ./Reports
		""" 
	
}

fastp_output.into { fastp_output_fastqc ; fastp_output_SPAdes ; fastp_output_skesa }
fastp_id_name.into { fastp_id_name_fastqc ; fastp_id_name_SPAdes ; fastp_id_name_skesa }

process fastqc_post_fastp_trim {
	
	publishDir './Results/'+folder_name+'/FastQC_post_trim', mode: 'copy'
	tag "${sampleID}"
	
	input:
	file sample from fastp_output_fastqc
	val sampleID from fastp_id_name_fastqc
	
	output:
	file "${sampleID}/*" into fastqc_trim_fastp
	
	when:
	params.fastp_trim_qc == true
	
	//script:
	//if ( params.fastp_trim_qc  == true || (params.fastp_trim_qc == true && params.fastp_qc == true) )
	
	"""
	mkdir ${sampleID}
	fastqc ${sample[0]} ${sample[1]} --outdir ${sampleID}
	

	"""
}


process trimmomatic{

	publishDir './Results/Trimmed_w_Trimmomatic/Trimmomatic/', mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from raw_data_for_trimmomatic 
	file f_r from trimmomatic_input
	val settings from trimmomatic_setting
	
	output:
	file "${sampleID}/*_paired.fq.gz" into trimmomatic_output 
	val "${sampleID}" into id_name_trimmomatic
	
	when:
	params.trimmomatic == true
	
	"""
	mkdir ${sampleID}
	cd ${sampleID}
	lead="\$(head -n 1 ../${f_r})"
	trail="\$(tail -n 1 ../${f_r})"
	
	trimmomatic PE ../${reads[0]} ../${reads[1]} ${sampleID}_R1_paired.fq.gz ${sampleID}_R1_unpaired.fq.gz ${sampleID}_R2_paired.fq.gz ${sampleID}_R2_unpaired.fq.gz LEADING:\$lead TRAILING:\$trail ${settings}
	"""

}

trimmomatic_output.into { trimmomatic_output_for_SPAdes ; trimmomatic_output_for_FastQC ; trimmomatic_output_for_skesa}
id_name_trimmomatic.into { id_name_for_SPAdes ; id_name_for_FastQC ; id_name_for_skesa}


process fastqc_post_trimmomatic {
	
	publishDir './Results/Trimmed_w_Trimmomatic/FastQC_post_trim', mode: 'copy'
	tag "${sampleID}"
	
	input:
	file sample from trimmomatic_output_for_FastQC
	val sampleID from id_name_for_FastQC
	
	output:
	file "${sampleID}/*" into fastqc_raw
	
	
	"""
	mkdir ${sampleID}
	fastqc ${sample[0]} ${sample[1]} --outdir ${sampleID}

	"""
}

process spades_no_trim{
	
	publishDir './Results/No_trimming/SPAdes', mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from raw_data_for_spades
	val settings from no_trim_settings
	val script_folder from script_spades_no_trim
	
	output:
	file "${sampleID}/${sampleID}*" into spades_output
	file "${sampleID}/*" into spades_all
	val "$sampleID" into spades_id_name 
	
	when:
	params.no_trim == true

	"""
	spades.py -1 ${reads[0]} -2 ${reads[1]} -o $sampleID ${settings}
	mv ${sampleID}/scaffolds.fasta ${sampleID}/${sampleID}_spades_scaffolds.fasta
	cd ${sampleID}
	python3 ${script_folder}/remove_contaminants_spades.py ${sampleID}_spades_scaffolds 300
	"""

}


process spades_after_fastp{
	
	publishDir './Results/'+folder_name+'/SPAdes', mode: 'copy'
	tag "${sample}"
	
	input:
	file sampleID from fastp_output_SPAdes
	val sample from fastp_id_name_SPAdes
	val script_folder from script_spades_fastp
	val settings from trim_fastp_setting

	output:
	file "${sample}/${sample}*" into spades_output_fastp
	file "${sample}/*" into spades_all_fastp
	val "$sample" into spades_id_name_fastp 
 	
 	when:
 	params.fastp_trim_qc == true
 	
 	"""
	spades.py -1 ${sampleID[0]} -2 ${sampleID[1]} -o $sample ${settings}
	mv ${sample}/scaffolds.fasta ${sample}/${sample}_trimmed_fastp_spades_scaffolds.fasta
	cd ${sample}
	python3 ${script_folder}/remove_contaminants_spades.py ${sample}_trimmed_fastp_spades_scaffolds 300
	"""

}

process spades_after_trimmomatic{
	
	publishDir './Results/Trimmed_w_Trimmomatic/SPAdes', mode: 'copy'
	tag "${sample}"
	
	input:
	file sampleID from trimmomatic_output_for_SPAdes
	val sample from id_name_for_SPAdes
	val script_folder from script_spades_trimmomatic
	val settings from trim_trimmomatic_setting

	output:
	file "${sample}/${sample}*" into spades_output_trimmomatic
	file "${sample}/*" into spades_all_trimmomatic
	val "$sample" into spades_id_name_trimmomatic
 	
 	when:
 	params.trimmomatic==true

	"""
	spades.py -1 ${sampleID[0]} -2 ${sampleID[1]} -o $sample ${settings}
	mv ${sample}/scaffolds.fasta ${sample}/${sample}_trimmed_trimmomatic_spades_scaffolds.fasta
	cd ${sample}
	python3 ${script_folder}/remove_contaminants_spades.py ${sample}_trimmed_trimmomatic_spades_scaffolds 300
	"""


}

spades_output.into { spades_output_for_quast_no_trim ; spades_output_for_pilon_no_trim ; spades_no_trim_sorting }
//spades_polished_output.into { spades_polished_quast_no_trim ; spades_polished_pilon_no_trim }
spades_output_fastp.into {spades_output_for_quast_trim ; spades_output_fastp_for_pilon_trim ; spades_fastp_sorting }
spades_output_trimmomatic.into { spades_output_quast_trimmomatic ; spades_trimmomatic_sorting}

spades_id_name.into { spades_id_pilon ; spades_id_quast }

process skesa_no_trim{
	
	publishDir './Results/No_trimming/SKESA', mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from raw_data_for_skesa
	val settings from no_trim_settings
	val script_folder from script_skesa_no_trim
	
	output:
	file "${sampleID}/${sampleID}*" into skesa_output_no_trim
	val "$sampleID" into skesa_id_name_no_trim 
	
	when:
	params.no_trim == true

	"""
	mkdir ${sampleID}
	skesa --reads ${reads[0]},${reads[1]} --use_paired_ends > ${sampleID}/${sampleID}_skesa_contigs.fasta
	cd ${sampleID}
	python3 ${script_folder}/remove_contaminants.py ${sampleID}_skesa_contigs 300
	"""

}

process skesa_after_fastp{
	
	publishDir './Results/'+folder_name+'/SKESA', mode: 'copy'
	tag "${sample}"
	
	input:
	file sampleID from fastp_output_skesa
	val sample from fastp_id_name_skesa
	val script_folder from script_skesa_fastp
	
	output:
	file "${sample}/${sample}*" into skesa_output_fastp
	val "$sample" into skesa_id_name_fastp 
	
	when:
	params.fastp_trim_qc == true

	"""
	echo $sample
	mkdir ${sample}
	skesa --reads ${sampleID[0]},${sampleID[1]} --use_paired_ends > ${sample}/${sample}_trimmed_fastp_skesa_contigs.fasta
	cd ${sample}
	python3 ${script_folder}/remove_contaminants.py ${sample}_trimmed_fastp_skesa_contigs 300
	"""
}

process skesa_after_trimmomatic{
	
	publishDir './Results/Trimmed_w_Trimmomatic/SKESA', mode: 'copy'
	
	tag "${sample}"
	
	input:
	file sampleID from trimmomatic_output_for_skesa
	val sample from id_name_for_skesa
	val script_folder from script_skesa_trimmomatic
	
	output:
	file "${sample}/${sample}*" into skesa_output_trimmomatic
	val "${sample}" into skesa_id_name_trimmomatic
	
	when:
	params.trimmomatic == true

	"""
	mkdir ${sample}
	skesa --reads ${sampleID[0]},${sampleID[1]} --use_paired_ends > ${sample}/${sample}_trimmed_trimmomatic_skesa_contigs.fasta
	cd ${sample}
	python3 ${script_folder}/remove_contaminants.py ${sample}_trimmed_trimmomatic_skesa_contigs 300
	"""
}

skesa_output_no_trim.into { skesa_output_quast_no_trim ; skesa_no_trim_sorting }
skesa_output_fastp.into { skesa_output_quast_fastp ; skesa_fastp_sorting }
skesa_output_trimmomatic.into { skesa_output_quast_trimmomatic ; skesa_trimmomatic_sorting }


//merge spades id with spades output files 
spades_merge = spades_id_pilon.merge(spades_output_for_pilon_no_trim)
//join the raw data for bowtie with the spades results by using sampleID as key
spades_raw_data_join = spades_merge.join(raw_data_for_bowtie2)

process pilon_improvement_no_trim {

	publishDir './Results/No_trimming/Pilon', mode: 'copy'
	
	input:
	tuple sampleID, file(scaffold), file(scaffold_polish), file(reads) from spades_raw_data_join
	
	output:
	file "${sampleID}/${scaffold}_improved.fasta*" into pilon_output
	
	script:
	index_base = scaffold.toString() - ~/.fasta/
	index_base_polish = scaffold_polish.toString() - ~/.fasta/
	"""
	mkdir ${sampleID}
	mkdir ${sampleID}/index
	mkdir ${sampleID}/index_polish
	bowtie2-build ${scaffold} ${sampleID}/index/$index_base
	bowtie2 -x ${sampleID}/index/${index_base} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}/${sampleID}_sorted_alignment.bam
	samtools index ${sampleID}/${sampleID}_sorted_alignment.bam
	
	pilon --genome ${scaffold} --frags ${sampleID}/${sampleID}_sorted_alignment.bam --output ${sampleID}/${scaffold}_improved
	rm *.bam
	rm *.bai
	"""

}



process quast_no_trimming{
	
	publishDir './Results/No_trimming/Quast', mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	file scaffold from spades_output_for_quast_no_trim
	file scaffold_skesa from skesa_output_quast_no_trim
	val reference from ref_for_quast_no_trim
	val sampleID from spades_id_quast
	val sampleID_skesa from skesa_id_name_no_trim 
	
	output:
	file "SPAdes/${sampleID}/*" into quast_output_spades_no_trim
	file "SKESA/${sampleID_skesa}/*" into quast_output_skesa_no_trim
	val "$sampleID" into quast_id_name
	
	"""
	quast.py ${scaffold[0]} -o SPAdes/$sampleID/Not_polished -r $reference
	quast.py ${scaffold[1]} -o SPAdes/$sampleID/Polished -r $reference
	quast.py ${scaffold_skesa[0]} -o SKESA/$sampleID_skesa/Not_polished -r $reference
	quast.py ${scaffold_skesa[1]} -o SKESA/$sampleID_skesa/Polished -r $reference
	"""
}

process quast_after_fastp{
	
	publishDir './Results/'+folder_name+'/Quast', mode: 'copy'
	tag "${sampleID}"
	
	input:
	file scaffold from spades_output_for_quast_trim
	file scaffold_skesa from skesa_output_quast_fastp
	val reference from ref_for_quast_fastp_trim
	val sampleID from spades_id_name_fastp 
	val sampleID_skesa from skesa_id_name_fastp 
	
	output:
	file "SPAdes/${sampleID}/*" into quast_output_spades_fastp
	file "SKESA/${sampleID_skesa}/*" into quast_output_skesa_fastp
	val "$sampleID" into quast_id_name_fastp
	
	"""
	quast.py ${scaffold[0]} -o SPAdes/$sampleID/Not_polished -r $reference
	quast.py ${scaffold[1]} -o SPAdes/$sampleID/Polished -r $reference
	quast.py ${scaffold_skesa[0]} -o SKESA/$sampleID_skesa/Not_polished -r $reference
	quast.py ${scaffold_skesa[1]} -o SKESA/$sampleID_skesa/Polished -r $reference
	"""
}

process quast_after_trimmomatic{
	
	publishDir './Results/Trimmed_w_Trimmomatic/Quast', mode: 'copy'
	tag "${sampleID}"
	
	input:
	file scaffold from spades_output_quast_trimmomatic
	file scaffold_skesa from skesa_output_quast_trimmomatic
	val reference from ref_for_quast_trimmomatic_trim
	val sampleID from spades_id_name_trimmomatic
	val sampleID_skesa from skesa_id_name_trimmomatic
	
	output:
	file "SPAdes/${sampleID}/*" into quast_output_trimmomatic
	file "SKESA/${sampleID_skesa}/*" into quast_output_trimmomatic_skesa
	val "$sampleID" into quast_id_name_trimmomatic
	
	"""
	quast.py ${scaffold[0]} -o SPAdes/$sampleID/Not_polished -r $reference
	quast.py ${scaffold[1]} -o SPAdes/$sampleID/Polished -r $reference
	quast.py ${scaffold_skesa[0]} -o SKESA/$sampleID_skesa/Not_polished -r $reference
	quast.py ${scaffold_skesa[1]} -o SKESA/$sampleID_skesa/Polished -r $reference
	"""
}


//result.view { it }
