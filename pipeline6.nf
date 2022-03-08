#!/usr/bin/env nextflow

genomes = Channel.fromFilePairs(params.path_to_reads).ifEmpty { error "Missing input reads, directory is either empty or missing read pair" }

genomes.into { raw_data_for_Fastp_trim; raw_data_for_Fastp_raw ; raw_data_for_trimmomatic; raw_data_for_spades; raw_data_for_bowtie2_spades_no_trim ; raw_data_for_bowtie2_skesa_no_trim ; raw_data_for_bowtie2_spades_fastp ; raw_data_for_bowtie2_skesa_fastp ; raw_data_for_bowtie2_spades_trimmomatic ; raw_data_for_bowtie2_skesa_trimmomatic ; raw_data_for_fastqc ; raw_data_for_skesa}

reference = Channel.value(params.path_to_reference).ifEmpty { error "No reference genome. Add reference genome to Raw_data directory" }
reference.into { ref_for_quast_no_trim ; ref_for_quast_fastp_trim ; ref_for_quast_trimmomatic_trim }

scripts_folder = Channel.value(params.path_to_scripts)
scripts_folder.into { script_fastqc ; script_spades_no_trim ; script_skesa_no_trim ; script_spades_fastp ; script_skesa_fastp ; script_spades_trimmomatic ; script_skesa_trimmomatic}

trimmomatic_setting = Channel.value(params.trimmomatic_set)
spades_setting = Channel.value(params.spades_set)
spades_setting.into { no_trim_settings ; trim_trimmomatic_setting ; trim_fastp_setting }
filter_setting = Channel.value(params.filter_set)
filter_setting.into { no_trim_skesa_filter_setting ; fastp_skesa_filter_setting ; trimmomatic_skesa_filter_setting ; no_trim_spades_filter_setting ; fastp_spades_filter_setting ; trimmomatic_spades_filter_setting }


out_dir = ""
filter = params.filter_set
spades = params.spades_set

if ( params.no_trim == true )
	out_dir = out_dir+"No_trimming_"+spades+"_"
if ( params.trimmomatic == true )
	out_dir = out_dir+"After_trimmomatic_"+spades+"_"
if ( params.fastp_trim_qc == true )
	out_dir = out_dir+"After_fastp_"+spades+"_"
if ( params.filter_contigs == true )
	out_dir = out_dir+filter+"filter_"
if ( params.assembly_improvement == true )
	out_dir = out_dir+"Improved_pilon"

println { out_dir }

process fastqc_raw {
	
	publishDir "./Results/${out_dir}/${sampleID}/FastQC_Raw_reads", mode: 'copy'
	
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

process fastp_qc_raw{
	
	publishDir "./Results/${out_dir}/${sampleID}/Fastp_QC_Raw_reads", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from raw_data_for_Fastp_raw
	
	output:
	file "Reports/*" into fastp_reports_raw_reads
	
	"""
	mkdir Reports
	fastp -A -L -Q -G -i ${reads[0]} -I ${reads[1]} -o ${sampleID}_R1.fq.gz -O ${sampleID}_R2.fq.gz
	mv fastp* ./Reports
	""" 
	
}
	
process fastp_trimming{
	
	publishDir "./Results/${out_dir}/${sampleID}/Fastp/", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from raw_data_for_Fastp_trim
	
	output:
	tuple sampleID, file("${sampleID}*") into fastp_output
	file "Reports/*" into fastp_reports

	when:
	params.fastp_trim_qc == true
	
	"""
	mkdir Reports
	fastp -i ${reads[0]} -I ${reads[1]} -o ${sampleID}_R1_trimmed_fastp.fq.gz -O ${sampleID}_R2_trimmed_fastp.fq.gz
	mv fastp* Reports/
	"""	
}



fastp_output.into { fastp_output_fastqc ; fastp_output_SPAdes ; fastp_output_skesa }

process fastqc_post_fastp_trim {
	
	publishDir "./Results/${out_dir}/${sampleID}/FastQC_post_trim", mode: 'copy'
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

	publishDir "./Results/${out_dir}/${sampleID}/Trimmomatic/", mode: 'copy'
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

trimmomatic_output.into { trimmomatic_output_for_SPAdes ; trimmomatic_output_for_FastQC ; trimmomatic_output_for_fastp ; trimmomatic_output_for_skesa}

process fastqc_post_trimmomatic {
	
	publishDir "./Results/${out_dir}/${sampleID}/FastQC_post_trim", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from trimmomatic_output_for_FastQC

	output:
	file "*" into fastqc_raw
	
	
	"""
	fastqc ${reads[0]} ${reads[1]} --outdir .

	"""
}

process fastp_qc_post_trimmomatic{
	
	publishDir "./Results/${out_dir}/${sampleID}/Fastp_QC/", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from trimmomatic_output_for_fastp 
	
	output:
	file "Reports/*" into fastp_reports_trimmomatic

	when:
	params.trimmomatic == true
	
	"""
	mkdir Reports
	fastp -i ${reads[0]} -I ${reads[1]} -o ${sampleID}_R1_trimmed_fastp.fq.gz -O ${sampleID}_R2_trimmed_fastp.fq.gz
	mv fastp* Reports/
	"""	
}


process spades_no_trim{
	
	publishDir "./Results/${out_dir}/${sampleID}/SPAdes", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from raw_data_for_spades
	val settings from no_trim_settings
	val script_folder from script_spades_no_trim
	val filter_settings from no_trim_spades_filter_setting
	
	output:
	tuple sampleID, file("${sampleID}*") into spades_output
	file "*" into spades_all
	
	when:
	params.no_trim == true
	
	script:
	if ( params.filter_contigs == false )
	
		"""
		spades.py -1 ${reads[0]} -2 ${reads[1]} -o . ${settings}
		mv scaffolds.fasta ./${sampleID}_spades_${settings}_scaffolds.fasta
		"""
	else
	
		"""
		echo ${filter_settings}
		spades.py -1 ${reads[0]} -2 ${reads[1]} -o . ${settings}
		mv scaffolds.fasta ./${sampleID}_spades_${settings}_scaffolds.fasta
		python3 ${script_folder}/remove_contaminants_spades.py ${sampleID}_spades_${settings}_scaffolds ${filter_settings}
		rm ${sampleID}_spades_${settings}_scaffolds.fasta
		"""

}



process spades_after_fastp{
	
	publishDir "./Results/${out_dir}/${sampleID}/SPAdes", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from fastp_output_SPAdes
	val script_folder from script_spades_fastp
	val settings from trim_fastp_setting
	val filter_settings from fastp_spades_filter_setting
	
	output:
	tuple sampleID, file("${sampleID}*") into spades_output_fastp
	file "*" into spades_all_fastp

 	
 	when:
	params.fastp_trim_qc == true
	
	script:
	if ( params.filter_contigs == false )
	
		"""
		spades.py -1 ${reads[0]} -2 ${reads[1]} -o . ${settings}
		mv scaffolds.fasta ./${sampleID}_spades_${settings}_scaffolds.fasta
		"""
	else
	
		"""
		spades.py -1 ${reads[0]} -2 ${reads[1]} -o . ${settings}
		mv scaffolds.fasta ./${sampleID}_spades_${settings}_scaffolds.fasta
		python3 ${script_folder}/remove_contaminants_spades.py ${sampleID}_spades_${settings}_scaffolds ${filter_settings}
		rm ${sampleID}_spades_${settings}_scaffolds.fasta
		"""

}

process spades_after_trimmomatic{
	
	publishDir "./Results/${out_dir}/${sampleID}/SPAdes", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from trimmomatic_output_for_SPAdes
	val script_folder from script_spades_trimmomatic
	val settings from trim_trimmomatic_setting
	val filter_settings from trimmomatic_spades_filter_setting
	
	output:
	tuple sampleID, file("${sampleID}*") into spades_output_trimmomatic
	file "*" into spades_all_trimmomatic

 	when:
	params.trimmomatic == true
	
	script:
	if ( params.filter_contigs == false )
	
		"""
		spades.py -1 ${reads[0]} -2 ${reads[1]} -o . ${settings}
		mv scaffolds.fasta ./${sampleID}_spades_${settings}_scaffolds.fasta
		"""
	else
	
		"""
		spades.py -1 ${reads[0]} -2 ${reads[1]} -o . ${settings}
		mv scaffolds.fasta ./${sampleID}_spades_${settings}_scaffolds.fasta
		python3 ${script_folder}/remove_contaminants_spades.py ${sampleID}_spades_${settings}_scaffolds ${filter_settings}
		rm ${sampleID}_spades_${settings}_scaffolds.fasta
		"""


}

spades_output.into { spades_output_for_quast_no_trim ; spades_output_for_pilon_no_trim ; spades_no_trim_sorting }
spades_output_fastp.into {spades_output_for_quast_fastp ; spades_output_fastp_for_pilon_trim ; spades_fastp_sorting }
spades_output_trimmomatic.into { spades_output_quast_trimmomatic ; spades_output_trimmomatic_for_pilon_trim ; spades_trimmomatic_sorting}


process skesa_no_trim{
	
	publishDir "./Results/${out_dir}/${sampleID}/SKESA", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from raw_data_for_skesa
	val settings from no_trim_settings
	val script_folder from script_skesa_no_trim
	val filter_settings from no_trim_skesa_filter_setting
		
	output:
	tuple sampleID, file("${sampleID}*") into skesa_output_no_trim
	
	when:
	params.no_trim == true
	
	script:
	
	if ( params.filter_contigs == false )
	
		"""
		skesa --reads ${reads[0]},${reads[1]} --use_paired_ends > ${sampleID}_skesa_contigs.fasta
	
		"""
	else
		"""
		skesa --reads ${reads[0]},${reads[1]} --use_paired_ends > ${sampleID}_skesa_contigs.fasta
		python3 ${script_folder}/remove_contaminants.py ${sampleID}_skesa_contigs ${filter_settings}
		rm ${sampleID}_skesa_contigs.fasta
		"""

}

process skesa_after_fastp{
	
	publishDir "./Results/${out_dir}/${sampleID}/SKESA", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from fastp_output_skesa
	val script_folder from script_skesa_fastp
	val filter_settings from fastp_skesa_filter_setting
		
	output:
	tuple sampleID, file("${sampleID}*") into skesa_output_fastp
	
	when:
	params.fastp_trim_qc == true
	
	script:
	
	if ( params.filter_contigs == false )
	
		"""
		skesa --reads ${reads[0]},${reads[1]} --use_paired_ends > ${sampleID}_skesa_contigs.fasta
	
		"""
	else
		"""
		skesa --reads ${reads[0]},${reads[1]} --use_paired_ends > ${sampleID}_skesa_contigs.fasta
		python3 ${script_folder}/remove_contaminants.py ${sampleID}_skesa_contigs ${filter_settings}
		rm ${sampleID}_skesa_contigs.fasta
		"""
}

process skesa_after_trimmomatic{
	
	publishDir "./Results/${out_dir}/${sampleID}/SKESA", mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(reads) from trimmomatic_output_for_skesa
	val script_folder from script_skesa_trimmomatic
	val filter_settings from trimmomatic_skesa_filter_setting
		
	output:
	tuple sampleID, file("${sampleID}*") into skesa_output_trimmomatic
	
	when:
	params.trimmomatic == true
	
	script:
	
	if ( params.filter_contigs == false )
	
		"""
		skesa --reads ${reads[0]},${reads[1]} --use_paired_ends > ${sampleID}_skesa_contigs.fasta
	
		"""
	else
		"""
		skesa --reads ${reads[0]},${reads[1]} --use_paired_ends > ${sampleID}_skesa_contigs.fasta
		python3 ${script_folder}/remove_contaminants.py ${sampleID}_skesa_contigs ${filter_settings}
		rm ${sampleID}_skesa_contigs.fasta
		"""
}

skesa_output_no_trim.into { skesa_output_quast_no_trim ; skesa_no_trim_sorting ; skesa_no_trim_pilon }
skesa_output_fastp.into { skesa_output_quast_fastp ; skesa_fastp_sorting ; skesa_fastp_pilon }
skesa_output_trimmomatic.into { skesa_output_quast_trimmomatic ; skesa_trimmomatic_sorting ; skesa_trimmomatic_pilon }

//join the raw data for bowtie with the spades results by using sampleID as key
spades_raw_data_join = spades_output_for_pilon_no_trim.join(raw_data_for_bowtie2_spades_no_trim)
spades_raw_trimmomatic_join = spades_output_trimmomatic_for_pilon_trim.join(raw_data_for_bowtie2_spades_trimmomatic)
spades_raw_fastp_join = spades_output_fastp_for_pilon_trim.join(raw_data_for_bowtie2_spades_fastp)

//join the raw data for bowtie with skesa results by using sampleID as key
skesa_raw_data_join = skesa_no_trim_pilon.join(raw_data_for_bowtie2_skesa_no_trim)
skesa_raw_trimmomatic_join = skesa_trimmomatic_pilon.join(raw_data_for_bowtie2_skesa_trimmomatic)
skesa_raw_fastp_join = skesa_fastp_pilon.join(raw_data_for_bowtie2_skesa_fastp)

process pilon_post_spades_no_trim {

	publishDir "./Results/${out_dir}/${sampleID}/Pilon/SPAdes", mode: 'copy'
	
	input:
	tuple sampleID, file(scaffold), file(reads) from spades_raw_data_join
	
	output:
	tuple sampleID, file("${sampleID}*") into pilon_output_spades_no_trim
	
	when:
	params.assembly_improvement==true && params.no_trim==true
	
	script:
	index_base = scaffold[0].toString() - ~/.fasta/
	
	"""
	mkdir index
	bowtie2-build ${scaffold[0]} index/$index_base
	bowtie2 -x index/${index_base} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_sorted_alignment.bam
	samtools index ${sampleID}_sorted_alignment.bam
	
	pilon --genome ${scaffold[0]} --frags ${sampleID}_sorted_alignment.bam --output ${index_base}_improved
	
	rm *.bam
	rm *.bai
	"""
	

}

process pilon_post_spades_fastp {

	publishDir "./Results/${out_dir}/${sampleID}/Pilon/SPAdes", mode: 'copy'
	
	input:
	tuple sampleID, file(scaffold), file(reads) from spades_raw_fastp_join
	
	output:
	tuple sampleID, file("${sampleID}*") into pilon_output_spades_fastp
	
	when:
	params.assembly_improvement==true && params.fastp_trim_qc==true
	
	script:
	index_base = scaffold[0].toString() - ~/.fasta/
	
	"""
	mkdir index
	bowtie2-build ${scaffold[0]} index/$index_base
	bowtie2 -x index/${index_base} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_sorted_alignment.bam
	samtools index ${sampleID}_sorted_alignment.bam
	
	pilon --genome ${scaffold[0]} --frags ${sampleID}_sorted_alignment.bam --output ${index_base}_improved
	
	rm *.bam
	rm *.bai
	"""
	

}

process pilon_post_spades_trimmomatic {

	publishDir "./Results/${out_dir}/${sampleID}/Pilon/SPAdes", mode: 'copy'
	
	input:
	tuple sampleID, file(scaffold), file(reads) from spades_raw_trimmomatic_join
	
	output:
	tuple sampleID, file("${sampleID}*") into pilon_output_spades_trimmomatic
	
	when:
	params.assembly_improvement==true && params.trimmomatic==true
	
	script:
	index_base = scaffold[0].toString() - ~/.fasta/
	
	"""
	mkdir index
	bowtie2-build ${scaffold[0]} index/$index_base
	bowtie2 -x index/${index_base} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_sorted_alignment.bam
	samtools index ${sampleID}_sorted_alignment.bam
	
	pilon --genome ${scaffold[0]} --frags ${sampleID}_sorted_alignment.bam --output ${index_base}_improved
	
	rm *.bam
	rm *.bai
	"""
	

}

process pilon_post_skesa_no_trim {

	publishDir "./Results/${out_dir}/${sampleID}/Pilon/SKESA", mode: 'copy'
	
	input:
	tuple sampleID, file(scaffold), file(reads) from skesa_raw_data_join
	
	output:
	tuple sampleID, file("${sampleID}*") into pilon_output_skesa_no_trim
	
	when:
	(params.assembly_improvement==true && params.no_trim==true)
	
	script:
	index_base = scaffold[0].toString() - ~/.fasta/
	
	"""
	mkdir index
	bowtie2-build ${scaffold[0]} index/$index_base
	bowtie2 -x index/${index_base} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_sorted_alignment.bam
	samtools index ${sampleID}_sorted_alignment.bam
	
	pilon --genome ${scaffold[0]} --frags ${sampleID}_sorted_alignment.bam --output ${index_base}_improved
	
	rm *.bam
	rm *.bai
	"""
	

}

process pilon_post_skesa_fastp {

	publishDir "./Results/${out_dir}/${sampleID}/Pilon/SKESA", mode: 'copy'
	
	input:
	tuple sampleID, file(scaffold), file(reads) from skesa_raw_fastp_join
	
	output:
	tuple sampleID, file("${sampleID}*") into pilon_output_skesa_fastp
	
	when:
	params.assembly_improvement==true && params.fastp_trim_qc==true
	
	script:
	index_base = scaffold[0].toString() - ~/.fasta/
	
	"""
	mkdir index
	bowtie2-build ${scaffold[0]} index/$index_base
	bowtie2 -x index/${index_base} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_sorted_alignment.bam
	samtools index ${sampleID}_sorted_alignment.bam
	
	pilon --genome ${scaffold[0]} --frags ${sampleID}_sorted_alignment.bam --output ${index_base}_improved
	
	rm *.bam
	rm *.bai
	"""
	

}

process pilon_post_skesa_trimmomatic {

	publishDir "./Results/${out_dir}/${sampleID}/Pilon/SKESA", mode: 'copy'
	
	input:
	tuple sampleID, file(scaffold), file(reads) from skesa_raw_trimmomatic_join
	
	output:
	tuple sampleID, file("${sampleID}*") into pilon_output_skesa_trimmomatic
	
	when:
	params.assembly_improvement==true && params.trimmomatic==true
	
	script:
	index_base = scaffold[0].toString() - ~/.fasta/
	
	"""
	mkdir index
	bowtie2-build ${scaffold[0]} index/$index_base
	bowtie2 -x index/${index_base} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}_sorted_alignment.bam
	samtools index ${sampleID}_sorted_alignment.bam
	
	pilon --genome ${scaffold[0]} --frags ${sampleID}_sorted_alignment.bam --output ${index_base}_improved
	
	rm *.bam
	rm *.bai
	"""
	

}

process quast_no_trimming{
	
	publishDir "./Results/${out_dir}/${sampleID}/Quast", mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from spades_output_for_quast_no_trim
	tuple sampleID_skesa, file(assembly_skesa) from skesa_output_quast_no_trim
	val reference from ref_for_quast_no_trim
	
	output:
	file "SPAdes/*" into quast_output_spades_no_trim
	file "SKESA/*" into quast_output_skesa_no_trim
	file "Reports/*" into quast_report1

	
	"""
	mkdir Reports
	quast.py ${assembly[0]} -o SPAdes -r $reference
	cp SPAdes/report.tsv Reports/${sampleID}_SPAdes_${out_dir}.tsv
	quast.py ${assembly_skesa[0]} -o SKESA -r $reference
	cp SKESA/report.tsv Reports/${sampleID}_SKESA_${out_dir}.tsv
	"""
}

process quast_after_fastp{
	
	publishDir "./Results/${out_dir}/${sampleID}/Quast", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from spades_output_for_quast_fastp
	tuple sampleID_skesa, file(assembly_skesa) from skesa_output_quast_fastp
	val reference from ref_for_quast_no_trim
	
	output:
	file "SPAdes/*" into quast_output_spades_fastp
	file "SKESA/*" into quast_output_skesa_fastp
	file "Reports/*" into quast_report2
	
	"""
	mkdir Reports
	quast.py ${assembly[0]} -o SPAdes -r $reference
	cp SPAdes/report.tsv Reports/${sampleID}_SPAdes_${out_dir}.tsv
	quast.py ${assembly_skesa[0]} -o SKESA -r $reference
	cp SKESA/report.tsv Reports/${sampleID}_SKESA_${out_dir}.tsv
	"""
}

process quast_after_trimmomatic{
	
	publishDir "./Results/${out_dir}/${sampleID}/Quast", mode: 'copy'
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from spades_output_quast_trimmomatic
	tuple sampleID_skesa, file(assembly_skesa) from skesa_output_quast_trimmomatic
	val reference from ref_for_quast_no_trim
	
	output:
	file "SPAdes/*" into quast_output_spades_trimmomatic
	file "SKESA/*" into quast_output_skesa_trimmomatic
	file "Reports/*" into quast_report3
	
	"""
	mkdir Reports
	quast.py ${assembly[0]} -o SPAdes -r $reference
	cp SPAdes/report.tsv Reports/${sampleID}_SPAdes_${out_dir}.tsv
	quast.py ${assembly_skesa[0]} -o SKESA -r $reference
	cp SKESA/report.tsv Reports/${sampleID}_SKESA_${out_dir}.tsv
	"""
}

process quast_no_trimming_improved{
	
	publishDir "./Results/${out_dir}/${sampleID}/Quast", mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from pilon_output_spades_no_trim
	tuple sampleID_skesa, file(assembly_skesa) from pilon_output_skesa_no_trim
	val reference from ref_for_quast_no_trim
	
	output:
	file "SPAdes/*" into quast_output_spades_improved_no_trim
	file "SKESA/*" into quast_output_skesa_improved_no_trim
	file "Reports/*" into quast_report4
	
	"""
	mkdir Reports
	quast.py ${assembly[0]} -o SPAdes -r $reference
	cp SPAdes/report.tsv Reports/${sampleID}_SPAdes_${out_dir}.tsv
	quast.py ${assembly_skesa[0]} -o SKESA -r $reference
	cp SKESA/report.tsv Reports/${sampleID}_SKESA_${out_dir}.tsv
	"""
}

process quast_after_fastp_improved{
	
	publishDir "./Results/${out_dir}/${sampleID}/Quast", mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from pilon_output_spades_fastp
	tuple sampleID_skesa, file(assembly_skesa) from pilon_output_skesa_fastp
	val reference from ref_for_quast_no_trim
	
	output:
	file "SPAdes/*" into quast_output_spades_improved_fastp
	file "SKESA/*" into quast_output_skesa_improved_fastp
	file "Reports/*" into quast_report5
	
	"""
	mkdir Reports
	quast.py ${assembly[0]} -o SPAdes -r $reference
	cp SPAdes/report.tsv Reports/${sampleID}_SPAdes_${out_dir}.tsv
	quast.py ${assembly_skesa[0]} -o SKESA -r $reference
	cp SKESA/report.tsv Reports/${sampleID}_SKESA_${out_dir}.tsv
	"""
}

process quast_after_trimmomatic_improved{
	
	publishDir "./Results/${out_dir}/${sampleID}/Quast", mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from pilon_output_spades_trimmomatic
	tuple sampleID_skesa, file(assembly_skesa) from pilon_output_skesa_trimmomatic
	val reference from ref_for_quast_no_trim
	
	output:
	file "SPAdes/*" into quast_output_spades_improved_trimmomatic
	file "SKESA/*" into quast_output_skesa_improved_trimmomatic
	file "Reports/*" into quast_report6
	
	"""
	mkdir Reports
	quast.py ${assembly[0]} -o SPAdes -r $reference
	cp SPAdes/report.tsv Reports/${sampleID}_SPAdes_${out_dir}.tsv
	quast.py ${assembly_skesa[0]} -o SKESA -r $reference
	cp SKESA/report.tsv Reports/${sampleID}_SKESA_${out_dir}.tsv
	"""

}



//result.view { it }
