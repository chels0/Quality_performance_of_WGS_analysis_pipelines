#!/usr/bin/env nextflow

genomes = Channel.fromFilePairs(params.path_to_reads).ifEmpty { error "Missing input reads, directory is either empty or missing read pair" }

genomes.into { raw_data_for_Fastp_trim; raw_data_for_Fastp_raw ; raw_data_for_trimmomatic; raw_data_for_spades; raw_data_for_bowtie2_spades_no_trim ; raw_data_for_bowtie2_skesa_no_trim ; raw_data_for_bowtie2_spades_fastp ; raw_data_for_bowtie2_skesa_fastp ; raw_data_for_bowtie2_spades_trimmomatic ; raw_data_for_bowtie2_skesa_trimmomatic ; raw_data_for_fastqc ; raw_data_for_skesa}

reference = Channel.value(params.path_to_reference).ifEmpty { error "No reference genome. Add reference genome to Raw_data directory" }
reference.into { ref_for_quast_no_trim ; ref_for_quast_fastp_trim ; ref_for_quast_trimmomatic_trim }

scripts_folder = Channel.value(params.path_to_scripts)
scripts_folder.into { script_fastqc ; script_spades_no_trim ; script_skesa_no_trim ; script_spades_fastp ; script_skesa_fastp ; script_spades_trimmomatic ; script_skesa_trimmomatic}

trimmomatic_setting = Channel.value(params.trimmomatic_set)
adapters = Channel.value(params.adapter)
minority_adapter = Channel.value(params.not_used_first_adapter)


spades_setting = Channel.value(params.spades_set)
spades_setting.into { no_trim_settings ; trim_trimmomatic_setting ; trim_fastp_setting }
filter_setting = Channel.value(params.filter_set)
filter_setting.into { no_trim_skesa_filter_setting ; fastp_skesa_filter_setting ; trimmomatic_skesa_filter_setting ; no_trim_spades_filter_setting ; fastp_spades_filter_setting ; trimmomatic_spades_filter_setting }


filter = params.filter_set
spades_set = params.spades_set.toUpperCase() - ~/--/
length = spades_set.length()
character = spades_set.substring(0, length - (length-1))
assembler = params.assembler 

out_dir = ""
if ( params.no_trim == true )
	out_dir = out_dir+"N"
if ( params.trimmomatic == true )
	out_dir = out_dir+"T"
if ( params.fastp_trim_qc == true )
	out_dir = out_dir+"F"
if ( params.assembler == 'skesa' )
	out_dir = out_dir+"Ske"

if ( params.assembler == 'spades' )
	out_dir = out_dir+"Sp"+character	
if ( params.filter_contigs == true )
	out_dir = out_dir+filter+"f"
if ( params.assembly_improvement == true )
	out_dir = out_dir + "P"


//else
//	out_dir = ""
//	if ( params.no_trim == true )
//		out_dir = out_dir+"NSp"
//	if ( params.trimmomatic == true )
//		out_dir = out_dir+"TSp"
//	if ( params.fastp_trim_qc == true )
//		out_dir = out_dir+"FSp"
//	if ( params.filter_contigs == true )
//		out_dir = out_dir+filter+"f"
//	if ( params.assembly_improvement == true )
//		out_dir = out_dir + "P"
//	out_dir = out_dir+"_"+spades_set

//if ( params.no_trim == true )
//	out_dir = out_dir+"No_trimming_"+assembler+"_"
//if ( params.trimmomatic == true )
//	out_dir = out_dir+"After_trimmomatic_"+assembler+"_"
//if ( params.fastp_trim_qc == true )
//	out_dir = out_dir+"After_fastp_"+assembler+"_"
//if ( params.assembler == 'spades' )
//	out_dir = out_dir+spades_set+"_"
//if ( params.filter_contigs == true )
//	out_dir = out_dir+filter+"filter_"
//if ( params.assembly_improvement == true )
//	out_dir = out_dir+"Improved_pilon"
println { out_dir }

	
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
	val settings from trimmomatic_setting
	val adapter from adapters
	val lab from minority_adapter
	
	output:
	tuple sampleID, file("*_paired.fq.gz") into trimmomatic_output 
	
	when:
	params.trimmomatic == true
	
	script:
	if ( sampleID.contains(lab[0]) || sampleID.contains(lab[1]) )
		adapters = adapter[1]
	else
		adapters = adapter[0]
	
	"""
	trimmomatic PE ${reads[0]} ${reads[1]} ${sampleID}_R1_paired.fq.gz ${sampleID}_R1_unpaired.fq.gz ${sampleID}_R2_paired.fq.gz ${sampleID}_R2_unpaired.fq.gz ${adapters} ${settings}
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


assembler = params.assembler
process assembly_no_trim{
	
	publishDir "./Results/${out_dir}/${sampleID}/${assembler}", mode: 'copy'
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
	
	if ( params.assembler == 'skesa' )
		if ( params.filter_contigs == false )
			"""
			skesa --reads ${reads[0]},${reads[1]} --use_paired_ends > ${sampleID}_no_trim_skesa_contigs.fasta
		
			"""
		else
			"""
			skesa --reads ${reads[0]},${reads[1]} --use_paired_ends > ${sampleID}_no_trim_skesa_contigs.fasta
			python3 ${script_folder}/remove_contaminants.py ${sampleID}_no_trim_skesa_contigs ${filter_settings}
			rm ${sampleID}_no_trim_skesa_contigs.fasta
			"""
	
	else if ( params.assembler == 'spades')
		if ( params.filter_contigs == false )
	
			"""
			spades.py -1 ${reads[0]} -2 ${reads[1]} -o . ${settings}
			mv scaffolds.fasta ./${sampleID}_no_trim_spades_${settings}_scaffolds.fasta
			"""
		else
	
			"""
			spades.py -1 ${reads[0]} -2 ${reads[1]} -o . ${settings}
			mv scaffolds.fasta ./${sampleID}_no_trim_spades_${settings}_scaffolds.fasta
			python3 ${script_folder}/remove_contaminants_spades.py ${sampleID}_no_trim_spades_${settings}_scaffolds ${filter_settings}
			rm ${sampleID}_no_trim_spades_${settings}_scaffolds.fasta
			"""
			
	
		

}



process assembly_after_fastp{
	
	publishDir "./Results/${out_dir}/${sampleID}/${assembler}", mode: 'copy'
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
	if ( params.assembler == 'skesa' )
		if ( params.filter_contigs == false )
			"""
			skesa --reads ${reads[0]},${reads[1]} --use_paired_ends > ${sampleID}_trim_fastp_skesa_contigs.fasta
		
			"""
		else
			"""
			skesa --reads ${reads[0]},${reads[1]} --use_paired_ends > ${sampleID}_trim_fastp_skesa_contigs.fasta
			python3 ${script_folder}/remove_contaminants.py ${sampleID}_trim_fastp_skesa_contigs ${filter_settings}
			rm ${sampleID}_trim_fastp_skesa_contigs.fasta
			"""
	
	else if ( params.assembler == 'spades')
		if ( params.filter_contigs == false )
	
			"""
			spades.py -1 ${reads[0]} -2 ${reads[1]} -o . ${settings}
			mv scaffolds.fasta ./${sampleID}_trim_fastp_spades_${settings}_scaffolds.fasta
			"""
		else
	
			"""
			spades.py -1 ${reads[0]} -2 ${reads[1]} -o . ${settings}
			mv scaffolds.fasta ./${sampleID}_trim_fastp_spades_${settings}_scaffolds.fasta
			python3 ${script_folder}/remove_contaminants_spades.py ${sampleID}_trim_fastp_spades_${settings}_scaffolds ${filter_settings}
			rm ${sampleID}_trim_fastp_spades_${settings}_scaffolds.fasta
			"""
			

}

process assembly_after_trimmomatic{
	
	publishDir "./Results/${out_dir}/${sampleID}/${assembler}", mode: 'copy'
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
	if ( params.assembler == 'skesa' )
		if ( params.filter_contigs == false )
			"""
			skesa --reads ${reads[0]},${reads[1]} --use_paired_ends > ${sampleID}_trim_trimmomatic_skesa_contigs.fasta
		
			"""
		else
			"""
			skesa --reads ${reads[0]},${reads[1]} --use_paired_ends > ${sampleID}_trim_trimmomatic_skesa_contigs.fasta
			python3 ${script_folder}/remove_contaminants.py ${sampleID}_trim_trimmomatic_skesa_contigs ${filter_settings}
			rm ${sampleID}_trim_trimmomatic_skesa_contigs.fasta
			"""
	
	else if ( params.assembler == 'spades')
		if ( params.filter_contigs == false )
	
			"""
			spades.py -1 ${reads[0]} -2 ${reads[1]} -o . ${settings}
			mv scaffolds.fasta ./${sampleID}_trim_trimmomatic_spades_${settings}_scaffolds.fasta
			"""
		else
	
			"""
			spades.py -1 ${reads[0]} -2 ${reads[1]} -o . ${settings}
			mv scaffolds.fasta ./${sampleID}_trim_trimmomatic_spades_${settings}_scaffolds.fasta
			python3 ${script_folder}/remove_contaminants_spades.py ${sampleID}_trim_trimmomatic_spades_${settings}_scaffolds ${filter_settings}
			rm ${sampleID}_trim_trimmomatic_spades_${settings}_scaffolds.fasta
			"""
		

}

spades_output.into { spades_output_for_quast_no_trim ; spades_output_for_pilon_no_trim ; spades_no_trim_sorting }
spades_output_fastp.into {spades_output_for_quast_fastp ; spades_output_fastp_for_pilon_trim ; spades_fastp_sorting }
spades_output_trimmomatic.into { spades_output_quast_trimmomatic ; spades_output_trimmomatic_for_pilon_trim ; spades_trimmomatic_sorting}




//join the raw data for bowtie with the spades results by using sampleID as key
spades_raw_data_join = spades_output_for_pilon_no_trim.join(raw_data_for_bowtie2_spades_no_trim)
spades_raw_trimmomatic_join = spades_output_trimmomatic_for_pilon_trim.join(raw_data_for_bowtie2_spades_trimmomatic)
spades_raw_fastp_join = spades_output_fastp_for_pilon_trim.join(raw_data_for_bowtie2_spades_fastp)


process pilon_post_spades_no_trim {

	publishDir "./Results/${out_dir}/${sampleID}/Pilon/${assembler}", mode: 'copy'
	
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

	publishDir "./Results/${out_dir}/${sampleID}/Pilon/${assembler}", mode: 'copy'
	
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

	publishDir "./Results/${out_dir}/${sampleID}/Pilon/${assembler}", mode: 'copy'
	
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



process quast_no_trimming{
	
	publishDir "./Results/${out_dir}/${sampleID}/Quast", mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from spades_output_for_quast_no_trim
	val reference from ref_for_quast_no_trim
	
	output:
	file "${assembler}/*" into quast_output_spades_no_trim
	file "Reports/*" into quast_report1

	
	"""
	mkdir Reports
	quast.py ${assembly[0]} -o ${assembler} -r $reference
	cp ${assembler}/report.tsv Reports/${sampleID}_${assembler}_${out_dir}.tsv
	"""
}

process quast_after_fastp{
	
	publishDir "./Results/${out_dir}/${sampleID}/Quast", mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from spades_output_for_quast_fastp
	val reference from ref_for_quast_no_trim
	
	output:
	file "${assembler}/*" into quast_output_spades_fastp
	file "Reports/*" into quast_report2

	
	"""
	mkdir Reports
	quast.py ${assembly[0]} -o ${assembler} -r $reference
	cp ${assembler}/report.tsv Reports/${sampleID}_${assembler}_${out_dir}.tsv
	"""
}

process quast_after_trimmomatic{
	
	publishDir "./Results/${out_dir}/${sampleID}/Quast", mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from spades_output_quast_trimmomatic
	val reference from ref_for_quast_no_trim
	
	output:
	file "${assembler}/*" into quast_output_spades_trimmomatic
	file "Reports/*" into quast_report3

	
	"""
	mkdir Reports
	quast.py ${assembly[0]} -o ${assembler} -r $reference
	cp ${assembler}/report.tsv Reports/${sampleID}_${assembler}_${out_dir}.tsv
	"""
}

process quast_no_trimming_improved{
	
	publishDir "./Results/${out_dir}/${sampleID}/Quast", mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from pilon_output_spades_no_trim
	val reference from ref_for_quast_no_trim
	
	output:
	file "${assembler}/*" into quast_output_spades_improved_no_trim
	file "Reports/*" into quast_report4

	
	"""
	mkdir Reports
	quast.py ${assembly[0]} -o ${assembler} -r $reference
	cp ${assembler}/report.tsv Reports/${sampleID}_${assembler}_${out_dir}.tsv
	"""
}

process quast_after_fastp_improved{
	
	publishDir "./Results/${out_dir}/${sampleID}/Quast", mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from pilon_output_spades_fastp
	val reference from ref_for_quast_no_trim
	
	output:
	file "${assembler}/*" into quast_output_spades_improved_fastp
	file "Reports/*" into quast_report5

	
	"""
	mkdir Reports
	quast.py ${assembly[0]} -o ${assembler} -r $reference
	cp ${assembler}/report.tsv Reports/${sampleID}_${assembler}_${out_dir}.tsv
	"""
}

process quast_after_trimmomatic_improved{
	
	publishDir "./Results/${out_dir}/${sampleID}/Quast", mode: 'copy'
	
	tag "${sampleID}"
	
	input:
	tuple sampleID, file(assembly) from pilon_output_spades_trimmomatic
	val reference from ref_for_quast_no_trim
	
	output:
	file "${assembler}/*" into quast_output_improved_trimmomatic
	file "Reports/*" into quast_report6

	
	"""
	mkdir Reports
	quast.py ${assembly[0]} -o ${assembler} -r $reference
	cp ${assembler}/report.tsv Reports/${sampleID}_${assembler}_${out_dir}.tsv
	"""

}



//result.view { it }
