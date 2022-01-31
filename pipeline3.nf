#!/usr/bin/env nextflow


genomes = Channel.fromFilePairs(params.path_to_reads)
genomes.into { raw_data_for_Fastp; raw_data_for_trimmomatic; raw_data_for_spades; raw_data_for_bowtie2 ; raw_data_for_fastqc}

reference = Channel.value(params.path_to_reference)
reference.into { ref_for_quast_no_trim ; ref_for_quast_fastp_trim ; ref_for_quast_trimmomatic_trim }

reference_index = Channel.value(params.path_to_bowtie2_index)

if ( params.fastp_trim_qc == false )
	folder_name = 'No_trimming'
	
else if ( params.fastp_trim_qc  == true )
	folder_name = 'Trimmed_w_fastp'

folder = Channel.from(folder_name)


process fastqc_raw {
	
	publishDir './Results/FastQC_Raw_reads', mode: 'copy'
	
	input:
	tuple sampleID, file(reads) from raw_data_for_fastqc
	
	output:
	file "${sampleID}/*" into fastqc_raw_output
	
	
	"""
	mkdir ${sampleID}
	fastqc ${reads[0]} ${reads[1]} --outdir ${sampleID}

	"""
}
	
process fastp{
	
	publishDir './Results/'+folder_name+'/Fastp/', mode: 'copy'
	
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

fastp_output.into { fastp_output_fastqc ; fastp_output_SPAdes }
fastp_id_name.into { fastp_id_name_fastqc ; fastp_id_name_SPAdes }

process fastqc_post_fastp_trim {
	
	publishDir './Results/${folders[0]/FastQC_post_trim', mode: 'copy'
	
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
	
	input:
	tuple sampleID, file(reads) from raw_data_for_trimmomatic 
	
	output:
	file "${sampleID}/*_paired.fq.gz" into trimmomatic_output 
	val "${sampleID}" into id_name_trimmomatic
	
	when:
	params.trimmomatic == true
	
	"""
	mkdir ${sampleID}
	cd ${sampleID}
	trimmomatic PE ../${reads[0]} ../${reads[1]} ${sampleID}_R1_paired.fq.gz ${sampleID}_R1_unpaired.fq.gz ${sampleID}_R2_paired.fq.gz ${sampleID}_R2_unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	"""

}

trimmomatic_output.into { trimmomatic_output_for_SPAdes ; trimmomatic_output_for_FastQC }
id_name_trimmomatic.into { id_name_for_SPAdes ; id_name_for_FastQC }

process fastqc_post_trimmomatic {
	
	publishDir './Results/Trimmed_w_Trimmomatic/FastQC_post_trim', mode: 'copy'
	
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

	input:
	tuple sampleID, file(reads) from raw_data_for_spades
	
	output:
	file "${sampleID}/*_scaffolds.fasta" into spades_output
	file "${sampleID}/*" into spades_all
	val "$sampleID" into spades_id_name 
	
	when:
	params.no_trim == true

	"""
	spades.py -1 ${reads[0]} -2 ${reads[1]} --only-assembler -o $sampleID
	mv ${sampleID}/scaffolds.fasta ${sampleID}/${sampleID}_scaffolds.fasta
	"""

}


process spades_after_fastp{
	
	publishDir './Results/${folders[0]}/SPAdes', mode: 'copy'

	input:
	file sampleID from fastp_output_SPAdes
	val sample from fastp_id_name_SPAdes

	output:
	file "${sample}/*_scaffolds.fasta" into spades_output_fastp
	file "${sample}/*" into spades_all_fastp
	val "$sample" into spades_id_name_fastp 
 	
 	when:
 	params.fastp_trim_qc == true
 	
 	"""
	spades.py -1 ${sampleID[0]} -2 ${sampleID[1]} --only-assembler -o $sample
	mv ${sample}/scaffolds.fasta ${sample}/${sample}_trimmed_fastp_scaffolds.fasta
	"""

}

spades_output.into { spades_output_for_quast_no_trim ; spades_output_for_pilon_no_trim }
spades_output_fastp.into {spades_output_for_quast_trim ; spades_output_fastp_for_pilon_trim }

process spades_after_trimmomatic{
	
	publishDir './Results/Trimmed_w_Trimmomatic/SPAdes', mode: 'copy'

	input:
	file sampleID from trimmomatic_output_for_SPAdes
	val sample from id_name_for_SPAdes

	output:
	file "${sample}/*_scaffolds.fasta" into spades_output_trimmomatic
	file "${sample}/*" into spades_all_trimmomatic
	val "$sample" into spades_id_name_trimmomatic
 	
 	when:
 	params.trimmomatic==true

	"""
	spades.py -1 ${sampleID[0]} -2 ${sampleID[1]} --only-assembler -o $sample 
	mv ${sample}/scaffolds.fasta ${sample}/${sample}_trimmed_trimmomatic_scaffolds.fasta
	"""


}


process bowtie2_no_trimming {

	publishDir './Results/No_trimming/Bowtie2', mode: 'copy'
	
	input:
	val reference from reference_index.collect()
	tuple sampleID, file(reads) from raw_data_for_bowtie2
	
	output:
	file "${sampleID}/*" into bowtie_all
	val "${sampleID}" into bowtie_id_name
	stdout result2
	
	when:
	params.assembly_improvement == true
	
	script:
	// Tar bort allt som har med bowtie index att göra på ett av namnen så endast namnet är med
	index_base = reference[0].toString() - ~/.rev.\d.bt2|.\d.bt2/
	
	"""
	echo $reference
	echo ${reference[0]}
	mkdir ${sampleID}
	bowtie2 -x ${index_base} -1 ${reads[0]} -2 ${reads[1]} -p 4 | samtools sort -@ 4 -o ${sampleID}/${sampleID}_sorted_alignment.bam
	samtools index ${sampleID}/${sampleID}_sorted_alignment.bam
	"""
}

process pilon_no_trimming {

	publishDir './Results/No_trimming/Pilon', mode: 'copy'
	
	input:
	file scaffold from spades_output_for_pilon_no_trim
	file alignment from bowtie_all.collect()
	val sampleID from bowtie_id_name
	
	
	output:
	file "${sampleID}/*" into pilon_output
	stdout result
	
	when:
	params.assembly_improvement == true
	
	script:
	// ta bort bam ändelsen
	bam_file = alignment[0].toString() - ~/.bam*/
	
	"""
	echo $bam_file
	mkdir ${sampleID}
	pilon --genome ${scaffold} --frags ${bam_file}.bam --output ${sampleID}/${sampleID}
	"""

}


process quast_no_trimming{
	
	publishDir './Results/No_trimming/Quast', mode: 'copy'

	input:
	file scaffold from spades_output_for_quast_no_trim
	val reference from ref_for_quast_no_trim
	val sampleID from spades_id_name
	
	output:
	file "$sampleID" into quast_output
	val "$sampleID" into quast_id_name
	
	"""
	quast.py $scaffold -o $sampleID -r $reference
	"""
}

process quast_after_fastp{
	
	publishDir './Results/Trimmed_w_Fastp/Quast', mode: 'copy'

	input:
	file scaffold from spades_output_for_quast_trim
	file reference from ref_for_quast_fastp_trim
	val sampleID from spades_id_name_fastp 
	
	output:
	file "$sampleID" into quast_output_fastp
	val "$sampleID" into quast_id_name_fastp
	
	"""
	quast.py $scaffold -o $sampleID -r ${reference}
	"""
}

process quast_after_trimmomatic{
	
	publishDir './Results/Trimmed_w_Trimmomatic/Quast', mode: 'copy'

	input:
	file scaffold from spades_output_trimmomatic
	file reference from ref_for_quast_trimmomatic_trim
	val sampleID from spades_id_name_trimmomatic
	
	output:
	file "$sampleID" into quast_output_trimmomatic
	val "$sampleID" into quast_id_name_trimmomatic
	
	"""
	quast.py $scaffold -o $sampleID -r ${reference}
	"""
}


result.view { it }
