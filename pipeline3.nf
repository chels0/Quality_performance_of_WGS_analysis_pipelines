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

trimmomatic_setting = Channel.value(params.trimmomatic_set)
spades_setting = Channel.value(params.spades_set)
spades_setting.into { no_trim_settings ; trim_trimmomatic_setting ; trim_fastp_setting }

process fastqc_raw {
	
	publishDir './Results/FastQC_Raw_reads', mode: 'copy'
	
	input:
	tuple sampleID, file(reads) from raw_data_for_fastqc
	val script_folder from scripts_folder

	output:
	file "${sampleID}/*" into fastqc_raw_output
	file "${sampleID}/result.txt" into trimmomatic_input
	
	script:
	if ( params.trimmomatic == true )
	
		"""
		mkdir ${sampleID}
		fastqc ${reads[0]} ${reads[1]} --outdir ${sampleID}
		cd ${sampleID}
		unzip '*.zip'
		cd ${sampleID}_R1_fastqc
		csplit fastqc_data.txt '/>>/' '{*}'
		awk '{ print \$1,\$2 }' xx03 > xx03_relevant
		tail -n +2 xx03_relevant > xx03.csv
		
		echo "\$(python3.6 $script_folder/fastqc_trim_beginning.py)" > lead_R1.txt
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
		
		l_R1="\$(cat lead_R1.txt)"
		l_R2="\$(cat lead_R2.txt)"
		t_R1="\$(cat trail_R1.txt)"
		t_R2="\$(cat trail_R2.txt)"
		
		echo "\$(python3.6 $script_folder/check_trim.py \$l_R1 \$l_R2)" > result.txt
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
	val settings from no_trim_settings
	
	output:
	file "${sampleID}/*_scaffolds.fasta" into spades_output
	file "${sampleID}/*" into spades_all
	val "$sampleID" into spades_id_name 
	
	when:
	params.no_trim == true

	"""
	spades.py -1 ${reads[0]} -2 ${reads[1]} -o $sampleID ${settings}
	mv ${sampleID}/scaffolds.fasta ${sampleID}/${sampleID}_scaffolds.fasta
	"""

}


process spades_after_fastp{
	
	publishDir './Results/'+folder_name+'/SPAdes', mode: 'copy'

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

process skesa_no_trim{
	
	publishDir './Results/No_trimming/SKESA', mode: 'copy'

	input:
	tuple sampleID, file(reads) from raw_data_for_skesa
	val settings from no_trim_settings
	
	output:
	file "${sampleID}/*_contigs.fasta" into skesa_output_no_trim
	val "$sampleID" into skesa_id_name_no_trim 
	
	when:
	params.no_trim == true

	"""
	mkdir ${sampleID}
	skesa --reads ${reads[0]},${reads[1]} --use_paired_ends > ${sampleID}/${sampleID}_contigs.fasta
	"""

}

process skesa_after_fastp{
	
	publishDir './Results/'+folder_name+'/SKESA', mode: 'copy'

	input:
	file sampleID from fastp_output_skesa
	val sample from fastp_id_name_skesa
	
	output:
	file "${sample}/*_contigs.fasta" into skesa_output_fastp
	val "$sample" into skesa_id_name_fastp 
	
	when:
	params.fastp_trim_qc == true

	"""
	dir=\$(echo ${sampleID} | sed 's/_.*//g')
	mkdir \$dir
	skesa --reads ${sampleID[0]},${sampleID[1]} --use_paired_ends > \${dir}/\${dir}_contigs.fasta
	"""
}

process skesa_after_trimmomatic{
	
	publishDir './Results/Trimmed_w_Trimmomatic/SKESA', mode: 'copy'

	input:
	file sampleID from trimmomatic_output_for_skesa
	val sample from id_name_for_skesa
	
	output:
	file "${sample}/*_contigs.fasta" into skesa_output_trimmomatic
	val "${sample}" into skesa_id_name_trimmomatic
	
	when:
	params.trimmomatic == true

	"""
	mkdir ${sample}
	skesa --reads ${sampleID[0]},${sampleID[1]} --use_paired_ends > ${sample}/${sample}_contigs.fasta
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
	file scaffold_skesa from skesa_output_no_trim
	val reference from ref_for_quast_no_trim
	val sampleID from spades_id_name
	
	output:
	file "SPAdes/$sampleID" into quast_output
	file "SKESA/$sampleID" into quast_output_skesa
	val "$sampleID" into quast_id_name
	
	"""
	mkdir SPAdes
	mkdir SKESA
	quast.py $scaffold -o SPAdes/$sampleID -r $reference
	quast.py $scaffold_skesa -o SKESA/$sampleID -r $reference
	"""
}

process quast_after_fastp{
	
	publishDir './Results/'+folder_name+'/Quast', mode: 'copy'

	input:
	file scaffold from spades_output_for_quast_trim
	file scaffold_skesa from skesa_output_fastp
	val reference from ref_for_quast_fastp_trim
	val sampleID from spades_id_name_fastp 
	
	output:
	file "SPAdes/$sampleID" into quast_output_fastp
	file "SKESA/$sampleID" into quast_output_fastp_skesa
	val "$sampleID" into quast_id_name_fastp
	
	"""
	mkdir SPAdes
	mkdir SKESA
	quast.py $scaffold -o SPAdes/$sampleID -r $reference
	quast.py $scaffold_skesa -o SKESA/$sampleID -r $reference
	"""
}

process quast_after_trimmomatic{
	
	publishDir './Results/Trimmed_w_Trimmomatic/Quast', mode: 'copy'

	input:
	file scaffold from spades_output_trimmomatic
	file scaffold_skesa from skesa_output_trimmomatic
	val reference from ref_for_quast_trimmomatic_trim
	val sampleID from spades_id_name_trimmomatic
	
	output:
	file "SPAdes/$sampleID" into quast_output_trimmomatic
	file "SKESA/$sampleID" into quast_output_trimmomatic_skesa
	val "$sampleID" into quast_id_name_trimmomatic
	
	"""
	mkdir SPAdes
	mkdir SKESA
	quast.py $scaffold -o SPAdes/$sampleID -r $reference
	quast.py $scaffold_skesa -o SKESA/$sampleID -r $reference
	"""
}


//result.view { it }
