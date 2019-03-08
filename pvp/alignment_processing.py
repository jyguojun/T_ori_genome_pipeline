#!/usr/bin/env python

import os
from pvp_utils import *

def sort_sam(sample_name, temp_directory, out): # To Sort a SAM or BAM file by "COORDINATE"
	#java -Xmx2g -XX:+UseSerialGC -jar $picardtools SortSam INPUT=$name.bam OUTPUT=$name.sorted.bam SORT_ORDER=coordinate
	output_ngm_bam = "%s.ngm.bam" %sample_name
	output_bwa_bam = "%s.bwa.bam" %sample_name
	output_ngm_sorted_bam = "%s.ngm.sorted.bam" %sample_name
	output_bwa_sorted_bam = "%s.bwa.sorted.bam" %sample_name
	INPUT_ngm = "INPUT=%s" %output_ngm_bam
	INPUT_bwa = "INPUT=%s" %output_bwa_bam
	OUTPUT_ngm = "OUTPUT=%s" %output_ngm_sorted_bam
	OUTPUT_bwa = "OUTPUT=%s" %output_bwa_sorted_bam
	if os.path.exists(os.path.join(out, output_ngm_sorted_bam)) and os.path.exists(os.path.join(out, output_bwa_sorted_bam)):
		return
	log = os.path.join(out,"genome_log.txt")
	with cd(temp_directory):
		copy_from_out(output_ngm_bam, out)
		print "%s: sortsam for ngm" %sample_name
		run_command(["picard", 
		"-Xmx2g", "-XX:+UseSerialGC",
		"SortSam", INPUT_ngm, OUTPUT_ngm, "SORT_ORDER=coordinate"], log)
		print "%s: sortsam for ngm successfully completed" %sample_name
		copy_to_results_dir(sample_name, output_ngm_sorted_bam, out)
		print "%s: sortsam for bwa" %sample_name
		copy_from_out(output_bwa_bam, out)
		run_command(["picard", 
		"-Xmx2g", "-XX:+UseSerialGC",
		"SortSam", INPUT_bwa, OUTPUT_bwa, "SORT_ORDER=coordinate"], log)
		print "%s: sortsam for bwa successfully completed" %sample_name
		copy_to_results_dir(sample_name, output_bwa_sorted_bam, out)


def mark_duplicates(sample_name, temp_directory, out): # To identify duplicate reads and indicate number of duplicates for both single and paired-end reads
	#java -Xmx2g -XX:+UseSerialGC -jar $picardtools MarkDuplicates INPUT=$name.sorted.bam OUTPUT=$name.dedup.bam METRICS_FILE=$name.dedup.metrics
	output_ngm_sorted_bam = "%s.ngm.sorted.bam" %sample_name
	output_bwa_sorted_bam = "%s.bwa.sorted.bam" %sample_name
	output_ngm_dedup_bam = "%s.ngm.dedup.bam" %sample_name
	output_bwa_dedup_bam = "%s.bwa.dedup.bam" %sample_name
	metrics_ngm = "%s.ngm.dedup.metrics" %sample_name
	metrics_bwa = "%s.bwa.dedup.metrics" %sample_name
	INPUT_ngm = "INPUT=%s" %output_ngm_sorted_bam
	INPUT_bwa = "INPUT=%s" %output_bwa_sorted_bam
	OUTPUT_ngm = "OUTPUT=%s" %output_ngm_dedup_bam
	OUTPUT_bwa = "OUTPUT=%s" %output_bwa_dedup_bam
	METRICS_ngm = "METRICS_FILE=%s" %metrics_ngm
	METRICS_bwa = "METRICS_FILE=%s" %metrics_bwa
	if os.path.exists(os.path.join(out, output_ngm_dedup_bam)) and os.path.exists(os.path.join(out, output_bwa_dedup_bam)):
		return
	log = os.path.join(out,"genome_log.txt")
	with cd(temp_directory):
		copy_from_out(output_ngm_sorted_bam, out)
		print "%s: marking duplicates for ngm" %sample_name
		run_command(["picard", 
		"-Xmx2g", "-XX:+UseSerialGC",
		"MarkDuplicates", INPUT_ngm, OUTPUT_ngm, METRICS_ngm,], log)
		print "%s: duplicates for ngm marked successfully" %sample_name
		copy_to_results_dir(sample_name, output_ngm_dedup_bam, out)
		copy_from_out(output_bwa_sorted_bam, out)
		print "%s: marking duplicates for bwa" %sample_name
		run_command(["picard", 
		"-Xmx2g", "-XX:+UseSerialGC",
		"MarkDuplicates", INPUT_bwa, OUTPUT_bwa, METRICS_bwa,], log)
		print "%s: duplicates for bwa marked successfully" %sample_name
		copy_to_results_dir(sample_name, output_bwa_dedup_bam, out)


def samtools_mpileup(sample_name, temp_directory, out): # To call Snps and short Indels
	output_bwa_dedup_bam = "%s.bwa.dedup.bam" %sample_name
	output_ngm_dedup_bam = "%s.ngm.dedup.bam" %sample_name
	output_bwa_pileup = "%s.bwa.pileup" %sample_name
	output_ngm_pileup = "%s.ngm.pileup" %sample_name
	if os.path.exists(os.path.join(out, output_bwa_pileup)) and os.path.exists(os.path.join(out, output_ngm_pileup)):
		return
	log = os.path.join(out,"genome_log.txt")
	with cd(temp_directory):
		copy_from_out(output_bwa_dedup_bam, out)
		print "%s: samtools mpileup for bwa" %sample_name
		run_po_command(["samtools", "mpileup", "-f", "ref.fa", output_bwa_dedup_bam], output_bwa_pileup, log)
		print "%s: samtools mpileup for bwa successfully completed" %sample_name
		copy_to_results_dir(sample_name, output_bwa_pileup, out)
		copy_from_out(output_ngm_dedup_bam, out)
		print "%s: samtools mpileup for ngm" %sample_name
		run_po_command(["samtools", "mpileup", "-f", "ref.fa", output_ngm_dedup_bam], output_ngm_pileup, log)
		print "%s: samtools mpileup for ngm successfully completed" %sample_name
		copy_to_results_dir(sample_name, output_ngm_pileup, out)
