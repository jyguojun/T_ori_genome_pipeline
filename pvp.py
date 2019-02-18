#!/usr/bin/env python

import sys
import os
import gzip
import subprocess
from contextlib import contextmanager
import shutil
from Bio import SeqIO
import vcf
import re
import multiprocessing
import math
from collections import defaultdict
import csv
import cPickle as pickle

@contextmanager
def cd(newdir):
	prevdir = os.getcwd()
	os.chdir(os.path.expanduser(newdir))
	try:
		yield
	finally:
		os.chdir(prevdir)

def get_sample_name(path):
	filename = os.path.basename(path.rstrip())
	sample_name, ext = os.path.splitext(filename)
	if ext == ".fastq" or ext == ".fq":
		return sample_name, False
	elif ext == ".gz":
		sample_name, ext = os.path.splitext(sample_name)
		if ext == ".fastq" or ext == ".fq":
			return sample_name, True
		else:
			print "Reads file: %s not recogized as fastq! Please ensure this is a fastq file and has the extenstion .fastq or .fq" % os.path.basename(path)
			sys.exit(1)
	else:
		print "Reads file: %s not recogized as fastq! Please ensure this is a fastq file and has the extenstion .fastq or .fq" % os.path.basename(path)
		sys.exit(1)

def unzip_reads(reads, gzipped):
	if gzipped == True:
		sample_name, null = get_sample_name(reads)
		print "%s: Reads file is gzipped. Unzipping..." % sample_name
		unzipped_reads, null = os.path.splitext(reads)
		with gzip.open(reads, 'rb') as f_in, open(unzipped_reads, 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
		print "%s: Unzipping completed." % sample_name
	else:
		unzipped_reads = reads
	return unzipped_reads

def run_command(command, outputfile):
	process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	while True:
		output = process.stdout.readline()
		if output == '' and process.poll() is not None:
			break
		if output:
			with open(outputfile,'a+') as outhandle:
				outhandle.write(output)
	if process.poll() == 1:
		print "\n\n###########    Error    ###########\n\nThe following command failed to complete successfully:\n\n%s\n\nOutput of this error can be found in:\n\n%s\n\n" % (" ".join(command), outputfile)
		sys.exit(1)


def run_po_command(command, outputfile, errorfile):
	process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output, error = stdout, stderr = process.communicate()
	if output:
		with open(outputfile,'w') as outhandle:
			outhandle.write(output)
	if error:
		with open(errorfile,'a+') as errorhandle:
			errorhandle.write(error)
	if process.returncode == 1:
		print "\n\n###########    Error    ###########\n\nThe following command failed to complete successfully:\n\n%s\n\nOutput of this error can be found in:\n\n%s\n\n" % (" ".join(command), errorfile)
		sys.exit(1)


def index_reference(ref, host_ref, temp_directory): #indexing to create samtools, bwa, picard and bowtie2 (cattle host genome) index
	ref_gb = SeqIO.parse(ref, "gb")
	SeqIO.write(ref_gb, os.path.join(temp_directory, "ref.fa"), "fasta")
	# host_ref_gb = SeqIO.parse(host_ref, "gb") # removed as host genome is already in fasta format
	# SeqIO.write(host_ref_gb, os.path.join(temp_directory, "host_ref.fa"), "fasta")
	log = os.path.join(out,"genome_log.txt")
	with cd(temp_directory): 
		run_command(["bwa", "index", "ref.fa" ], log)
		print "Created bwa index"
		run_command(["samtools", "faidx", "ref.fa"], log)
		print "Created samtools index"
		run_command(["picard", "-Xmx2g", "-XX:+UseSerialGC", "CreateSequenceDictionary", "REFERENCE=ref.fa", "OUTPUT=ref.dict" ], log)
		print "Created picard index"
		run_command(["bowtie2-build", host_ref, "host_ref_db"], log)
		print "Created bowtie2 cattle host genome index"


def copy_from_out(filename, out):
	if os.path.exists(filename) == False:
		shutil.copyfile(os.path.join(out, filename),filename)


def trim_reads(reads, sample_name, temp_directory, out, gzipped):  # Using trimmomatic to trim reads 
	paired_output1 = "%s.trimmed_paired1.fastq" % sample_name
	paired_output2 = "%s.trimmed_paired2.fastq" % sample_name
	if os.path.exists(os.path.join(out, paired_output1)) and os.path.exists(os.path.join(out, paired_output2)):
		return
	unpaired_output1 = "%s.trimmed_unpaired1.fastq" % sample_name
	unpaired_output2 = "%s.trimmed_unpaired2.fastq" % sample_name
	merged_trimmed = "%s.trimmed.fastq.gz" % sample_name
	trim_log = "%s.trimlog" % sample_name
	reads1 = "%s.1.fastq" % sample_name
	reads2 = "%s.2.fastq" % sample_name
	print "%s: Processing reads file for trimming" % sample_name 
	log = os.path.join(out,"genome_log.txt")
	with cd(temp_directory):
		print "%s: unmerge reads with reformat.sh " % sample_name
		run_command(["reformat.sh", "in=%s" % reads, "out=%s" % reads1, "out2=%s" % reads2], log)
		print "%s: unmerge complete." % sample_name
		print "%s: Trimming reads" % sample_name
		run_command(["trimmomatic",
		"-Xmx2g","-XX:+UseSerialGC",
		"PE","-phred33",
		"-threads","1",
		"-trimlog",trim_log,
		reads1, reads2,
		paired_output1,unpaired_output1,
		paired_output2,unpaired_output2, "ILLUMINACLIP:/shared/homes/12436414/parasite_pipeline/adapter.fasta:2:30:10",
		"LEADING:20","TRAILING:20",
		"SLIDINGWINDOW:4:15","MINLEN:70"],log)
		print "%s: Read trimming complete." % (sample_name)
		copy_to_results_dir(sample_name, paired_output1, out)
		copy_to_results_dir(sample_name, paired_output2, out)

# def biobloomcategorizer(reads, sample_name, temp_directory, out): #biobloomcategorizer to remove host (Bos taurus) DNA from trimmed reads
	# paired_output1 = "%s.trimmed_paired1.fastq" % sample_name
	# paired_output2 = "%s.trimmed_paired2.fastq" % sample_name
	# output_prefix = "%s.bbc" % sample_name
	# bbc_output1 = "%s_noMatch_1.fq" % output_prefix
	# bbc_output2 = "%s_noMatch_2.fq" % output_prefix
	# if os.path.exists(os.path.join(out, bbc_output1)) and os.path.exists(os.path.join(out, bbc_output2)):
		# return
	# log = os.path.join(out,"genome_log.txt")
	# with cd(temp_directory):
		# copy_from_out(paired_output1, out)
		# copy_from_out(paired_output2, out)
		# print "%s: removing host DNA" % sample_name 
		# run_command(["biobloomcategorizer",
		# "--fq","-e",
		# "-p", output_prefix,
		# "-f", sys.argv[6],
		# paired_output1, paired_output2], log)
		# print "%s: host DNA removal complete" % (sample_name)
		# copy_to_results_dir(sample_name, bbc_output1, out)
		# copy_to_results_dir(sample_name, bbc_output2, out)


def bowtie2_mapping(reads, sample_name, temp_directory, out): #bowtie2 mapping against host sequence
	paired_output1 = "%s.trimmed_paired1.fastq" % sample_name
	paired_output2 = "%s.trimmed_paired2.fastq" % sample_name
	output_prefix = "%s.mapped_and_unmapped.sam" % sample_name # contains both parasite and host reads in sam format.
	output_bam = "%s.mapped_and_unmapped.bam" % sample_name # contains both parasite and host reads in bam format.
	both_unmapped_bam = "%s.both_ends_unmapped.bam" % sample_name #reads that are not mapped to host genome (The parasite seqeunces that we want!)
	both_unmapped_sorted = "%s.both_ends_unmapped.sorted" % sample_name #sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
	both_unmapped_sorted_bam = "%s.both_ends_unmapped.sorted.bam" % sample_name #convert output from "samtool sort" to give it bam format to read into next command
	host_read_removed1 = "%s.host_read_removed_r1.fastq" % sample_name # read 1 of host removed sample
	host_read_removed2 = "%s.host_read_removed_r2.fastq" % sample_name # read 2 of host removed sample
	host_read_removed_merged = "%s.host_read_removed_merged.fastq" %sample_name # merged fastq file of host removed sample
	if os.path.exists(os.path.join(out, output_bam)) and os.path.exists(os.path.join(out, both_unmapped_bam)) and os.path.exists(os.path.join(out, both_unmapped_sorted_bam)) and os.path.exists(os.path.join(out, output_prefix)) and os.path.exists(os.path.join(out, host_read_removed1)) and os.path.exists(os.path.join(out, host_read_removed2)):
		return
	log = os.path.join(out,"genome_log.txt")
	with cd(temp_directory):
		copy_from_out(paired_output1, out)
		copy_from_out(paired_output2, out)
		print "%s:bowtie2 mapping seqeunces to host." %sample_name
		run_command(["bowtie2", "-x", "host_ref_db", "-1", paired_output1, "-2", paired_output2, "-S", output_prefix], log) #bowtie2 mapping against host sequence database, keep both mapped and unmapped reads (paired-end reads)
		print "%s: bowtie2 mapping with host sequence completed" %sample_name
		copy_to_results_dir(sample_name, output_prefix, out)
		print "%s: samtools converting mapped and unmapped sequences from .sam format to .bam" %sample_name
		copy_from_out(output_prefix, out)
		run_po_command(["samtools", "view", "-bS", output_prefix], output_bam, log) #convert file .sam to .bam
		print "%s: conversion to .bam completed" %sample_name
		copy_to_results_dir(sample_name, output_bam, out)
		print "%s:samtools filter unmapped reads" %sample_name #SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
		copy_from_out(output_bam, out)
		run_po_command(["samtools", "view", "-b", "-f", "12", "-F", "256", output_bam], both_unmapped_bam, log)
		print "%s:samtools filter of unmapped reads (parasite reads) completed" %sample_name
		copy_to_results_dir(sample_name, both_unmapped_bam, out)
		print "%s:samtools sort unmapped reads.bam (parasite reads)" %sample_name
		copy_from_out(both_unmapped_bam, out)
		run_command(["samtools", "sort", "-n", "-o", both_unmapped_sorted_bam, both_unmapped_bam], log)
		print "%s:samtools sort unmapped reads.bam (parasite reads) completed" %sample_name
		copy_to_results_dir(sample_name, both_unmapped_sorted_bam, out)
		print "%s:bedtools convert sorted bam files to fastq" %sample_name
		copy_from_out(both_unmapped_sorted_bam, out)
		run_command(["bedtools", "bamtofastq", "-i", both_unmapped_sorted_bam, "-fq", host_read_removed1, "-fq2", host_read_removed2], log)
		print "%s:bedtools conversion of sorted bam files to fastq completed" %sample_name
		copy_to_results_dir(sample_name, host_read_removed1, out)
		copy_to_results_dir(sample_name, host_read_removed2, out)
		copy_from_out(host_read_removed1, out)
		copy_from_out(host_read_removed2, out)
		print "%s: Merging host read removed fastq files." %sample_name
		run_command(["reformat.sh", "in1=%s" %host_read_removed1, "in2=%s" %host_read_removed2, "out=%s" %host_read_removed_merged], log)
		copy_to_results_dir(sample_name, host_read_removed_merged, out)
		copy_from_out(host_read_removed_merged, out)
		run_command(["seqkit", "stats", host_read_removed_merged], log)
		print "%s: Merging of host removed fastq file completed." %sample_name


def bwa(sample_name, temp_directory, out): # using bwa (burrows-wheeler aligner) to align sequence data with reference seq
	host_read_removed1 = "%s.host_read_removed_r1.fastq" % sample_name
	host_read_removed2 = "%s.host_read_removed_r2.fastq" % sample_name
	output_sam = "%s.bwa.sam" % sample_name
	output_bam = "%s.bwa.bam" % sample_name
	if os.path.exists(os.path.join(out, output_bam)):
		return
	log = os.path.join(out,"genome_log.txt")
	with cd(temp_directory):
		copy_from_out(host_read_removed1, out)
		copy_from_out(host_read_removed2, out)
		print "%s: aligning sequences with bwa" %sample_name
		run_po_command(["bwa",
		"mem", "ref.fa",
		host_read_removed1, host_read_removed2], output_sam, log)
		run_po_command(["samtools", "view", "-bh", "-q", "10", output_sam], output_bam, log)
		print "%s: sequence alignment with bwa completed" %sample_name
		os.remove(output_sam)
		copy_to_results_dir(sample_name, output_bam, out)

def ngm(sample_name, temp_directory, out): # using ngm (next gen map) to align sequence data with reference seq
	#$ngm -r ref.fa -1 $name.ec1.fastq -2 $name.ec2.fastq -o $name.pair.sam -p
	host_read_removed1 = "%s.host_read_removed_r1.fastq" % sample_name
	host_read_removed2 = "%s.host_read_removed_r2.fastq" % sample_name
	output_sam = "%s.ngm.sam" % sample_name
	output_bam = "%s.ngm.bam" % sample_name
	if os.path.exists(os.path.join(out, output_bam)):
		return
	log = os.path.join(out,"genome_log.txt")
	with cd(temp_directory):
		copy_from_out(host_read_removed1, out)
		copy_from_out(host_read_removed2, out)
		print "%s: aligning sequences with ngm" %sample_name
		run_command(["ngm",
		"-r", "ref.fa", "-1", host_read_removed1, "-2", host_read_removed2, "-o", output_sam, "-p"], log)
		run_po_command(["samtools", "view", "-bh", "-q", "10", output_sam], output_bam, log)
		print "%s:sequence alignment with ngm completed" %sample_name
		os.remove(output_sam)
		copy_to_results_dir(sample_name, output_bam, out)


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
	#samtools mpileup -f /home/jerald/Documents/jerald_genome/theileria_orientalis_shintoku.fa $name.dedup.bam > $name.pileup
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


def varscan(sample_name, temp_directory, out): # To call variants that meet desired threshhold
	#java -Xmx2g -XX:+UseSerialGC -jar $varscan mpileup2snp $name.pileup --min-coverage 10 --p-value 0.05 --min-avg-qual 20 --min-reads2 4 --min-var-freq 0.01 --min-freq-for-hom 0.9 --strand-filter 0 --output-vcf 1 > $name.raw.vcf 
	output_bwa_pileup = "%s.bwa.pileup" %sample_name
	output_ngm_pileup = "%s.ngm.pileup" %sample_name
	output_bwa_raw_vcf = "%s.bwa.raw.vcf" %sample_name
	output_ngm_raw_vcf = "%s.ngm.raw.vcf" %sample_name
	if os.path.exists(os.path.join(out, output_bwa_raw_vcf)) and os.path.exists(os.path.join(out, output_ngm_raw_vcf)):
		return
	log = os.path.join(out,"genome_log.txt")
	with cd(temp_directory):
		copy_from_out(output_bwa_pileup, out)
		print "%s: counting variants with varscan for bwa" %sample_name
		run_po_command(["varscan", "-Xmx2g", "-XX:+UseSerialGC", 
		"mpileup2snp", output_bwa_pileup, "-p-value", "0.05", "--min-avg-qual", "25", "--min-reads2", "2", "--strand-filter", "1",
		"--output-vcf", "1"], output_bwa_raw_vcf, log)
		print "%s: varscan for bwa successfully completed" %sample_name
		copy_to_results_dir(sample_name, output_bwa_raw_vcf, out)
		copy_from_out(output_ngm_pileup, out)
		print "%s: counting variants with varscan for ngm" %sample_name
		run_po_command(["varscan", "-Xmx2g", "-XX:+UseSerialGC", 
		"mpileup2snp", output_ngm_pileup, "-p-value", "0.05", "--min-avg-qual", "25", "--min-reads2", "2", "--strand-filter", "1",
		"--output-vcf", "1"], output_ngm_raw_vcf, log)
		print "%s: varscan for ngm successfully completed" %sample_name
		copy_to_results_dir(sample_name, output_ngm_raw_vcf, out)

def create_empty_rows(snp_info_d, chrom, position, sample_name, ref_allele, snp_sample_info_d): # to create empty dictonaries for data entry from samples
	if not (chrom, position) in snp_info_d["pass"]:
		snp_info_d["chromosome"][chrom,position] = chrom
		snp_info_d["position"][chrom,position] = position
		snp_info_d["files"][chrom,position] = []
		snp_info_d["total_depth"][chrom,position] = 0
		snp_info_d["alt_depth"][chrom,position] = 0
		snp_info_d["ten_reads"][chrom,position] = False
		snp_info_d["pass"][chrom,position] = True
		snp_info_d["biallelic"][chrom,position] = False
		snp_info_d["cumulative_allele_freq"][chrom,position] = 0
		snp_info_d["multiple_allele"][chrom,position] = []
		snp_info_d["fail_description"][chrom,position] = []
		snp_info_d["ref_allele"][chrom,position] = []
	if not (chrom, position, sample_name) in snp_sample_info_d["heterozygous"]:
		snp_sample_info_d["heterozygous"][chrom,position,sample_name] = False #3 item tuple
		snp_sample_info_d["sample_chromosome"][chrom,position,sample_name] = chrom
		snp_sample_info_d["sample_position"][chrom,position,sample_name] = position
		snp_sample_info_d["sample_name"][chrom,position,sample_name] = sample_name
		snp_sample_info_d["sample_allele"][chrom,position,sample_name] = ref_allele
	return snp_info_d, snp_sample_info_d

def snp_data_entry(snp_info_d, chrom, position, sample_name, sample_info, alt_allele, ref_allele, vcf_file, snp_sample_info_d): # Inserting sample data from vcf files into different dictionaries
	snp_info_d["multiple_allele"][chrom, position] += alt_allele
	snp_sample_info_d["sample_allele"][chrom,position,sample_name] = alt_allele[0]
	snp_info_d["files"][chrom,position].append(vcf_file)
	# snp_info_d["total_depth"][chrom,position] += int(sample_info["DP"])
	snp_info_d["alt_depth"][chrom,position] += int(sample_info["AD"])
	snp_info_d["cumulative_allele_freq"][chrom,position] += float(sample_info["FREQ"][:-1]) / 100
	snp_info_d["ref_allele"][chrom,position] = ref_allele
	if int(sample_info["DP"]) >= 10: # if reads are greater than 10 
		snp_info_d["ten_reads"][chrom,position] = True
	alt_depth_int = int(sample_info["AD"])
	ref_depth_int = int(sample_info["DP"]) - alt_depth_int
	if ref_depth_int >= 2 and alt_depth_int >= 2:
		snp_info_d["biallelic"][chrom,position] = True
		snp_sample_info_d["heterozygous"][chrom,position,sample_name] = True
	return snp_info_d, snp_sample_info_d


def cumulative_snp_dictionary(reads_list, aligner): 
	snp_info_d = {} #dictionary that contains all info on snps, 2 item tuple
	snp_sample_info_d = {} #dictionary that contains all info on snps, 3 item tuple
	snp_info_d["chromosome"] = {}
	snp_info_d["position"] = {}
	snp_info_d["files"] = {} #create dictonaries for all these variables
	snp_info_d["total_depth"] = {}
	snp_info_d["alt_depth"] = {}
	snp_info_d["ten_reads"] = {}
	snp_info_d["pass"] = {}
	snp_info_d["biallelic"] = {}
	snp_info_d["cumulative_allele_freq"] = {}
	snp_info_d["ref_allele"] = {}
	snp_sample_info_d["heterozygous"] = {} #3 item tuple
	snp_sample_info_d["sample_allele"] = {}
	snp_sample_info_d["sample_chromosome"] = {} 
	snp_sample_info_d["sample_position"] = {} 
	snp_sample_info_d["sample_name"] = {} 
	snp_info_d["multiple_allele"] = {} # to record alt allele base across all vcf files and check for snps that have multiple mutations.
	snp_info_d["fail_description"] = {} #to record why snps fail filter 
	print "all dictionaries created successfully"
	for reads in reads_list:
		sample_name, null_output = get_sample_name(reads) 
		print "importing variants for sample %s" % sample_name
		vcf_file = "%s.%s.raw.vcf" %(sample_name, aligner)
		vcf_file = os.path.join(out, vcf_file)
		vcf_reader = vcf.Reader(open(vcf_file, 'r')) # use vcf reader to read vcf file
		for index, snp in enumerate(vcf_reader): 
			chrom = snp.CHROM
			position = snp.POS
			ref_allele = snp.REF
			alt_allele = snp.ALT
			snp_info_d, snp_sample_info_d = create_empty_rows(snp_info_d, chrom, position, sample_name, ref_allele, snp_sample_info_d)
			sample_info = snp.samples[0]
			snp_info_d, snp_sample_info_d = snp_data_entry(snp_info_d, chrom, position, sample_name, sample_info, alt_allele, ref_allele, vcf_file, snp_sample_info_d)
		print "imported {} variants from sample {}".format(index+1, sample_name)
	for reads in reads_list:
		sample_name, null_output = get_sample_name(reads)
		for chrom, position in snp_info_d["pass"]:
			if not (chrom, position, sample_name) in snp_sample_info_d["sample_allele"]:
				snp_sample_info_d["sample_allele"][chrom,position,sample_name] = snp_info_d["ref_allele"][chrom,position]
	print "Removing multi-allele snps"
	snp_info_d = remove_multi_allele_snps(snp_info_d)
	return snp_info_d, snp_sample_info_d


def remove_multi_allele_snps(snp_info_d): # To remove multiples allele and non-nuclear snps 
	for chrom, position in snp_info_d["multiple_allele"].keys():
		allele_list = snp_info_d["multiple_allele"][chrom, position]
		base_d = {'A':0, 'T':0, 'G':0, 'C':0}
		for allele in allele_list:
			allele = str(allele)
			if allele in base_d:
				base_d[allele] = 1
			else:
				print "Error, {} allele found. Only A,T,C,G allowed.".format(allele)
		if sum(base_d.values()) > 1:
			snp_info_d["pass"][chrom, position] = False
			snp_info_d["fail_description"][chrom, position].append("multiple_allele")
		if chrom.startswith("NW"):
			snp_info_d["pass"][chrom, position] = False
			snp_info_d["fail_description"][chrom, position].append("Non-nuclear")
	print "Multiple allele and non-nuclear snps removed. %s SNPs failed filter, %s SNPs remain." % ((len(snp_info_d["pass"])-sum(snp_info_d["pass"].values())),sum(snp_info_d["pass"].values()))
	return snp_info_d

 
def hyperheterozygousity_filter(snp_info_d, reads_list, snp_sample_info_d, aligner):
	snp_info_d = average_minor_allele_frequency(reads_list, snp_info_d)
	hetero_prob_d = calculate_heterozygous_probability(reads_list, snp_info_d, snp_sample_info_d, aligner)
	snp_info_d, snp_sample_info_d = calculate_outcome_probability(snp_info_d, reads_list, hetero_prob_d, snp_sample_info_d)
	snp_info_d = pseudo_likelihood_score(snp_info_d, reads_list, snp_sample_info_d)
	outfile = os.path.join(out, "pseudo_likelihood_score")
	out_handle = open(outfile, "w")
	for chrom, position in snp_info_d["pseudo_likelihood_score"].keys():
		print "pseudo_likelihood_score, chromosome %s, position %d, pseudo_likelihood_score %.4f" %(chrom, position, snp_info_d["pseudo_likelihood_score"][(chrom , position)])
		outline = "%s\t%d\t%.4f\n" %(chrom, position, snp_info_d["pseudo_likelihood_score"][(chrom , position)])
		out_handle.write(outline)
		print "write into pseudo_likelihood_score text file"
	out_handle.close()
	return snp_info_d, snp_sample_info_d


def calculate_outcome_probability(snp_info_d, reads_list, hetero_prob_d, snp_sample_info_d):
	snp_sample_info_d["outcome_probability"] = {}
	for chrom, position in snp_info_d["pass"].keys():
		maf_bin = snp_info_d["maf_bin"][(chrom, position)]
		for reads in reads_list:
			sample_name, null_output = get_sample_name(reads) 
			if snp_sample_info_d["heterozygous"][(chrom, position, sample_name)] == True:
				snp_sample_info_d["outcome_probability"][(chrom, position, sample_name)] = hetero_prob_d[maf_bin]
			elif snp_sample_info_d["heterozygous"][(chrom, position, sample_name)] == False:
				snp_sample_info_d["outcome_probability"][(chrom, position, sample_name)] = 1 - hetero_prob_d[maf_bin]
			else:
				print "error %s %s %s not typable, cannot calculate outcome probabilty" %(chrom, position, sample_name)
	return snp_info_d, snp_sample_info_d


def product(list):
    p = 1
    for i in list:
        p *= i
    return p

def pseudo_likelihood_score(snp_info_d, reads_list, snp_sample_info_d):
	snp_info_d["pseudo_likelihood_score"] = {}
	for chrom , position in snp_info_d["pass"].keys():
		outcome_probability_list = []
		for reads in reads_list:
			sample_name, null_output = get_sample_name(reads) 
			outcome_probability_list.append(snp_sample_info_d["outcome_probability"][(chrom, position, sample_name)])
		outcome_probability_product = product(outcome_probability_list)
		outcome_probability_product_log = math.log10(outcome_probability_product)
		snp_info_d["pseudo_likelihood_score"][(chrom, position)] = (-1 * outcome_probability_product_log) / len(reads_list)
	return snp_info_d


def calculate_heterozygous_probability(reads_list, snp_info_d, snp_sample_info_d, aligner):
	hetero_prob_d = {}
	hetero_sum_per_bin_d = {}
	typable_sum_per_bin = {}
	total_depth_d = read_from_pickle_file("total_depth_d", aligner, "acceptable_sample_depth")
	for maf_bin in xrange(5, 51, 5):
		hetero_sum_per_bin_d[maf_bin] = 0
		typable_sum_per_bin[maf_bin] = 0
	for chrom, position in snp_info_d["pass"].keys():
		if snp_info_d["pass"][(chrom, position)] == False:
			continue
		if snp_info_d["biallelic"][(chrom, position)] == False:
			continue
		maf_bin = snp_info_d["maf_bin"][(chrom, position)]
		for reads in reads_list:
			sample_name, null_output = get_sample_name(reads) 
			if (chrom, position, sample_name) in snp_sample_info_d["heterozygous"]:
				if snp_sample_info_d["heterozygous"][(chrom, position, sample_name)] == True:
					hetero_sum_per_bin_d[maf_bin] += 1
			if total_depth_d["acceptable_sample_depth"][(chrom, position)] == True:
			# if total_depth_d[(chrom, position)] >= 5:
				typable_sum_per_bin[maf_bin] += 1
	print "hetero_sum_per_bin_d: %s" % hetero_sum_per_bin_d
	print "typable_sum_per_bin: %s" % typable_sum_per_bin
	print "hetero_sum_per_bin_d[maf_bin]: %d" % hetero_sum_per_bin_d[maf_bin]
	print "typable_sum_per_bin[maf_bin]: %d" % typable_sum_per_bin[maf_bin]
	for maf_bin in xrange(5, 51, 5):
		hetero_prob_d[maf_bin] = hetero_sum_per_bin_d[maf_bin] / typable_sum_per_bin[maf_bin]
	print "hetero_prob_d[maf_bin]: %d" % hetero_prob_d[maf_bin]
	hetero_prob_d[100] = sum(snp_sample_info_d["heterozygous"].values()) / (len(snp_info_d["pass"]) * len(reads_list)) # Check if len is the right function to return number of heterozygous snps and number of snps passing previous filters (typeable per sample)
	print "hetero_prob_d[100]: %s" % hetero_prob_d[100]
	print "Heterogenous probabilities...\n"
	print "MAF bin\theterogenous probability"
	for maf_bin in sorted(hetero_prob_d.keys()):
		print "{}\t{}".format(maf_bin,hetero_prob_d[maf_bin])
	return hetero_prob_d


def average_minor_allele_frequency(reads_list, snp_info_d): # to allocate a minor allele frequency bin for a snp position
	snp_info_d["maf_bin"] = {}
	for chrom, position in snp_info_d["pass"].keys():
		if snp_info_d["pass"][(chrom, position)] == False:
			continue
		if snp_info_d["biallelic"][(chrom, position)] == False:
			continue
		average_allele_freq = snp_info_d["cumulative_allele_freq"][(chrom, position)] / len(reads_list)
		if average_allele_freq > 0.5:
			minor_allele_freq = 1 - average_allele_freq
		elif average_allele_freq <= 0.5:
			minor_allele_freq = average_allele_freq
		else: 
			print "error %s %s average allele freq less than 0 or greater than 1" %(chrom, position)
			sys.exit(1)
		if minor_allele_freq >= 0 and minor_allele_freq <= 0.05:
			snp_info_d["maf_bin"][(chrom, position)] = 5
		elif minor_allele_freq > 0.05 and minor_allele_freq <= 0.1:
			snp_info_d["maf_bin"][(chrom, position)] = 10
		elif minor_allele_freq > 0.1 and minor_allele_freq <= 0.15:
			snp_info_d["maf_bin"][(chrom, position)] = 15
		elif minor_allele_freq > 0.15 and minor_allele_freq <= 0.2:
			snp_info_d["maf_bin"][(chrom, position)] = 20
		elif minor_allele_freq > 0.2 and minor_allele_freq <= 0.25:
			snp_info_d["maf_bin"][(chrom, position)] = 25
		elif minor_allele_freq > 0.25 and minor_allele_freq <= 0.3:
			snp_info_d["maf_bin"][(chrom, position)] = 30
		elif minor_allele_freq > 0.3 and minor_allele_freq <= 0.35:
			snp_info_d["maf_bin"][(chrom, position)] = 35
		elif minor_allele_freq > 0.35 and minor_allele_freq <= 0.4:
			snp_info_d["maf_bin"][(chrom, position)] = 40
		elif minor_allele_freq > 0.4 and minor_allele_freq <= 0.45:
			snp_info_d["maf_bin"][(chrom, position)] = 45
		elif minor_allele_freq > 0.45 and minor_allele_freq <= 0.5:
			snp_info_d["maf_bin"][(chrom, position)] = 50
		else:
			print "error %s %s is not in any bin" %(chrom, position)
			sys.exit(1)
	return snp_info_d


def rare_allele_filter(snp_info_d): # To remove positions with alleles at extremely low freq that might have resulted from alignment errors
	for chrom, position in snp_info_d["total_depth"].keys(): 
		if snp_info_d["pass"][(chrom, position)] == False:
			continue
		elif snp_info_d["alt_depth"][(chrom, position)]/snp_info_d["total_depth"][(chrom, position)] >= 0.01: 
			continue
		elif snp_info_d["ten_reads"][(chrom, position)] == True:
			continue
		else: 
			snp_info_d["pass"][(chrom, position)] = False
			snp_info_d["fail_description"][chrom, position].append("rare_allele_filter failed")
	return snp_info_d

# def coverage_filter(ref, snp_pass_d):
	# reference = SeqIO.to_dict(SeqIO.parse(ref, "gb")) # create reference list from genbank reference
	# for chrom, position in snp_pass_d.keys():
		# if snp_pass_d[(chrom, position)] == False:
			# continue
		# for feature in reference[chrom].features:
			# if feature.type == 'CDS':
				# py_snp_position = position - 1
				# if py_snp_position in feature:
					# break
		# else:
			# snp_pass_d[(chrom, position)] = False


# def gatk(temp_directory, sample_name, aligner, out): 
	# outputname = "%s.%s.doc" %(sample_name, aligner)
	# output_dedup_bam = "%s.%s.dedup.bam" %(sample_name, aligner)
	# log = os.path.join(out, "genome_log.txt")
	# with cd(temp_directory):
		# copy_from_out(output_dedup_bam, out)
		# bam_list = open("bam.list", "w")
		# bam_list.write(output_dedup_bam)
		# bam_list.write("\n")
		# bam_list.close()
		# print "%s running gatk" %sample_name
		# run_command(["gatk", "-Xmx2g", "-XX:+UseSerialGC", 
		# "-T", "DepthOfCoverage", "-R", "ref.fa", "-o", outputname, "-I", "bam.list"], log)
		# print "%s gatk completed" %sample_name
		# copy_to_results_dir(sample_name, outputname, out)


def report_progress(number, reporting_number):
	if number % reporting_number == 0:
		print "%d bases scanned" % number

# def outlier_coverage_filter(reads_list, ref, aligner, out, snp_pass_d, total_per_locus_depth):
	# reference = SeqIO.to_dict(SeqIO.parse(ref, "gb")) 
	# cds_d = {}
	# total_per_cds_locus_depth = {}
	# print "Generating empty cds dictionary..."
	# for key in total_per_locus_depth.keys():
		# cds_d[key] = False
	# print "Importing coding sequence locations..."
	# for chrom in reference.keys():
		# if chrom.startswith("NW"):
			# continue
		# for feature in reference[chrom].features:
			# if feature.type == 'CDS':
				# for locus in feature.location:
					# cds_d[(chrom, (locus+1))] = True
	# print "Generating coding sequence depth dictionary..."
	# for key in total_per_locus_depth.keys():
		# if cds_d[key] == True:
			# total_per_cds_locus_depth[key] = total_per_locus_depth[key]

	# coding_bases = len(total_per_cds_locus_depth)
	# print "number on coding bases in theileria genome = %d" % coding_bases
	# percent_85 = 0.85 *coding_bases 
	# print "85 percentile = %.2f" % percent_85
	# percent_15 = 0.15 *coding_bases
	# print "15 percentile = %.2f" % percent_15
	# percent_85_roundup = int(math.ceil(percent_85))
	# print "85 roundup percentile = %d" % percent_85_roundup
	# percent_15_rounddown = int(percent_15)
	# print "15 rounddown percentile = %d" % percent_15_rounddown
	# per_locus_depth_list = total_per_cds_locus_depth.values()
	# sortedlist = sorted(per_locus_depth_list)
	# print "per locus depth length %d" % len(per_locus_depth_list)
	# coverage_85 = sortedlist[percent_85_roundup]
	# print "coverage at 85 percentile = %d" % coverage_85
	# coverage_15 = sortedlist[percent_15_rounddown]
	# print "coverage at 15 percentile = %d" % coverage_15
	# print "running coverage filter to check if SNPs are within range"
	# for chrom, position in snp_pass_d.keys():
		# try:
			# if total_per_cds_locus_depth[(chrom, position)] >= coverage_85 or total_per_cds_locus_depth[(chrom, position)] <= coverage_15:
				# snp_pass_d[(chrom, position)] == False
		# except KeyError:
			# snp_pass_d[(chrom, position)] == False
	# return snp_pass_d

def create_empty_missingness_d(reference, reads_list):
	total_depth_d = {}
	total_depth_d["cumulative_depth"] = {}
	total_depth_d["acceptable_sample_depth"] = {}
	#Create entries for all reference bases in total_per_locus_depth dictionary
	for chrom in reference.keys():
		chrom_length = len(reference[chrom])
		for position in xrange(1, chrom_length+1):
			total_depth_d["cumulative_depth"][(chrom, position)] = 0
			total_depth_d["acceptable_sample_depth"][(chrom, position)] = True
	return total_depth_d

def sort_pileup_filter_and_report(reference, total_depth_d, reads_list, aligner, out, snp_info_d): # Missingness filter - to determine acceptable coverage depth and output to snp dictionary
	for reads in reads_list:
		sample_name, null_output = get_sample_name(reads)
		temp_depth_d = {}
		for chrom in reference.keys():
			chrom_length = len(reference[chrom])
			for position in xrange(1, chrom_length+1):
				temp_depth_d[(chrom, position)] = 0
		print "analysing sample {}".format(sample_name)
		outputname = "%s.%s.pileup" %(sample_name, aligner)
		outputname = os.path.join(out, outputname)
		depth_handle = open(outputname)
		depth_count = 0
		total_count = 0
		sum_of_depth = 0 
		print "Importing depth data from pileup..."
		for line in depth_handle:
			data1 = line.split('\t')
			sum_of_depth += int(data1[3])
			temp_depth_d[(data1[0], int(data1[1]))] += int(data1[3])
			total_depth_d["cumulative_depth"][(data1[0], int(data1[1]))] += int(data1[3])
			if (data1[0], int(data1[1])) in snp_info_d["total_depth"]:
				snp_info_d["total_depth"][(data1[0], int(data1[1]))] = total_depth_d["cumulative_depth"][(data1[0], int(data1[1]))]
		depth_handle.close()
		print "Determining acceptable coverage depth..."
		for chrom, position in temp_depth_d.keys():
			if temp_depth_d[(chrom, position)] < 5:
				total_depth_d["acceptable_sample_depth"][(chrom, position)] = False
			else:
				depth_count += 1
			total_count += 1 
		print "average coverage depth = {:.2f}x".format(float(sum_of_depth)/len(total_depth_d["cumulative_depth"]))
		for chrom, position in snp_info_d["pass"].keys():
			if total_depth_d["acceptable_sample_depth"][(chrom, position)] == False:
				snp_info_d["pass"][(chrom, position)] = False
				snp_info_d["fail_description"][chrom, position].append("Missingness filter: {} ".format(sample_name))
		print "sample %s number of bases with greater than 5 times coverage in reference genome = %d" %(sample_name,depth_count)
		print "sample %s total number of bases in reference genome = %d" %(sample_name,total_count)
		percentage_count = (float(depth_count) / total_count) * 100
		percentage_acceptable_genome = (float(sum(total_depth_d["acceptable_sample_depth"].values()))/len(total_depth_d["acceptable_sample_depth"]))*100
		print "sample %s percentage of reference genome with coverage greater than 5 times: %.2f" %(sample_name,percentage_count)
		print "percentage of reference genome acceptable for snp calling (not filtered by missingness_filter) = {:.2f}%".format(percentage_acceptable_genome)
	output_to_pickle_file(total_depth_d, "total_depth_d", "acceptable_sample_depth", aligner)
	return total_depth_d, snp_info_d

def create_empty_outlier_coverage_d(total_depth_d): #to create and empty CDS dictionary
	cds_d = {}
	total_depth_d["cds_locus_depth"] = {}
	print "Generating empty cds dictionary..."
	for key in total_depth_d["cumulative_depth"].keys():
		cds_d[key] = False
	return total_depth_d, cds_d

def bases_within_cds(reference, cds_d): #to insert snps within CDS data into CDS dictionary 
	for chrom in reference.keys():
		if chrom.startswith("NW"):
			continue
		for feature in reference[chrom].features:
			if feature.type == 'CDS':
				for locus in feature.location:
					cds_d[(chrom, (locus+1))] = True
	return cds_d

def calculate_filter_report_coverage_percentiles(total_depth_d, cds_d, snp_info_d): # Coverage filtering: to filter snps that are above the 85% percentile or below the 15% percentile and not in CDS. snps are acceptable if they are within the range and CDS. 
	for key in total_depth_d["cumulative_depth"].keys():
		if cds_d[key] == True:
			total_depth_d["cds_locus_depth"][key] = total_depth_d["cumulative_depth"][key]
	coding_bases = len(total_depth_d["cds_locus_depth"])
	print "number on coding bases in theileria genome = %d" % coding_bases
	percent_85 = 0.85 *coding_bases 
	print "85 percentile = %.2f" % percent_85
	percent_15 = 0.15 *coding_bases
	print "15 percentile = %.2f" % percent_15
	percent_85_roundup = int(math.ceil(percent_85))
	print "85 roundup percentile = %d" % percent_85_roundup
	percent_15_rounddown = int(percent_15)
	print "15 rounddown percentile = %d" % percent_15_rounddown
	per_locus_depth_list = total_depth_d["cds_locus_depth"].values()
	sortedlist = sorted(per_locus_depth_list)
	print "per locus depth length %d" % len(per_locus_depth_list)
	coverage_85 = sortedlist[percent_85_roundup]
	print "coverage at 85 percentile = %d" % coverage_85
	coverage_15 = sortedlist[percent_15_rounddown]
	print "coverage at 15 percentile = %d" % coverage_15
	print "running coverage filter to check if SNPs are within range"
	coverage_85_count = 0
	coverage_15_count = 0
	coverage_non_cds = 0
	for chrom, position in snp_info_d["pass"].keys():
		if cds_d[chrom, position] == True:
			if total_depth_d["cds_locus_depth"][(chrom, position)] >= coverage_85:  
				snp_info_d["pass"][(chrom, position)] = False
				coverage_85_count += 1
				snp_info_d["fail_description"][chrom, position].append("outlier coverage filter: snp is greater than 85 percentile")
			if total_depth_d["cds_locus_depth"][(chrom, position)] <= coverage_15:
				snp_info_d["pass"][(chrom, position)] = False
				coverage_15_count += 1
				snp_info_d["fail_description"][chrom, position].append("outlier coverage filter: snp is lower than 85 percentile")
		else:
			snp_info_d["pass"][(chrom, position)] = False
			coverage_non_cds += 1
			snp_info_d["fail_description"][chrom, position].append("outlier coverage filter: snp in non CDS")
	print "number of snps greater than 85 percentile = {}, number of snps lower than 15 percentile = {}, number of snps not in CDS = {}".format(coverage_85_count, coverage_15_count, coverage_non_cds)
	return snp_info_d, total_depth_d

def missingness_filter(snp_info_d, reads_list, aligner, out): # missingness_filter to filter snps that have certain criteria missing. A compilation of different functions. Also includes coverage filtering
	reference = SeqIO.to_dict(SeqIO.parse(ref, "gb"))
	total_depth_d = create_empty_missingness_d(reference, reads_list)
	total_depth_d, snp_info_d = sort_pileup_filter_and_report(reference, total_depth_d, reads_list, aligner, out, snp_info_d)
	print "Missingness filter complete: %s SNPs failed filter, %s SNPs remain." % ((len(snp_info_d["pass"])-sum(snp_info_d["pass"].values())),sum(snp_info_d["pass"].values()))
	print "running outlier coverage filter"
#outlier coverage filter
	total_depth_d, cds_d = create_empty_outlier_coverage_d(total_depth_d)
	print "Importing coding sequence locations..."
	cds_d = bases_within_cds(reference, cds_d)
	print "Generating coding sequence depth dictionary..."
	snp_info_d, total_depth_d = calculate_filter_report_coverage_percentiles(total_depth_d, cds_d, snp_info_d)
	return snp_info_d, total_depth_d


def find_iter(count, slice_sequence, search_sequence):
	for match in re.finditer(slice_sequence, search_sequence):
		count += 1
		if count > 1:
			return count
	return count



def unique_window_match(reference, str_slice, str_rev_slice): 
	count = 0
	for search_chrom in reference.keys():
		search_sequence = str(reference[search_chrom].seq)
		count = find_iter(count, str_slice, search_sequence)
		if count > 1:
			return True
		count = find_iter(count, str_rev_slice, search_sequence)
		if count > 1:
			return True
	return False


def calculate_uniqueness_score(index, reference, chrom, window_length, window_step, attempted_window_list): # To calculate uniqueness score of each individual snp.
	while True:
		slice_start = index
		slice_end = slice_start + window_length
		iteration = 1
		if window_length in attempted_window_list:
			return attempted_window_list
		attempted_window_list.append(window_length)
		while slice_end >= index and slice_start >= 0: # while loop for uniquess filter 
			# print slice_start
			# print slice_end
			if slice_end > len(reference[chrom]):
				slice_end = len(reference[chrom])
				slice_start = slice_end - window_length 
			slice = reference[chrom].seq[slice_start:slice_end] #defining the slice
			str_slice = str(slice)
			rev_slice = slice.reverse_complement() #defining the reverse search motif
			str_rev_slice = str(rev_slice)
				# print "Found %d matches in chromosome %s" % (matches, search_chrom)
				# print "Found %d reverse complement matches in chromosome %s" % (rev_matches, search_chrom)
			if unique_window_match(reference, re.compile(str_slice), re.compile(str_rev_slice)) == True: # if count is greater than 1, increase window size
				# print "uniqueness_score not identified for window_length %s, chromosome %s, position %d. Iteration %d has %d matches." %(window_length, reference[chrom].id, index, iteration, count)
				break
				#print "uniqueness_score not identified for window_length %s, chromosome %s, position %d" %(window_length, reference[chrom].id, index)
			slice_start -= 1 # move slice start position down by 1
			slice_end -= 1 # move slice end position down by 1
			iteration += 1
		else:
			return attempted_window_list
		window_length += window_step # reducing window length by 1



def generate_reference_uniqueness_file(argument): 
	reference, chrom, base_no = argument
	attempted_window_list = [0]
	attempted_window_list = calculate_uniqueness_score(base_no, reference, chrom, 16, 16, attempted_window_list)
	# print "First iteration complete, chromosome %s, position %d of %d total bases, uniqueness_score less than %d" %(reference[chrom].id, index, len(reference[chrom]), max(attempted_window_list))
	attempted_window_list.remove(max(attempted_window_list))
	# print attempted_window_list
	attempted_window_list = calculate_uniqueness_score(base_no, reference, chrom, (max(attempted_window_list) + 4), 4, attempted_window_list)
	# print "Second iteration complete, chromosome %s, position %d of %d total bases, uniqueness_score less than %d" %(reference[chrom].id, index, len(reference[chrom]), max(attempted_window_list))
	attempted_window_list.remove(max(attempted_window_list))
	# print attempted_window_list
	attempted_window_list = calculate_uniqueness_score(base_no, reference, chrom, (max(attempted_window_list) + 1), 1, attempted_window_list)
	# print "Third iteration complete, chromosome %s, position %d of %d total bases, uniqueness_score less than %d" %(reference[chrom].id, index, len(reference[chrom]), max(attempted_window_list))
	# print attempted_window_list
	uniqueness_score = max(attempted_window_list)
	print "Uniqueness_score found, chromosome %s, position %d of %d total bases, uniqueness_score %d" %(reference[chrom].id, base_no, len(reference[chrom]), uniqueness_score)
	return (chrom, base_no, uniqueness_score)



def uniqueness_filter(ref, out, snp_info_d, aligner, threads):
	reference = SeqIO.to_dict(SeqIO.parse(ref, "gb")) # create reference list from genbank reference
	outfile = "%s.uniqueness_score" %(aligner)
	outfile = os.path.join(out, outfile)
	if not os.path.isfile(outfile):
		out_handle = open(outfile, "w")
		argument_list = []
		for chrom, base_no in snp_info_d["pass"].keys(): 
			if snp_info_d["pass"][(chrom, base_no)] == False:
				continue
			argument_list.append((reference, chrom, base_no))
		results_list = []
		print "Calculating uniqueness scores..."
		# for argument in argument_list:
			# results_list.append(generate_reference_uniqueness_file(argument))
		p = multiprocessing.Pool(threads)
		result = p.map_async(generate_reference_uniqueness_file, argument_list, chunksize=1)
		prev_jobs_remaining = len(argument_list)
		while not result.ready():
			jobs_remaining = result._number_left
			if jobs_remaining != prev_jobs_remaining:
				print "%d of %d bases remaining to be processed" % (result._number_left, len(argument_list))
			prev_jobs_remaining = jobs_remaining
		results_list = result.get(999999999999999)
		p.close()
		p.join()
		print "Writing results to file..."
		output_dict = defaultdict(dict)
		for chrom, base_no, uniqueness_score in results_list:
			output_dict[chrom][base_no] = uniqueness_score
		for chrom in sorted(output_dict.keys()):
			sub = output_dict[chrom]
			for index in sorted(sub.keys()):
				if output_dict[chrom][index] > 100:
					snp_info_d["pass"][(chrom, index)] = False
					snp_info_d["fail_description"][chrom, index].append("uniqueness filter: score greater than 100")
				outline = "%s\t%d\t%d\n" % (chrom, index, output_dict[chrom][index])
				out_handle.write(outline)
		out_handle.close()
	else:
		out_handle = open(outfile)
		output_dict = defaultdict(dict)
		for line in out_handle:
			data = line.strip().split('\t')
			chrom = data[0]
			base_no = int(data[1])
			uniqueness_score = int(data[2])
			output_dict[chrom][base_no] = uniqueness_score
		for chrom, base_no in snp_info_d["pass"].keys():
			if snp_info_d["pass"][(chrom, base_no)] == True:
				if output_dict[chrom][base_no] > 100:
					snp_info_d["pass"][(chrom, base_no)] = False
					snp_info_d["fail_description"][chrom, base_no].append("uniqueness filter: score greater than 100")
	return snp_info_d

def import_reads(reads_list):
	with open(reads_list, "r") as reads_handle:
	# reads_handle=open(reads_list, "r")
	# reads_handle.close()
		reads_iterable=[]
		for line in reads_handle:
			stripline=line.strip()
			if stripline != '':
				reads_iterable.append(stripline)
	return reads_iterable

def replace_dir(dir):
	if os.path.exists(dir):
		shutil.rmtree(dir)
	os.mkdir(dir)

def make_output_directory(out, force):
	if force=="true":
		replace_dir(out)
	elif force=="false":
		if os.path.exists(out) == False:
			os.mkdir(out)
	else:
		print "error, can't identify force variable" 
		sys.exit(1)

def copy_to_results_dir(sample_name, source, destination_dir):
	hold_dir = os.path.join(destination_dir,".copy", sample_name)
	if os.path.isdir(os.path.join(destination_dir,".copy")) == False:
		os.mkdir(os.path.join(destination_dir,".copy"))
	if os.path.isdir(hold_dir) == False:
		os.mkdir(hold_dir)
	if os.path.isdir(source):
		shutil.copytree(source,os.path.join(hold_dir,os.path.basename(source)))
		os.rename(os.path.join(hold_dir,os.path.basename(source)), os.path.join(destination_dir,os.path.basename(source)))
		shutil.rmtree(hold_dir)
	elif os.path.isfile(source):
		shutil.copyfile(source,os.path.join(hold_dir,os.path.basename(source)))
		os.rename(os.path.join(hold_dir,os.path.basename(source)), os.path.join(destination_dir,os.path.basename(source)))
		shutil.rmtree(hold_dir)
	else:
		print "Error copying object: %s cannot determine if directory or file." % source
		shutil.rmtree(hold_dir)
		sys.exit(1)

def generate_alignments(reads):
	reads = reads.rstrip()
	sample_name, gzipped = get_sample_name(reads)
	trim_reads(reads, sample_name, temp_directory, out, gzipped)
	# biobloomcategorizer(reads, sample_name, temp_directory, out)
	bowtie2_mapping(reads, sample_name, temp_directory, out)
	bwa(sample_name, temp_directory, out)
	ngm(sample_name, temp_directory, out)
	sort_sam(sample_name, temp_directory, out)
	mark_duplicates(sample_name, temp_directory, out)
	samtools_mpileup(sample_name, temp_directory, out)
	varscan(sample_name, temp_directory, out)
	# for aligner in ["bwa", "ngm"]:
		# gatk(temp_directory, sample_name, aligner, out)


def save_to_csv(dictionary, name, filter, aligner):
	out_dict = defaultdict(dict)
	base_path = os.path.dirname(__file__)
	temp_directory = os.path.join(base_path,"temp_directory")
	file_name = "{}.{}.{}.tsv".format(name, aligner, filter)
	path = os.path.join(temp_directory,file_name)
	for column in dictionary.keys():
		for row in dictionary[column].keys():
			out_dict[row][column] = dictionary[column][row]
	with open(path, 'w') as fp:
		column_names = list(sorted(dictionary.keys()))
		writer = csv.DictWriter(fp, fieldnames = column_names, delimiter='\t')
		writer.writeheader()
		for row in sorted(out_dict.keys()):
			writer.writerow(out_dict[row])


def output_to_pickle_file(dictionary, name, filter, aligner): # write python dict to a file
	# dictionary = defaultdict(dict)
	base_path = os.path.dirname(__file__)
	file_name = "{}.{}.{}.pkl".format(name, aligner, filter)
	temp_directory = os.path.join(base_path,"temp_directory")
	path = os.path.join(out,file_name)
	with open(path, 'wb') as f:
		pickle.dump(dictionary, f, protocol=pickle.HIGHEST_PROTOCOL)

# write python dict to a file
# mydict = {'a': 1, 'b': 2, 'c': 3}
# output = open('myfile.pkl', 'wb')
# pickle.dump(mydict, output)
# output.close()

def read_from_pickle_file(name, aligner, filter):  # read python dict back from the file
	# dictionary = defaultdict(dict)
	base_path = os.path.dirname(__file__)
	file_name = "{}.{}.{}.pkl".format(name, aligner, filter)
	temp_directory = os.path.join(base_path,"temp_directory")
	path = os.path.join(out,file_name)
	with open(path, 'rb') as f:
		return pickle.load(f)

# read python dict back from the file
# pkl_file = open('myfile.pkl', 'rb')
# mydict2 = pickle.load(pkl_file)
# pkl_file.close()

def str2bool(v):
	return v.lower() in ("yes", "true", "t", "1")


def assign_correct_data_type_on_import(dictionary, column, row):
	if column == "alt_depth":
		dictionary[column][(row["chromosome"],int(row["position"]))] = int(row[column]) #This column is integers, changing to integer entry
	if column == "biallelic":
		dictionary[column][(row["chromosome"],int(row["position"]))] = str2bool(row[column]) #boolean entry
	if column == "chromosome":
		dictionary[column][(row["chromosome"],int(row["position"]))] = row[column] #string entry
	if column == "cumulative_allele_freq":
		dictionary[column][(row["chromosome"],int(row["position"]))] = float(row[column]) #floating point number entry
	if column == "files":
		dictionary[column][(row["chromosome"],int(row["position"]))] = row[column][1:-1].split(', ') #list entry
	if column == "multiple_allele":
		dictionary[column][(row["chromosome"],int(row["position"]))] = row[column][1:-1].split(', ')
	if column == "pass":
		dictionary[column][(row["chromosome"],int(row["position"]))] = str2bool(row[column])
	if column == "position":
		dictionary[column][(row["chromosome"],int(row["position"]))] = int(row[column])
	if column == "ten_reads":
		dictionary[column][(row["chromosome"],int(row["position"]))] = str2bool(row[column])
	if column == "total_depth":
		dictionary[column][(row["chromosome"],int(row["position"]))] = int(row[column])
	if column == "maf_bin":
		dictionary[column][(row["chromosome"],int(row["position"]))] = int(row[column])
	if column == "pseudo_likelihood_score":
		dictionary[column][(row["chromosome"],int(row["position"]))] = float(row[column])
	if column == "heterozygous":
		dictionary[column][(row["chromosome"],row["position"], row["sample_name"])] = str2bool(row[column])
	if column == "sample_chromosome":
		dictionary[column][(row["chromosome"],row["position"], row["sample_name"])] = row[column]
	if column == "sample_position":
		dictionary[column][(row["chromosome"],row["position"], row["sample_name"])] = int(row[column])
	if column == "sample_name":
		dictionary[column][(row["chromosome"],row["position"], row["sample_name"])] = row[column]
	if column == "outcome_probability":
		dictionary[column][(row["chromosome"],row["position"], row["sample_name"])] = float(row[column])
	if column == "fail_description":
		dictionary[column][(row["chromosome"],int(row["position"]))] = row[column][1:-1].split(', ') #list entry
	return dictionary


def open_csv(name, aligner, filter):
	dictionary = defaultdict(dict)
	base_path = os.path.dirname(__file__)
	file_name = "{}.{}.{}.tsv".format(name, aligner, filter)
	temp_directory = os.path.join(base_path,"temp_directory")
	path = os.path.join(out,file_name)
	with open(path, 'r') as fp:
		reader = csv.DictReader(fp, delimiter='\t')
		for row in reader:
			for column in row.keys():
				dictionary = assign_correct_data_type_on_import(dictionary, column, row)
	return dictionary


def filter_alignments(reads_list, aligners, threads):
	base_path = os.path.dirname(__file__)
	temp_directory = os.path.join(base_path,"temp_directory")
	for aligner in aligners:
		filename = "{}.{}.{}.tsv".format("snp_info_d", aligner, "cumulative_snp_dictionary")
		if not os.path.isfile(os.path.join(out, filename)):
			print "filtering reads generated by %s" % aligner
			print "generating cumulative SNP dictionary."
			snp_info_d, snp_sample_info_d = cumulative_snp_dictionary(reads_list, aligner)
			save_to_csv(snp_info_d, "snp_info_d", "cumulative_snp_dictionary", aligner)
			output_to_pickle_file(snp_info_d, "snp_info_d", "cumulative_snp_dictionary", aligner)
			save_to_csv(snp_sample_info_d, "snp_sample_info_d", "cumulative_snp_dictionary", aligner)
			output_to_pickle_file(snp_sample_info_d, "snp_sample_info_d", "cumulative_snp_dictionary", aligner)
			os.rename(os.path.join(temp_directory,filename),os.path.join(out,filename))
			filename = "{}.{}.{}.tsv".format("snp_sample_info_d", aligner, "cumulative_snp_dictionary")
			os.rename(os.path.join(temp_directory,filename),os.path.join(out, filename))
		filename = "{}.{}.{}.tsv".format("snp_info_d", aligner, "missingness_filter")
		if not os.path.isfile(os.path.join(out, filename)):
			# snp_info_d = open_csv("snp_info_d", aligner, "cumulative_snp_dictionary")
			snp_info_d = read_from_pickle_file("snp_info_d", aligner, "cumulative_snp_dictionary")
			print "Running missingness filter"
			snp_info_d, total_depth_d = missingness_filter(snp_info_d, reads_list, aligner, out)
			print "outlier coverage filter complete: %s SNPs failed filter, %s SNPs remain." % ((len(snp_info_d["pass"])-sum(snp_info_d["pass"].values())),sum(snp_info_d["pass"].values()))
			save_to_csv(snp_info_d, "snp_info_d", "missingness_filter", aligner)
			output_to_pickle_file(snp_info_d, "snp_info_d", "missingness_filter", aligner)
			os.rename(os.path.join(temp_directory,filename),os.path.join(out,filename))
		filename = "{}.{}.{}.tsv".format("snp_info_d", aligner, "rare_allele_filter")
		if not os.path.isfile(os.path.join(out, filename)):
			# snp_info_d = open_csv("snp_info_d", aligner, "missingness_filter")
			snp_info_d = read_from_pickle_file("snp_info_d", aligner, "missingness_filter")
			print "Running rare allele filter"
			snp_info_d = rare_allele_filter(snp_info_d)
			print "Rare allele filter complete: %s SNPs failed filter, %s SNPs remian " % ((len(snp_info_d["pass"])-sum(snp_info_d["pass"].values())),sum(snp_info_d["pass"].values()))
			save_to_csv(snp_info_d, "snp_info_d", "rare_allele_filter", aligner)
			output_to_pickle_file(snp_info_d, "snp_info_d", "rare_allele_filter", aligner)
			os.rename(os.path.join(temp_directory,filename),os.path.join(out,filename))
		filename = "{}.{}.{}.tsv".format("snp_info_d", aligner, "uniqueness_filter")
		if not os.path.isfile(os.path.join(out, filename)):
			# snp_info_d = open_csv("snp_info_d", aligner, "rare_allele_filter")
			snp_info_d = read_from_pickle_file("snp_info_d", aligner, "rare_allele_filter")
			print "Running uniqueness filter for %s." %(aligner)
			snp_info_d = uniqueness_filter(ref, out, snp_info_d, aligner, threads)
			print "%s uniqueness filter complete: %s SNPs failed filter, %s SNPs remain." % (aligner,(len(snp_info_d["pass"])-sum(snp_info_d["pass"].values())),sum(snp_info_d["pass"].values()))
			save_to_csv(snp_info_d, "snp_info_d", "uniqueness_filter", aligner)
			output_to_pickle_file(snp_info_d, "snp_info_d", "uniqueness_filter", aligner)
			os.rename(os.path.join(temp_directory,filename),os.path.join(out,filename))
		# filename = "{}.{}.{}.tsv".format("snp_info_d", aligner, "hyperheterozygousity_filter")
		# if not os.path.isfile(os.path.join(out, filename)):
			snp_info_d = open_csv("snp_info_d", aligner, "uniqueness_filter")
			# snp_info_d = read_from_pickle_file("snp_info_d", aligner, "uniqueness_filter")
			# snp_sample_info_d = read_from_pickle_file("snp_sample_info_d", aligner, "cumulative_snp_dictionary")
			# print "Running hyperheterozygousity filter"
			# snp_info_d, snp_sample_info_d = hyperheterozygousity_filter(snp_info_d, reads_list, snp_sample_info_d, aligner)
			# print "hyperheterozygousity filter complete: %s SNPs failed filter, %s SNPs remain." % ((len(snp_info_d["pass"])-sum(snp_info_d["pass"].values())),sum(snp_info_d["pass"].values()))
			# save_to_csv(snp_info_d, "snp_info_d", "hyperheterozygousity_filter", aligner)
			# save_to_csv(snp_sample_info_d, "snp_sample_info_d", "hyperheterozygousity_filter", aligner)
			# output_to_pickle_file(snp_info_d, "snp_info_d", "hyperheterozygousity_filter", aligner)
			# output_to_pickle_file(snp_sample_info_d, "snp_sample_info_d", "hyperheterozygousity_filter", aligner)
			# os.rename(os.path.join(temp_directory,filename),os.path.join(out,filename))
			# filename = "{}.{}.{}.tsv".format("snp_sample_info_d", aligner, "hyperheterozygousity_filter")
			# os.rename(os.path.join(temp_directory,filename),os.path.join(out,filename))

def run_pipeline_local(threads, reads_list):
	# for read_data in reads_list:
		# generate_alignments(read_data)
	p = multiprocessing.Pool(int(threads))
	p.map_async(generate_alignments,reads_list).get(9999999)
	p.close()
	p.join()
	filter_alignments(reads_list, ["bwa", "ngm"], threads)


reads_list = import_reads(os.path.abspath(sys.argv[1]))
ref = os.path.abspath(sys.argv[2])
out = os.path.abspath(sys.argv[3])
force = sys.argv[4]
threads = int(sys.argv[5])
host_ref = os.path.abspath(sys.argv[6])
base_path = os.path.dirname(__file__)

make_output_directory(out, force)
temp_directory = os.path.join(base_path,"temp_directory")
replace_dir(temp_directory)
index_reference(ref, host_ref, temp_directory)

run_pipeline_local(threads, reads_list)
