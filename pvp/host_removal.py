#!/usr/bin/env python

import os, shutil
from pvp_utils import *
from Bio import SeqIO

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

def index_reference(out, ref, host_ref, temp_directory, threads): #indexing to create samtools, bwa, picard and bowtie2 (cattle host genome) index
	ref_gb = SeqIO.parse(ref, "gb")
	SeqIO.write(ref_gb, os.path.join(temp_directory, "ref.fa"), "fasta")
	# host_ref_gb = SeqIO.parse(host_ref, "gb") # removed as host genome is already in fasta format
	# SeqIO.write(host_ref_gb, os.path.join(temp_directory, "host_ref.fa"), "fasta")
	log = os.path.join(out,"genome_log.txt")
	ref_dir = os.path.join(out, "references")
	if not os.path.isdir(ref_dir):
		os.mkdir(ref_dir)
	if not os.path.exists(os.path.join(ref_dir, 'host_ref_db.1.bt2')):
		with cd(temp_directory): 
			run_command(["bwa", "index", "ref.fa" ], log)
			print "Created bwa index"
			run_command(["samtools", "faidx", "ref.fa"], log)
			print "Created samtools index"
			run_command(["picard", "-Xmx2g", "-XX:+UseSerialGC", "CreateSequenceDictionary", "REFERENCE=ref.fa", "OUTPUT=ref.dict" ], log)
			print "Created picard index"
			run_command(["bowtie2-build", host_ref, "host_ref_db", "--threads", str(threads)], log)
			print "Created bowtie2 cattle host genome index"
			for file in os.listdir(temp_directory):
				shutil.copyfile(os.path.join(ref_dir, file), os.path.join(temp_directory, file))
	else:
		for file in os.listdir(ref_dir):
			shutil.copyfile(os.path.join(ref_dir, file), os.path.join(temp_directory, file))

