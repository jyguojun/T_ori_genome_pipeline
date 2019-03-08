#!/usr/bin/env python

import os
from pvp_utils import *

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
