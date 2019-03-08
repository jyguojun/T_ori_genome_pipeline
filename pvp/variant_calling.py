#!/usr/bin/env python

import os
from pvp_utils import *

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
