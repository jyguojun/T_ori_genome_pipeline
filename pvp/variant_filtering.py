#!/usr/bin/env python

import os, sys, gzip, subprocess, shutil, csv, re, math, multiprocessing, vcf
from collections import defaultdict
from pvp_utils import *


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
	if not uniqueness:
		uniqueness = os.path.join(out, "%s.uniqueness_score" %(aligner))
	if not os.path.isfile(uniqueness):
		out_handle = open(uniqueness, "w")
		argument_list = []
		for chrom, base_no in snp_info_d["pass"].keys(): 
			if snp_info_d["pass"][(chrom, base_no)] == False:
				continue
			argument_list.append((reference, chrom, base_no))
		results_list = []
		print "Calculating uniqueness scores..."
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
		out_handle = open(uniqueness)
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


def filter_alignments(reads, aligners, threads):
	for aligner in aligners:
		filename = "{}.{}.{}.tsv".format("snp_info_d", aligner, "cumulative_snp_dictionary")
		if not os.path.isfile(os.path.join(out, filename)):
			print "filtering reads generated by %s" % aligner
			print "generating cumulative SNP dictionary."
			snp_info_d, snp_sample_info_d = cumulative_snp_dictionary(reads, aligner)
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
			snp_info_d, total_depth_d = missingness_filter(snp_info_d, reads, aligner, out)
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
