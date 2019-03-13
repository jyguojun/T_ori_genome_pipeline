#!/usr/bin/env python

import os, sys, gzip, subprocess, shutil, csv, re, math, multiprocessing, vcf
from contextlib import contextmanager
import cPickle as pickle
from collections import defaultdict

@contextmanager
def cd(newdir):
	prevdir = os.getcwd()
	os.chdir(os.path.expanduser(newdir))
	try:
		yield
	finally:
		os.chdir(prevdir)

def create_log_entry(log_file, log_entry):
	with open(log_file, 'a+') as log_handle:
		log_handle.write(log_entry + '\n')

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


def copy_from_out(filename, out):
	if os.path.exists(filename) == False:
		shutil.copyfile(os.path.join(out, filename),filename)


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
		print "deleting {}".format(dir)
		shutil.rmtree(dir)
	os.mkdir(dir)


def make_output_directory(out, force):
	if force==True:
		replace_dir(out)
	elif force==False:
		if os.path.exists(out) == False:
			os.mkdir(out)
	else:
		print "error, can't identify force variable" 
		sys.exit(1)


def save_to_csv(temp_directory, dictionary, name, filter, aligner):
	out_dict = defaultdict(dict)
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

			
def output_to_pickle_file(out, dictionary, name, filter, aligner): # write python dict to a file
	# dictionary = defaultdict(dict)
	file_name = "{}.{}.{}.pkl".format(name, aligner, filter)
	path = os.path.join(out,file_name)
	with open(path, 'wb') as f:
		pickle.dump(dictionary, f, protocol=pickle.HIGHEST_PROTOCOL)


def read_from_pickle_file(out, name, aligner, filter):  # read python dict back from the file
	# dictionary = defaultdict(dict)
	file_name = "{}.{}.{}.pkl".format(name, aligner, filter)
	path = os.path.join(out,file_name)
	with open(path, 'rb') as f:
		return pickle.load(f)


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
	file_name = "{}.{}.{}.tsv".format(name, aligner, filter)
	path = os.path.join(out,file_name)
	with open(path, 'r') as fp:
		reader = csv.DictReader(fp, delimiter='\t')
		for row in reader:
			for column in row.keys():
				dictionary = assign_correct_data_type_on_import(dictionary, column, row)
	return dictionary

