"""
File:		full_pipeline.py
Author:		Somayeh Kafaie
Date:		September 2020


Purpose:	To extract the neighborhood of an AMR gene from an assembly graph,
			annotate it and compare with the reference genome(s)

To run:
	python full_pipeline.py

NOTE: all parameters can be set in params.py
NOTE: if use_RGI = TRUE, make sure either RGI has been installed system-wide or
	you already are in the environment RGI installed in!
"""

import sys
import os
import errno
import copy
import gfapy
import re
import argparse
import difflib
import datetime
import csv
from csv import DictReader
import collections
from Bio import SeqIO
from gfapy.sequence import rc
import enum
import subprocess
import random
import shutil
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import json
import logging
from functools import partial
from multiprocessing.pool import Pool
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import math
import yaml

from sarand.params import Pipeline_tasks, Assembler_name
from sarand.extract_neighborhood import neighborhood_sequence_extraction
from sarand.annotation_visualization import visualize_annotation
from sarand.utils import initialize_logger, str2bool, verify_file_existence,\
			retrieve_AMR, extract_files, create_fasta_file,\
			annotate_sequence, split_up_down_info,\
			seqs_annotation_are_identical, similar_seq_annotation_already_exist,\
			amr_name_from_comment,\
			amr_name_from_title, retreive_original_amr_name, extract_name_from_file_name,\
			restricted_amr_name_from_modified_name, extract_info_from_overlap_file,\
			read_path_info_from_align_file, read_path_info_from_align_file_with_multiple_amrs,\
			extract_path_info_for_amrs, compare_two_sequences,\
			delete_lines_started_with, validate_task_values, validate_print_parameters_tools

ASSEMBLY_FILE = 'assembly_graph_with_scaffolds.gfa'
#ALL_AMR_SEQUENCES ='nucleotide_fasta_protein_homolog_model_without_efflux_without_space.fasta'
AMR_DIR_NAME = 'AMR_info'
AMR_SEQ_DIR = 'sequences'
AMR_ALIGN_DIR = 'alignments'
AMR_OVERLAP_FILE = 'overlaps.txt'
SUBGRAPH_DIR_NAME = 'subgraphs'
SEQ_DIR_NAME = 'sequences_info'
SEQ_NAME_PREFIX = 'ng_sequences_'
ANNOTATION_DIR = 'annotations'
EVAL_DIR = 'evaluation'
NOT_FOUND_FILE = 'not_found_amrs_in_graph.txt'

#To run the code for a list of sequence neighborhood length rather than just one length
#the default seq length is 1000
MULTIPLE_SEQ_LENGTH = False
seq_eval_ref_subject = False

def write_info_in_annotation_file(annotation_writer, visual_annotation_writer,
								gene_info, use_RGI, found, len_seq = None):
	"""
	To write annotation details into files
	Parameters:
		annotation_writer:	annotation file containing all annotations
		visual_annotation_writer: annotation file containing unique annotations
		seq_description: a small description of the sequence used for naming
		seq: the extracted dna sequence that has been annotated
		gene_info: annotation info
		contig_name: the name of contig matching this extracted sequence (if there is any contig)
		use_RGI: if True, RGI has been used to annotate AMRs
		found: if True, the annotation info has already found in other annotated sequences
	"""
	seq = gene_info['seq_value']
	seq_description = gene_info['seq_name']
	if len_seq is None:
		len_seq = len(seq)
	if use_RGI:
		annotation_writer.writerow([seq_description, seq, len_seq,
							gene_info['gene'], gene_info['prokka_gene_name'],
							gene_info['product'], gene_info['length'],
							gene_info['start_pos'], gene_info['end_pos'],
							gene_info['RGI_prediction_type'], gene_info['coverage'],
							gene_info['family'], gene_info['target_amr']])
		if not found:
			visual_annotation_writer.writerow([seq_description, seq, len_seq,
								gene_info['gene'], gene_info['prokka_gene_name'],
								gene_info['product'], gene_info['length'],
								gene_info['start_pos'], gene_info['end_pos'],
								gene_info['RGI_prediction_type'], gene_info['coverage'],
								gene_info['family'], gene_info['target_amr']])
	else:
		annotation_writer.writerow([seq_description, seq, len_seq,
							gene_info['gene'],
							gene_info['product'], gene_info['length'],
							gene_info['start_pos'], gene_info['end_pos'],
							gene_info['coverage'], gene_info['target_amr']])
		if not found:
			visual_annotation_writer.writerow([seq_description, seq, len_seq,
								gene_info['gene'],
								gene_info['product'], gene_info['length'],
								gene_info['start_pos'], gene_info['end_pos'],
								gene_info['coverage'], gene_info['target_amr']])

def seq_annotation_already_exist(seq_info_list, all_seq_info_lists, out_dir):
	"""
	To check if annotations found for the new sequence have already exists in the
	list of annotations extracted from other sequences.
	Parameters:
		seq_info_list: list of annotations of the new sequence
		all_seq_info_lists: list of annotations of the other sequences annotated so far
	Return:
		True if seq_info_list is identical to another list in all_seq_info_lists in
		terms of 'gene', 'length','start_pos' and 'end_pos' for all thei items
	"""
	found = False
	for seq_list in all_seq_info_lists:
		if seqs_annotation_are_identical(seq_info_list, seq_list, out_dir, 100):
			found = True
			break

	return found

def read_path_info_file(path_info_file):
	"""
	To read the csv file containing the list of nodes and their coverage representing
	in an extracted sequence
	Parameter:
		path_info_file: the csv file
	Return:
		stored info in a list
	"""
	seq_info_list = []
	seq_info = []
	with open(path_info_file, 'r') as myfile:
		myreader = DictReader(myfile)
		old_seq = ''
		for row in myreader:
			node_info = {'node':row['node'], 'coverage':float(row['coverage']),
			'start':int(row['start']), 'end':int(row['end'])}
			cur_seq = row['sequence']
			if cur_seq!=old_seq:
				if (seq_info):
					seq_info_list.append(seq_info)
				seq_info = []
				old_seq = cur_seq
			seq_info.append(node_info)
		seq_info_list.append(seq_info)
	return seq_info_list

def find_gene_coverage(seq_info_list, path_info):
	"""
	Calculate the coverage of genes available in a sequence based on the coverage
	of the nodes representing them
	Parameter:
		seq_info_list: the list of genes (and their info) annotated in a sequence
		path_info:	the list of nodes (and their info) representing the sequence
	Return:
		the list of calculated gene-coverages
	"""
	coverage_list = []
	for seq_info in seq_info_list:
		sum_coverage = 0
		#minus 1 because I guess in what prokka returns the sequence starts from position 1
		start, end = seq_info['start_pos']-1, seq_info['end_pos']-1
		if start > end:
			start, end = end, start
		found = False
		for j, node_info in enumerate(path_info):
			if node_info['start']<= start and node_info['end']>= start:
				found = True
				#the annotated gene is represented by a single node
				if node_info['end'] >= end:
					sum_coverage = (end - start + 1) * node_info['coverage']
				else:
					#calculate the length of node representing the gene
					sum_coverage = (node_info['end'] - start + 1)*node_info['coverage']
					n_index = j + 1
					assert n_index < len(path_info), "wrong index calculated for path_info!!"
					while (n_index < len(path_info) and path_info[n_index]['start']<=end):
						if path_info[n_index]['end'] >= end:
							sum_coverage+=(end-path_info[n_index]['start']+1)*path_info[n_index]['coverage']
							break
						else:
							sum_coverage+=(path_info[n_index]['end']-path_info[n_index]['start']+1)*path_info[n_index]['coverage']
						n_index+=1
				coverage_list.append(sum_coverage / (end-start+1))
				break
		if not found:
			logging.error("ERROR: no nodes were found for this gene!!!")
			import pdb; pdb.set_trace()
			sys.exit()
	return coverage_list

def find_target_amr_in_seqvalue_and_return_coverage(seq_info):
	"""
	To find the target AMR in a sequence and its covearge:
	We first look at the extracted sequence as the part of the sequence in lower-case
	represents the AMR; then, we look at the list of annotated genes and find the
	gene that has the most overlap with the lower-case area and finally find its coverage
	Parameters:
		seq_info: the details of an extracted sequences and its annotation
	Return:
		the coverage of the AMR, its index in seq_info and False if we can find a
		gene in the lower-case area of the sequence; otherwise it returns True
	"""
	error = True
	amr_coverage = 0
	#find the indeces of lower case string (amr sequence) in the extracted sequence
	sequence = seq_info[0]['seq_value']
	amr_start = -1
	amr_end = -1
	index = 0
	while index < len(sequence):
		if sequence[index].islower() and amr_start==-1:
			amr_start = index
		elif sequence[index].isupper() and amr_start>-1:
			amr_end=index-1
			break
		index+=1
	#find the gene has the most overlap with the found range
	#we look at all found genes that their overlap with the sequence is more than initial value of most_overlap and chose the one with max
	most_overlap = 50
	amr_index = -1
	for i, gene_info in enumerate(seq_info):
		start, end = min(gene_info['start_pos'], gene_info['end_pos']), max(gene_info['start_pos'], gene_info['end_pos'])
		if end<amr_start or start> amr_end:
			continue
		else:
			# added by 1 because in string indecesstarts from 0
			diff = max((amr_start+1-start), 0)+max((end - (amr_end+1)), 0)
			if ((1-(float(diff)/(end-start)))*100)>most_overlap:
				most_overlap = (1-(float(diff)/(end-start)))*100
				amr_coverage = gene_info['coverage']
				amr_index = i
				error = False
	return amr_coverage, amr_index, error

def check_coverage_consistency_remove_rest_seq(seq_info_list_input,
								coverage_thr, amr_name, annotate_dir):
	"""
	To compare the coverage of each gene in the neighborhood with that of the AMR
	and remove the gene (and any gene before that in case of upstream and any gene
	after that in case of downstream) if the difference between this gene's coverage
	and ame coverage is more than covearge_thr
	Parameters:
		seq_info_list_input: the list of all sequences and their info including annotations
		coverage_thr:	the threshold used to compare the gene covearge and amr covearge
		amr_name:	the name of target AMR
		annotate_dir:	the directory in which annotation info are stored
	Return:
	"""
	seq_info_list = copy.deepcopy(seq_info_list_input)
	#extract amr info
	amr_coverages = []
	amr_indeces = []
	for seq_info in seq_info_list:
		found_amr = False
		for gene_counter, gene_info in enumerate(seq_info):
			if gene_info['coverage'] is None:
				logging.info("Coverage information are not available for "+amr_name)
				return "", 0
			coverage = round(gene_info['coverage'], 2)
			if gene_info['target_amr']=='yes':
				amr_coverages.append(coverage)
				amr_indeces.append(gene_counter)
				found_amr = True
				break
		if not found_amr:
			amr_coverage, amr_index, error = find_target_amr_in_seqvalue_and_return_coverage(seq_info)
			if error:
				logging.info("ERROR: no target amr was found for "+ str(seq_info)+" regarding "+amr_name)
				import pdb; pdb.set_trace()
				sys.exit()
			else:
				amr_coverages.append(amr_coverage)
				amr_indeces.append(amr_index)
	if len(amr_coverages)!=len(seq_info_list):
		logging.error("ERROR: inconsistency between the number of sequences and found amr-coverages!")
		import pdb; pdb.set_trace()
	#check coverage consistency by comparing its coverage with AMR coverage
	# and remove genes with inconsistent coverage and whatever comes before them if upstream OR after them if downstream
	remained_seqs = []
	for i, seq_info in enumerate(seq_info_list):
		#find the genes need to be removed
		to_be_removed_genes=[]
		for j, gene_info in enumerate(seq_info):
			if abs(gene_info['coverage'] - amr_coverages[i])>coverage_thr:
				if j<amr_indeces[i]:
					for k in range(j+1):
						if k not in to_be_removed_genes:
							to_be_removed_genes.append(k)
				elif j>amr_indeces[i]:
					for k in range(j, len(seq_info)):
						if k not in to_be_removed_genes:
							to_be_removed_genes.append(k)
					break
		for j in reversed(range(len(seq_info))):
			if j in to_be_removed_genes:
				del seq_info[j]
		#check if the remained sequence already exists in the seq_info_list
		if seq_info and not similar_seq_annotation_already_exist(seq_info, remained_seqs, annotate_dir):
			remained_seqs.append(seq_info)

	#Initialize coverage file
	coverage_annotation = os.path.join(annotate_dir,'coverage_annotation_'+str(coverage_thr)+'_'+amr_name+'.csv')
	with open(coverage_annotation,'w') as fd:
		writer = csv.writer(fd)
		writer.writerow(['seq_name', 'seq_value', 'seq_length', 'gene',
							'coverage', 'length', 'start_pos', 'end_pos', 'target_amr'])
		#write extracted sequences with consistent coverage
		for seq_info in remained_seqs:
			for gene_info in seq_info:
				writer.writerow([gene_info['seq_name'], gene_info['seq_value'],
							len(gene_info['seq_value']),
							gene_info['gene'], gene_info['coverage'],
							gene_info['length'], gene_info['start_pos'],
							gene_info['end_pos'], gene_info['target_amr']])

	return coverage_annotation, len(remained_seqs)

def extract_seq_annotation(annotate_dir, prokka_prefix, use_RGI, RGI_include_loose,
								seq_pair):
	"""
	The function used in parallel anntation to call the function for annotating a sequence
	Parameters:
		annotate_dir: the directory to store annotation output
		prokka_prefix: the first part of prokka command probably is nonempty only when docker is used
		use_RGI: if True we want to call RGI for annotating AMRs
		RGI_include_loose: if True loose mode is used
		seq_pair: the index and the value of a sequence to be annotated
	Return:
		the list of annotated genes
	"""
	counter, ext_seq = seq_pair
	seq_description = 'extracted'+str(counter)
	seq_info_list = annotate_sequence(ext_seq, seq_description, annotate_dir,
									prokka_prefix, use_RGI, RGI_include_loose)
	return seq_info_list

def extract_graph_seqs_annotation_parallel(amr_name, path_info_file, neighborhood_seq_file,
					annotate_dir, core_num, prokka_prefix,
					use_RGI, RGI_include_loose,
					annotation_writer, trimmed_annotation_writer, gene_file,
					product_file, error_file):
	"""
	To annotate neighborhood sequences of AMR extracted from the graph in parallel
	Parameters:
		amr_name: the name of AMR
		path_info_file: the information of nodes representing the extracted sequence
		neighborhood_seq_file: the file containing extracted sequences
		annotate_dir: the directory to sore annotation results
		core_num: the number of core for parallel processing
		prokka_prefix: the prefix in prokka command mostly used in docker
		use_RGI: if True RGI is used for AMR annotation
		RGI_include_loose: if True use loose mode in RGI
		annotation_writer: the file to store annotation results
		trimmed_annotation_writer: the file to store unique annotation results
		gene_file: the file to store gene nams in annotation
		product_file: the file to store product name in annotation
		error_file: the file to store errors
	Return:
		the list of annotated genes and their details
	"""
	error_writer = open(error_file, 'a')
	#Read path_info from the file
	path_info_list = []
	if path_info_file!=-1:
		path_info_list = read_path_info_file(path_info_file)
	#find the list of all extracted sequences
	logging.info('Reading '+ neighborhood_seq_file + ' for '+ amr_name)
	sequence_list = []
	counter = 1
	with open(neighborhood_seq_file, 'r') as read_obj:
		for line in read_obj:
			if line.startswith('>') or line.startswith('Path') or line.startswith('The'):
				continue
			sequence_list.append((counter, line))
			counter+=1
	#Parallel annotation
	p_annotation = partial(extract_seq_annotation, annotate_dir, prokka_prefix,
							use_RGI, RGI_include_loose)
	with Pool(core_num) as p:
		seq_info_list = p.map(p_annotation, sequence_list)
	#Further processing of result of parallel annotation
	all_seq_info_list =[]
	# if amr_name=='GES-4' or amr_name=='GES-21':
	# 	import pdb; pdb.set_trace()
	for i, seq_pair in enumerate(sequence_list):
		counter, line = seq_pair
		seq_description = 'extracted'+str(counter)
		seq_info = seq_info_list[i]
		#extract amr from seq_info
		amr_found, amr_info, up_info, down_info, seq_info = split_up_down_info(line[:-1], seq_info)
		if not amr_found:
			logging.error("ERROR: no target amr was found in the extracted sequence")
			error_writer.write(amr_name+' annotation not found! '+" seq_info: "+str(seq_info)+'\n')
			continue
			#import pdb; pdb.set_trace()
		#calculate coverage for the genes available in the annotation
		coverage_list = []
		if path_info_list:
			coverage_list = find_gene_coverage(seq_info, path_info_list[counter-1])
		#Check if this annotation has already been found
		found = seq_annotation_already_exist(seq_info, all_seq_info_list, annotate_dir)
		#If it's a novel sequence correct annotation (if applicable) for cases that RGI doesn't have a hit but Prokka has
		if not found:
			all_seq_info_list.append(seq_info)
		myLine1 = myLine2 = seq_description +':\t'
		#write annotation onfo into the files
		for j, gene_info in enumerate(seq_info):
			coverage = coverage_list[j] if coverage_list else -1
			gene_info['coverage'] = coverage
			gene_info['seq_name'] = seq_description
			write_info_in_annotation_file(annotation_writer, trimmed_annotation_writer,
										gene_info, use_RGI, found)
			if gene_info['gene']=='':
				myLine1+='UNKNOWN---'
			else:
				myLine1+=gene_info['gene']+'---'
			myLine2+=gene_info['product']+'---'
		gene_file.write(myLine1[:-3]+'\n')
		product_file.write(myLine2[:-3]+'\n')
	if not all_seq_info_list:
		error_writer.write(amr_name+' no annotation was found in the graph.\n')
	error_writer.close()
	return all_seq_info_list

def neighborhood_annotation_parallel(amr_name, neighborhood_seq_file,
								path_info_file, seq_length,
								output_dir, prokka_prefix, use_RGI = True,
								RGI_include_loose = False, output_name ='',
								core_num = 4):
	"""
	To annotate reference genomes (a piece extracted around the AMR gene) as well as
		extracted neighborhood sequences from assembly graph, summarize the results
		in a couple of formats and visualize them
	Parameters:
		amr_name:	the name of target AMR
		neighborhood_seq_file:	the address of the file containing all extracted
		 	neighborhood sequences from assembly graph
		seq_length:	the length of neighborhood sequence extracted from each side of
			the AMR sequence (downstream and up stream)
		output_dir:	the path for the output directory
		prokka_prefix: it's used to run prokka properly via docker or conda
		use_RGI:	RGI annotations incorporated for AMR annotation
		RGI_include_loose: Whether to include loose annotaions in RGI
		output_name:the name used to distinguish different output files usually based on the name of AMR
	Return:
		the address of files stroing annotation information (annotation_detail_name,
			trimmed_annotation_info, gene_file_name, product_file_name, visual_annotation)
	"""
	logging.debug('Started annotation for '+amr_name)
	# initializing required files and directories
	annotate_dir = os.path.join(output_dir,ANNOTATION_DIR, ANNOTATION_DIR+'_'+str(seq_length),'annotation'+output_name+'_'+str(seq_length))
	if os.path.exists(annotate_dir):
		try:
			shutil.rmtree(annotate_dir)
		except OSError as e:
			logging.error("Error: %s - %s." % (e.filename, e.strerror))
	os.makedirs(annotate_dir)
	error_file =os.path.join(output_dir, ANNOTATION_DIR, ANNOTATION_DIR+'_'+str(seq_length), "not_found_annotation_amrs_in_graph.txt")
	annotation_detail_name = os.path.join(annotate_dir,'annotation_detail'+output_name+'.csv')
	trimmed_annotation_info_name = os.path.join(annotate_dir, 'trimmed_annotation_info'+output_name+'.csv')
	annotation_detail = open(annotation_detail_name, mode='w', newline='')
	trimmed_annotation_info = open(trimmed_annotation_info_name, mode='w', newline='')
	annotation_writer = csv.writer(annotation_detail)
	trimmed_annotation_writer = csv.writer(trimmed_annotation_info)
	gene_info = {'seq_value':'seq_value', 'gene':'gene', 'prokka_gene_name':'prokka_gene_name',
				'product':'product', 'length':'length', 'start_pos':'start_pos',
				'end_pos':'end_pos', 'RGI_prediction_type':'RGI_prediction_type',
				'coverage':'coverage', 'family':'family', 'seq_name':'seq_name',
				'target_amr':'target_amr'}
	write_info_in_annotation_file(annotation_writer, trimmed_annotation_writer,
								gene_info, use_RGI, False, 'seq_length')
	gene_file_name = os.path.join(annotate_dir, 'seq_comparison_genes'+output_name+'.txt')
	gene_file = open(gene_file_name, 'w')
	product_file_name = os.path.join(annotate_dir,'seq_comparison_products'+output_name+'.txt')
	product_file = open(product_file_name, 'w')

	#annotate the sequences extraced from assembly graph
	all_seq_info_list = extract_graph_seqs_annotation_parallel(amr_name, path_info_file, neighborhood_seq_file,
									annotate_dir, core_num,
									prokka_prefix, use_RGI, RGI_include_loose,
									annotation_writer, trimmed_annotation_writer,
									gene_file, product_file, error_file)
	logging.info("NOTE: The comparison of neighborhood sequences are available in " +\
	 		annotation_detail_name+", "+gene_file_name+", "+product_file_name)
	annotation_detail.close()
	trimmed_annotation_info.close()
	gene_file.close()
	product_file.close()

	return all_seq_info_list, trimmed_annotation_info_name

def extract_graph_seqs_annotation(amr_name, path_info_file, neighborhood_seq_file,
					annotate_dir, prokka_prefix,
					use_RGI, RGI_include_loose,
					annotation_writer, trimmed_annotation_writer, gene_file, product_file):
	"""
	To annotate neighborhood sequences of AMR extracted from the graph
	Parameters:
		amr_name: the name of AMR
		path_info_file: the information of nodes representing the extracted sequence
		neighborhood_seq_file: the file containing extracted sequences
		annotate_dir: the directory to sore annotation results
		prokka_prefix: the prefix in prokka command mostly used in docker
		use_RGI: if True RGI is used for AMR annotation
		RGI_include_loose: if True use loose mode in RGI
		annotation_writer: the file to store annotation results
		trimmed_annotation_writer: the file to store unique annotation results
		gene_file: the file to store gene nams in annotation
		product_file: the file to store product name in annotation
		error_file: the file to store errors
	Return:
		the list of annotated genes and their details
	"""
	error_file ="not_found_annotation_amrs_in_graph.txt"
	error_writer = open(error_file, 'a')
	counter = 1
	all_seq_info_list =[]
	#Read path_info from the file
	path_info_list = []
	if path_info_file!=-1 and path_info_file!='':
		path_info_list = read_path_info_file(path_info_file)
	logging.info('Reading '+ neighborhood_seq_file + ' for '+ amr_name)
	with open(neighborhood_seq_file, 'r') as read_obj:
		for line in read_obj:
			if line.startswith('>') or line.startswith('Path') or line.startswith('The'):
				continue
			seq_description = 'extracted'+str(counter)
			seq_info = annotate_sequence(line, seq_description, annotate_dir,
											prokka_prefix, use_RGI, RGI_include_loose)
			#extract amr from seq_info
			amr_found, amr_info, up_info, down_info, seq_info = split_up_down_info(line[:-1], seq_info)
			if not amr_found:
				logging.error("ERROR: no target amr was found in the extracted sequence")
				error_writer.write(amr_name+' annotation not found! '+" seq_info: "+str(seq_info)+'\n')
				counter+=1
				continue
				#import pdb; pdb.set_trace()
			#calculate the coverage of annotated genes
			coverage_list = []
			if path_info_list:
				coverage_list = find_gene_coverage(seq_info, path_info_list[counter-1])
			#Check if this annotation has already been found
			found = seq_annotation_already_exist(seq_info, all_seq_info_list, annotate_dir)
			#If it's a novel sequence correct annotation if possible
			if not found:
				all_seq_info_list.append(seq_info)
			myLine1 = myLine2 = seq_description +':\t'
			#write annotation onfo into the files
			for j, gene_info in enumerate(seq_info):
				coverage = coverage_list[j] if coverage_list else -1
				gene_info['coverage'] = coverage
				gene_info['seq_name'] = seq_description
				write_info_in_annotation_file(annotation_writer, trimmed_annotation_writer,
											gene_info, use_RGI, found)
				if gene_info['gene']=='':
					myLine1+='UNKNOWN---'
				else:
					myLine1+=gene_info['gene']+'---'
				myLine2+=gene_info['product']+'---'
			gene_file.write(myLine1[:-3]+'\n')
			product_file.write(myLine2[:-3]+'\n')
			counter+=1
	if not all_seq_info_list:
		error_writer.write(amr_name+' no annotation was found in the graph.\n')
	error_writer.close()
	return all_seq_info_list

def neighborhood_annotation(amr_name, neighborhood_seq_file,
								path_info_file, seq_length,
								output_dir, prokka_prefix, use_RGI = True,
								RGI_include_loose = False, output_name =''):
	"""
	To annotate reference genomes (a piece extracted around the AMR gene) as well as
		extracted neighborhood sequences from assembly graph, summarize the results
		in a couple of formats and visualize them.
	Parameters:
		amr_name:	the name of target AMR
		neighborhood_seq_file:	the address of the file containing all extracted
		 	neighborhood sequences from assembly graph
		seq_length:	the length of neighborhood sequence extracted from each side of
			the AMR sequence (downstream and up stream)
		output_dir:	the path for the output directory
		prokka_prefix: it's used to run prokka properly via docker or conda
		use_RGI:	RGI annotations incorporated for AMR annotation
		RGI_include_loose: Whether to include loose annotaions in RGI
		output_name:the name used to distinguish different output files usually based on the name of AMR
	Return:
		the address of files storing annotation information (annotation_detail_name,
			trimmed_annotation_info, gene_file_name, product_file_name, visual_annotation)
	"""
	logging.debug('Started annotation for '+amr_name)
	# initializing required files and directories
	annotate_dir = os.path.join(output_dir,ANNOTATION_DIR, ANNOTATION_DIR+'_'+str(seq_length), 'annotation'+output_name+'_'+str(seq_length))
	if os.path.exists(annotate_dir):
		try:
			shutil.rmtree(annotate_dir)
		except OSError as e:
			logging.error("Error: %s - %s." % (e.filename, e.strerror))
	os.makedirs(annotate_dir)
	annotation_detail_name = os.path.join(annotate_dir,'annotation_detail'+output_name+'.csv')
	trimmed_annotation_info_name = os.path.join(annotate_dir,'trimmed_annotation_info'+output_name+'.csv')
	annotation_detail = open(annotation_detail_name, mode='w', newline='')
	trimmed_annotation_info = open(trimmed_annotation_info_name, mode='w', newline='')
	#annotation_writer = csv.writer(annotation_detail, delimiter=',', quoting=csv.QUOTE_MINIMAL)
	annotation_writer = csv.writer(annotation_detail)
	trimmed_annotation_writer = csv.writer(trimmed_annotation_info)
	gene_info = {'seq_value':'seq_value', 'gene':'gene', 'prokka_gene_name':'prokka_gene_name',
				'product':'product', 'length':'length', 'start_pos':'start_pos',
				'end_pos':'end_pos', 'RGI_prediction_type':'RGI_prediction_type',
				'coverage':'coverage', 'family':'family', 'seq_name':'seq_name',
				'target_amr':'target_amr'}
	write_info_in_annotation_file(annotation_writer, trimmed_annotation_writer,
								gene_info, use_RGI, False, 'seq_length')
	gene_file_name = os.path.join(annotate_dir,'seq_comparison_genes'+output_name+'.txt')
	gene_file = open(gene_file_name, 'w')
	product_file_name = os.path.join(annotate_dir,'seq_comparison_products'+output_name+'.txt')
	product_file = open(product_file_name, 'w')
	#annotate the sequences extraced from assembly graph
	all_seq_info_list = extract_graph_seqs_annotation(amr_name, path_info_file, neighborhood_seq_file,
					annotate_dir, prokka_prefix,
					use_RGI, RGI_include_loose,
					annotation_writer, trimmed_annotation_writer, gene_file, product_file)
	annotation_detail.close()
	trimmed_annotation_info.close()
	gene_file.close()
	product_file.close()
	logging.info("NOTE: The comparison of neighborhood sequences are available in " +\
	 		annotation_detail_name+", "+gene_file_name+", "+product_file_name)

	return all_seq_info_list, trimmed_annotation_info_name

def is_there_amr_in_graph(gfa_file, output_dir, bandage_path, threshold, amr_file):
	"""
	To call bandage+blast and check if the amr sequence can be found in the assembly graph
	Parameters:
		amr_file: the address of the query file
		amr_name: the name of the AMR sequence
		gfa_file: the address of the assembly graph
		output_dir: the address of the output directory
		bandage_path: the address of bandage executable file
		threshold: the threshold for coverage and identity
	Return:
		a boolean value which is True if amr_file was found in gfa_file and the
		list of paths returned by bandage+blast in which the coverage and identiry
		are greater/equal than/to threshold
	"""
	amr_name = extract_name_from_file_name(amr_file)
	logging.info('Checking if AMR "'+amr_name+'" exists in the assembly graph...')
	output_name=os.path.join(output_dir,amr_name+'_align_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
	if os.path.isfile(output_name+'.tsv'):
		os.remove(output_name+'.tsv')
	bandage_command = subprocess.run([bandage_path, "querypaths", gfa_file, amr_file,
						output_name, "--pathnodes", "50", "--minpatcov",
						str((threshold-1)/100.0), "--minmeanid", str((threshold-1)/100.0),
						"--minhitcov", str((threshold-1)/100.0)],
						stdout=subprocess.PIPE, check= True )
	logging.info(bandage_command.stdout.decode('utf-8'))
	# command = bandage_path +' querypaths '+gfa_file+' '+amr_file+' '+output_name + ' --pathnodes 50'
	# os.system(command)

	found, paths_info = read_path_info_from_align_file(output_name+".tsv", threshold)
	if not found:
		os.remove(output_name+'.tsv')
		logging.debug(amr_name+' not found in the graph!')
	else:
		logging.debug(amr_name+' found!')
	return found, paths_info, output_name+".tsv"

def amr_path_overlap(found_amr_paths, new_paths, new_amr_len, overlap_percent = 95):
	"""
	To check if all paths found for the new AMR seq overlap significantly (greater/equal
	 than/to overlap percent) with the already found paths for other AMRs
	Parameters:
	 	found_amr_paths:  	the paths already found for AMR genes
		new_paths: 			the paths found for the new AMR gene
		overlap_percent:	the threshold for overlap
	Return:
		False only if every paths in new_paths have overlap with at least one path in found_amr_paths
		True if we can find at least one path that is unique and not available in found_amr_paths
		Also, it returns the list of indeces from found_amr_paths that had overlap with a path in new_paths
	"""
	id_list = []
	for new_path in new_paths:
		found = False
		for i, paths in enumerate(found_amr_paths):
			for path in paths:
				if path['nodes'] == new_path['nodes'] and path['orientations']==new_path['orientations']:
					# for now we just check overlaps when they are in the same node(s)
					diff_length = max(path['start_pos']-new_path['start_pos'],0)+max(new_path['end_pos']-path['end_pos'],0)
					percent = (1 - (float(diff_length)/(new_amr_len)))*100
					if percent >= overlap_percent:
						found = True
						if i not in id_list:
							id_list.append(i)
						break
			if found:
				break
	if len(id_list)==len(new_paths):
		return True, id_list
	return False, None

def are_there_amrs_in_graph(gfa_file, output_dir, bandage_path, threshold, amr_object):
	"""
	To call bandage+blast and check if the amr sequence can be found in the assembly graph
	Parameters:
		amr_file: the address of the query file
		amr_name: the name of the AMR sequence
		gfa_file: the address of the assembly graph
		output_dir: the address of the output directory
		bandage_path: the address of bandage executable file
		threshold: the threshold for coverage and identity
	Return:
		a boolean value which is True if amr_file was found in gfa_file and the
		list of paths returned by bandage+blast in which the coverage and identiry
		are greater/equal than/to threshold
	"""
	cat_file, amr_files = amr_object
	amr_names = [extract_name_from_file_name(e) for e in amr_files]
	logging.info('Checking if AMRs "'+str(amr_names)+'" exists in the assembly graph...')
	output_name=os.path.join(output_dir, extract_name_from_file_name(cat_file)+'_align_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
	if os.path.isfile(output_name+'.tsv'):
		os.remove(output_name+'.tsv')
	bandage_command = subprocess.run([bandage_path, "querypaths", gfa_file, cat_file,
						output_name, "--pathnodes", "50", "--minpatcov",
						str((threshold-1)/100.0), "--minmeanid", str((threshold-1)/100.0),
						"--minhitcov", str((threshold-1)/100.0)],
						stdout=subprocess.PIPE, check= True )
	logging.info(bandage_command.stdout.decode('utf-8'))
	# command = bandage_path +' querypaths '+gfa_file+' '+cat_file+' '+output_name + ' --pathnodes 50'
	# os.system(command)

	paths_info_list = read_path_info_from_align_file_with_multiple_amrs(output_name+".tsv", threshold)

	return paths_info_list

def process_amr_group_and_find(gfa_file, align_dir, output_dir, bandage_path,
										amr_threshold, amr_object):
	"""
	Read a group of AMRs, write them into a single file and call bandage+blast for it
	to identify them in the graph
	Parameters:
		gfa_file: the file containing the assembly graph
		align_dir: the directory for storing alignment info
		output_dir: the directory to store the list of AMRs in a single file
		bandage_path: the path for bandage
		amr_threshor: the threshold for identity and coverage used in alignment
		amr_object: the list of AMRs and their ids
	Return:
		the alignment info for AMRs
	"""
	g_id, amr_group = amr_object
	#read info of the group into a single file
	cat_file = os.path.join(output_dir,AMR_DIR_NAME, 'amr_group_'+str(g_id)+'.fasta')
	file_group = []
	with open(cat_file, 'w') as writer:
		for amr_info in amr_group:
			amr_seq, amr_title = amr_info
			writer.write(amr_title)
			writer.write(amr_seq)
			amr_name1 = amr_name_from_comment(amr_title)
			amr_file_name = restricted_amr_name_from_modified_name(amr_name1)
			file_group.append(amr_file_name+'.fasta')

	#Run Bandage+BLAST
	p_find_amr_align = are_there_amrs_in_graph(gfa_file, align_dir,
										bandage_path, amr_threshold, (cat_file, file_group))
	#Remove temporary AMR file
	if os.path.isfile(cat_file):
		os.remove(cat_file)
	return p_find_amr_align

def find_all_amr_in_graph_parallel(gfa_file, output_dir, amr_sequences_file,
									bandage_path, amr_threshold, core_num):
	"""
	To go over a list of AMR sequences (amr_sequences_file) and run bandage+blast
	to check if any of them exists in the assembly graph (gfa_file)
	Parameters:
		gfa_file: the address of the assembly graph
		output_dir: the address of the output directory
		amr_sequences_file: the address of the file containing the sequence of all AMRs from CARD
		bandage_path: the address of bandage executable file
		amr_threshold: the threshold for coverage and identity
		core_num: the number of used cores
	Return:

	"""
	align_dir = os.path.join(output_dir, AMR_DIR_NAME, AMR_ALIGN_DIR)
	if not os.path.exists(align_dir):
		os.makedirs(align_dir)

	#generate the groups and store the group of each amr
	group_num = 5
	amr_group_id = collections.defaultdict(list)
	amr_file_groups = [[] for i in range(group_num*core_num)]
	amr_title = ''
	amr_seq_title_list = []
	#Read AMR sequences one by one
	amr_counter = 0
	with open(amr_sequences_file) as fp:
		for line in fp:
			if line.startswith('>'):
				amr_title = line
				continue
			amr_name = amr_name_from_comment(amr_title[:-1])
			amr_seq_title_list.append((line, amr_title))
			id = amr_counter % (group_num*core_num)
			amr_file_groups[id].append((line, amr_title))
			amr_group_id[amr_name] = id
			amr_counter+=1

	amr_objects = [(i, e) for i,e in enumerate(amr_file_groups)]
	#parallel run Bandage+BLAST
	p_find_amr = partial(process_amr_group_and_find, gfa_file, align_dir,
							output_dir, bandage_path, amr_threshold)
	with Pool(core_num) as p:
		paths_info_group_list = p.map(p_find_amr, amr_objects)


	unique_amr_seqs = []
	unique_amr_infos = []
	unique_amr_paths = []
	#process the result of parallel processes
	for i, amr_object in enumerate(amr_seq_title_list):
		amr_name = amr_name_from_comment(amr_object[1])
		id = amr_group_id[amr_name]
		restricted_amr_name = restricted_amr_name_from_modified_name(amr_name)
		if restricted_amr_name in paths_info_group_list[id]:
			logging.debug(amr_name + ' was found!')
			path_info = paths_info_group_list[id][restricted_amr_name]
			overlap, amr_ids =  amr_path_overlap(unique_amr_paths, path_info,
													len(amr_object[0])-1)
			if not overlap:
				unique_amr_seqs.append(amr_object[0])
				amr_info = {'name':amr_object[1], 'overlap_list':[]}
				unique_amr_infos.append(amr_info)
				unique_amr_paths.append(path_info)
			else:
				if len(amr_ids)>1:
					logging.error("an AMR has overlap with more than one group")
					import pdb; pdb.set_trace()
				# add this AMR to the right group of AMRs all having overlaps
				for id in amr_ids:
					if amr_name not in unique_amr_infos[id]['overlap_list']:
						unique_amr_infos[id]['overlap_list'].append(amr_name)

	# write information (the sequence of found AMRs that don't have overlaped paths with others
	# + the list of groups in each all AMRs have overlaped paths) into files
	AMR_dir = os.path.join(output_dir, AMR_DIR_NAME, AMR_SEQ_DIR)
	if not os.path.exists(AMR_dir):
		os.makedirs(AMR_dir)
	overlap_file_name = os.path.join(output_dir, AMR_DIR_NAME, AMR_OVERLAP_FILE)
	overlap_file = open(overlap_file_name, 'w')
	unique_amr_files = []
	for i, seq in enumerate(unique_amr_seqs):
		amr_name = amr_name_from_comment(unique_amr_infos[i]['name'])
		restricted_amr_name = restricted_amr_name_from_modified_name(amr_name)
		amr_file = create_fasta_file(seq, AMR_dir, unique_amr_infos[i]['name'], restricted_amr_name)
		unique_amr_files.append(amr_file)
		overlap_file.write(amr_name+":")
		if unique_amr_infos[i]['overlap_list']:
			overlap_file.write(', '.join(e for e in unique_amr_infos[i]['overlap_list']))
			overlap_file.write("\n")
		else:
			overlap_file.write("\n")
	overlap_file.close()

	return unique_amr_files, unique_amr_paths

def find_all_amr_in_graph(gfa_file, output_dir, amr_sequences_file, bandage_path, amr_threshold):
	"""
	To go over a list of AMR sequences (amr_sequences_file) and run bandage+blast
	to check if any of them exists in the assembly graph (gfa_file)
	Parameters:
		gfa_file: the address of the assembly graph
		output_dir: the address of the output directory
		amr_sequences_file: the address of the file containing the sequence of all AMRs from CARD
		bandage_path: the address of bandage executable file
		amr_threshold: the threshold for coverage and identity
	Return:

	"""
	align_dir = os.path.join(output_dir,AMR_DIR_NAME, AMR_ALIGN_DIR)
	if not os.path.exists(align_dir):
		os.makedirs(align_dir)

	amr_name = ''
	#found_amr_names = []
	unique_amr_seqs = []
	unique_amr_infos = []
	unique_amr_paths = []
	#Read AMR sequences one by one
	with open(amr_sequences_file) as fp:
		for line in fp:
			if line.startswith('>'):
				amr_title = line
				continue
			#create a fasta file for it
			amr_file = create_fasta_file(line, output_dir)
			#Run Bandage+BLAST
			amr_name = amr_name_from_comment(amr_title)
			found, paths_info, tsv_file = is_there_amr_in_graph(gfa_file, align_dir,
											bandage_path, amr_threshold, amr_file)
			if found:
				#found_amr_names.append(amr_name)
				overlap, amr_ids =  amr_path_overlap(unique_amr_paths, paths_info,
										len(line)-1)
				if not overlap:
					unique_amr_seqs.append(line)
					amr_info = {'name':amr_title, 'overlap_list':[]}
					unique_amr_infos.append(amr_info)
					unique_amr_paths.append(paths_info)
				else:
					if len(amr_ids)>1:
						logging.error("an AMR has overlap with more than one group")
						import pdb; pdb.set_trace()
					# add this AMR to the right group of AMRs all having overlaps
					for id in amr_ids:
						if amr_name not in unique_amr_infos[id]['overlap_list']:
							unique_amr_infos[id]['overlap_list'].append(amr_name)
	# write information (the sequence of found AMRs that don't have overlaped paths with others
	# + the list of groups in each all AMRs have overlaped paths) into files
	AMR_dir = os.path.join(output_dir,AMR_DIR_NAME, AMR_SEQ_DIR)
	if not os.path.exists(AMR_dir):
		os.makedirs(AMR_dir)
	overlap_file_name = os.path.join(output_dir, AMR_DIR_NAME, AMR_OVERLAP_FILE)
	overlap_file = open(overlap_file_name, 'w')
	unique_amr_files = []
	for i, seq in enumerate(unique_amr_seqs):
		amr_name = amr_name_from_comment(unique_amr_infos[i]['name'])
		restricted_amr_name = restricted_amr_name_from_modified_name(amr_name)
		amr_file = create_fasta_file(seq, AMR_dir, unique_amr_infos[i]['name'], restricted_amr_name)
		unique_amr_files.append(amr_file)
		overlap_file.write(amr_name+":")
		if unique_amr_infos[i]['overlap_list']:
			#overlap_file.write(amr_name+":\n")
			overlap_file.write(', '.join(e for e in unique_amr_infos[i]['overlap_list']))
			overlap_file.write("\n")
		else:
			overlap_file.write("\n")
	overlap_file.close()

	return unique_amr_files, unique_amr_paths

def find_corrsponding_seq_path_file(amr_name, sequences_file_names, path_info_file_names, seq_length):
	"""
	To return the name of a sequence file (from sequences_file_names)
	and the name of a path file (from path_info_file_names)
	dedicated to sequences extracted with a given length (seq_length)
	from a given amr sequence (amr_name)
	"""
	seq_file = -1
	for file_name in sequences_file_names:
		if SEQ_NAME_PREFIX+amr_name+'_'+str(seq_length) in file_name:
			seq_file = file_name
	path_file = -1
	for file_name in path_info_file_names:
		if SEQ_NAME_PREFIX+amr_name+'_'+str(seq_length) in file_name:
			path_file = file_name
	return seq_file, path_file

def seq_annotation_trim_main(params, amr_files, all_seq_info_lists,
								annotation_files, visualize = False):
	"""
	The core function to filter and remove genes that their coverage difference
	from AMR coverage is above a threshold
	Prameters:
		params: the list of parameters imported from params.py
		amr_files: the list of files containing AMRs
		all_seq_info_lists: the list of annotations of neighborhood sequences extracted from the graph
		annotation_files: the files containing annotation info
		visualize: if True, visualize
	"""
	coverage_annotation_list = []
	for i, amr_file in enumerate(amr_files):
		restricted_amr_name = extract_name_from_file_name(amr_file)
		#remove some extracted sequences based on coverage consistency
		annotate_dir = os.path.join(params.output_dir,ANNOTATION_DIR, ANNOTATION_DIR+'_'+\
			str(params.seq_length), 'annotation_'+restricted_amr_name+'_'+str(params.seq_length))
		coverage_annotation = ''
		remained_seqs = []
		if params.coverage_thr>0:
			coverage_annotation, remained_seqs = check_coverage_consistency_remove_rest_seq(\
								all_seq_info_lists[i],
								params.coverage_thr, restricted_amr_name, annotate_dir)
		if visualize:
			# create an image presenting the annotations for all sequences
			if coverage_annotation!='':
				visual_annotation_csv = coverage_annotation
			else:
				visual_annotation_csv = annotation_files[i]
			visual_annotation = os.path.join(annotate_dir, 'gene_comparison_'+str(params.coverage_thr)+'_'+restricted_amr_name+'.png')
			visualize_annotation(visual_annotation_csv, output=visual_annotation)
		if params.coverage_thr>0:
			coverage_annotation_list.append(coverage_annotation)
	return coverage_annotation_list

def seq_annotation_main(params, seq_files, path_info_files, amr_files):
	"""
	The core function for annotation of neighborhood sequences of all AMRs
	Parameters:
		params: the list of parameters extracted from params.py
		seq_files: the list of neighborhood sequence files
		path_info_files: the list of files containing node info for all sequences
		amr_files: the list of files containing AMRs
	Return:

	"""
	logging.info("Neighborhood Annotation ...")
	if seq_files:
		neighborhood_files = seq_files
	else:
		neighborhood_files = extract_files(params.ng_seq_files, 'please provide the \
			address of the files containing all extracted sequences from AMR neighborhood \
			in the assembly graph')
	if path_info_files:
		nodes_info_files = path_info_files
	else:
		nodes_info_files = extract_files(params.ng_path_info_files, '')
	#extract ref neighborhood annotation from the file

	all_seq_info_lists = []
	annotation_files = []
	for amr_file in amr_files:
		restricted_amr_name = extract_name_from_file_name(amr_file)
		_, amr_name = retrieve_AMR(amr_file)
		neighborhood_file, nodes_info_file = find_corrsponding_seq_path_file(restricted_amr_name,
								neighborhood_files, nodes_info_files, params.seq_length)
		if neighborhood_file == -1:
			logging.error('no sequence file for '+ amr_file +' was found! We looked for a file like '+restricted_amr_name)
			import pdb; pdb.set_trace()
			sys.exit()
		if params.multi_processor:
			all_seq_info_list, annotation_file=\
				neighborhood_annotation_parallel(amr_name, neighborhood_file,
					nodes_info_file, params.seq_length,
					params.output_dir, params.PROKKA_COMMAND_PREFIX,params.use_RGI,
					params.RGI_include_loose, '_'+restricted_amr_name,
					params.core_num)
		else:
			all_seq_info_list, annotation_file =\
				neighborhood_annotation(amr_name, neighborhood_file,
					nodes_info_file, params.seq_length,
					params.output_dir, params.PROKKA_COMMAND_PREFIX,params.use_RGI,
					params.RGI_include_loose, '_'+restricted_amr_name)
		all_seq_info_lists.append(all_seq_info_list)
		annotation_files.append(annotation_file)

	return all_seq_info_lists, annotation_files

def sequence_neighborhood_main(params, gfa_file, amr_seq_align_info):
	"""
	The core function to extract the neighborhood of AMRs
	Parameters:
		params: the list pf parameters imported from params.py
		gfa_file: the file containing the assembly graph
		amr_seq_align_info: the alignment info (AMR alignment in the graph)
	Return:
		the list of files containing extracted sequence and the details of nodes representing them
	"""
	seq_files = []
	path_info_files = []
	logging.info("Extracting neighborhood sequences with length = %s", params.seq_length)

	#remove paths from GFA file
	delete_lines_started_with('P', gfa_file)
	sequence_dir = os.path.join(params.output_dir, SEQ_DIR_NAME, SEQ_DIR_NAME+'_'+str(params.seq_length))
	if not os.path.exists(sequence_dir):
		os.makedirs(sequence_dir)
	if params.multi_processor and params.ng_extraction_time_out<=0:
		p_extraction = partial(neighborhood_sequence_extraction, gfa_file, params.seq_length,
							sequence_dir, params.BANDAGE_PATH,
							params.amr_identity_threshold, SEQ_NAME_PREFIX,
							params.path_node_threshold ,
							params.max_kmer_size, params.ng_extraction_time_out, params.assembler)
		with Pool(params.core_num) as p:
			lists = p.map(p_extraction, amr_seq_align_info)
		seq_files, path_info_files = zip(*lists)
	else:
		for amr_file in amr_seq_align_info:
			seq_file, path_info_file = neighborhood_sequence_extraction(gfa_file, params.seq_length,
								sequence_dir, params.BANDAGE_PATH,
								params.amr_identity_threshold, SEQ_NAME_PREFIX,
								params.path_node_threshold ,
								params.max_kmer_size, params.ng_extraction_time_out,
								params.assembler, amr_file)
			if seq_file:
				path_info_files.append(path_info_file)
				seq_files.append(seq_file)

	return seq_files, path_info_files

def extract_amr_info(params, gfa_file):
	"""
	To read AMR info if available or generate them otherwise, including alignment of AMRs in the graph!
	Parameters:
		params: the list of parameters imported from params.py
		gfa_file: the file containing the assembly graph
	Return:
		the list of unique AMR files (the heads of the groups) and their alignment info
		as well as the list of AMRs not found in the graph
	"""
	not_found_amr_names = []
	#if no ref is available detect AMRs in the assembly graph itself
	if params.find_amr_genes:
		logging.info("Finding AMR genes in the assembly graph ...")
		if params.multi_processor:
			unique_amr_files, unique_amr_path_list =\
							find_all_amr_in_graph_parallel(gfa_file, params.output_dir,
							params.amr_db, params.BANDAGE_PATH,
							params.amr_identity_threshold, params.core_num)
		else:
			unique_amr_files, unique_amr_path_list =\
							find_all_amr_in_graph(gfa_file, params.output_dir,
							params.amr_db, params.BANDAGE_PATH,
							params.amr_identity_threshold)
	#if ref genomes are not available but the list of AMRs has already been found in the graph and is available
	else:
		align_dir = params.os.path.join(params.output_dir, AMR_DIR_NAME, AMR_ALIGN_DIR)
		overlap_file_name = os.path.join(params.output_dir, AMR_DIR_NAME, AMR_OVERLAP_FILE)
		#We assume that only unique AMRs (heads) are stored
		unique_amr_files = extract_files(params.amr_files, 'please provide the address of the AMR gene(s)')
		all_align_files = extract_files(align_dir, 'the alignments are not available!')
		heads, member_lists, unique_amr_list = extract_info_from_overlap_file(overlap_file_name)
		if len(unique_amr_list+heads)!=len(unique_amr_files):
			logging.error("inconsisteny between the list of provided unique AMRs and overlapfile info!")
			import pdb;pdb.set_trace()
		amr_count = len(heads) + len(unique_amr_list) + sum(len(list) for list in member_lists)
		#extract path_info for unique AMRs
		unique_amr_path_list = extract_path_info_for_amrs(all_align_files, unique_amr_files,
												amr_count, params.amr_identity_threshold)

	return unique_amr_files, not_found_amr_names, unique_amr_path_list

def full_pipeline_main(params):
	logging.info("Startting the pipeline ...")
	#Validate task values
	task_list = validate_task_values(params.task)
	path_info_files = []
	seq_files = []

	gfa_file = verify_file_existence(params.gfa_file, \
				'please provide the address of the file containing the assembly graph')

	#extract AMR and alignment information
	unique_amr_files, not_found_amr_names, unique_amr_path_list =\
			extract_amr_info(params, gfa_file)
	if not unique_amr_files:
		logging.info("No AMR gene was found!")
		import pdb; pdb.set_trace()
		sys.exit()
	send_amr_align_info = False
	if unique_amr_path_list and len(unique_amr_path_list)==len(unique_amr_files):
		send_amr_align_info = True
	else:
		logging.error("AMR alignment info are not available")
		import pdb; pdb.set_trace()

	#create pairs of seq and align info
	amr_seq_align_info = []
	for i, amr_file in enumerate(unique_amr_files):
		amr_seq_align_info.append((amr_file, unique_amr_path_list[i]))

	if Pipeline_tasks.sequence_neighborhood.value in task_list:
		if MULTIPLE_SEQ_LENGTH:
			for seq_len in [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]:
				params.seq_length = seq_len
				seq_files, path_info_files = sequence_neighborhood_main(params,
						gfa_file, amr_seq_align_info)
		else:
			seq_files, path_info_files = sequence_neighborhood_main(params,
					gfa_file, amr_seq_align_info)

	coverage_annotation_list = []
	if Pipeline_tasks.neighborhood_annotation.value in task_list:
		all_seq_info_lists, annotation_file_list =\
			seq_annotation_main(params, seq_files, path_info_files, unique_amr_files)
		#if MULTIPLE_COV_THR:
		if isinstance(params.coverage_thr, list) and len(coverage_thr)>1:
			coverage_thr_list = params.coverage_thr
			evaluation_dir = os.path.join(params.output_dir, EVAL_DIR, EVAL_DIR+'_'+str(params.seq_length))
			if not os.path.exists(evaluation_dir):
				try:
					os.makedirs(evaluation_dir)
				except OSError as exc:
					if exc.errno != errno.EEXIST:
						raise
					pass
			coverage_evaluation_file = os.path.join(evaluation_dir, 'coverage_evaluation_'+\
				datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.csv')
			with open(coverage_evaluation_file,'a') as fd:
				writer = csv.writer(fd)
				writer.writerow(['coverage_thr', 'precision', 'sensitivity'])
				data = []
				# for cov_thr in range(1, 50):
				# 	params.coverage_thr = cov_thr
				for cov_thr in coverage_thr_list:
					params.coverage_thr = cov_thr
					coverage_annotation_list = seq_annotation_trim_main(params, unique_amr_files,\
					 		all_seq_info_lists, annotation_file_list, False)
			df = pd.DataFrame(data, columns = ['cov_thr', 'value', 'type'])
			sns.scatterplot(data=df, x = 'cov_thr', y = 'value', hue='type', style='type')
			plt.show()
		else:
			if isinstance(params.coverage_thr, list):
				params.coverage_thr = params.coverage_thr[0]
			coverage_annotation_list = seq_annotation_trim_main(params, unique_amr_files,\
				all_seq_info_lists, annotation_file_list, True)

	logging.info("All Done!")

def update_full_pipeline_params(params, config):
	"""
	"""
	config_params = config.keys()
	main_dir_changed = False
	if 'main_dir' in config_params and os.path.realpath(config['main_dir'])!=os.path.realpath(params.main_dir):
		main_dir_changed = True
		#changes params variables dependant on main_dir only accessible through params.py
		params.output_dir = params.output_dir.replace(params.main_dir.rstrip(' /'), config['main_dir'].rstrip(' /'))
		params.amr_files = params.amr_files.replace(params.main_dir.rstrip(' /'), config['main_dir'].rstrip(' /'))
		params.ng_seq_files = params.ng_seq_files.replace(params.main_dir.rstrip(' /'), config['main_dir'].rstrip(' /'))
		params.ng_path_info_files = params.ng_path_info_files.replace(params.main_dir.rstrip(' /'), config['main_dir'].rstrip(' /'))
	if 'task' in config_params:
		params.task = config['task']
	if 'seq_length' in config_params:
		params.seq_length = config['seq_length']
	if 'gfa_file' in config_params:
		params.gfa_file = config['gfa_file']
	elif main_dir_changed:
		params.gfa_file = params.gfa_file.replace(params.main_dir.rstrip(' /'), config['main_dir'].rstrip(' /'))
	if 'find_amr_genes' in config_params:
		params.find_amr_genes = config['find_amr_genes']
	if 'amr_identity_threshold' in config_params:
		params.amr_identity_threshold = config['amr_identity_threshold']
	if 'coverage_thr' in config_params:
		params.coverage_thr = config['coverage_thr']
	if 'multi_processor' in config_params:
		params.multi_processor = config['multi_processor']
	if 'core_num' in config_params:
		params.core_num = config['core_num']
	if 'BANDAGE_PATH' in config_params:
		params.BANDAGE_PATH = config['BANDAGE_PATH']
	if 'PROKKA_COMMAND_PREFIX' in config_params:
		params.PROKKA_COMMAND_PREFIX = config['PROKKA_COMMAND_PREFIX']
	if 'amr_db' in config_params and config['amr_db']!='':
		params.amr_db = config['amr_db']
	if 'max_kmer_size' in config_params:
		params.max_kmer_size = config['max_kmer_size']
	if 'path_node_threshold' in config_params:
		params.path_node_threshold = config['path_node_threshold']
	if 'ng_extraction_time_out' in config_params:
		params.ng_extraction_time_out = config['ng_extraction_time_out']
	if 'use_RGI' in config_params:
		params.use_RGI = config['use_RGI']
	if 'RGI_include_loose' in config_params:
		params.RGI_include_loose = config['RGI_include_loose']
	if main_dir_changed:
		params.main_dir = config['main_dir']

	return params

def create_arguments(params, parser):
    """
    To create all required qrguments
    """
    parser.add_argument('--task', type=int, default=params.task,
            help="which task would you like to do?\
            For the entire pipeline choose "+str(Pipeline_tasks.all.value)+"; otherwise\
            provide a number representing one of the following tasks:\n \
            sequence_neighborhood = "+str(Pipeline_tasks.sequence_neighborhood.value)+\
            "\nneighborhood_annotation = "+str(Pipeline_tasks.neighborhood_annotation.value))
    parser.add_argument('--ref_genome_files', nargs="+", default=params.ref_genome_files,
            help = 'the ddress of reference genomes that AMR genome will be inserted in them')
    parser.add_argument('--main_dir', '-m', type = str, default=params.main_dir,
            help = 'the main dir to retrieve required files')
    parser.add_argument('--read_length',type=int, default=params.read_length,
            help = 'the length of simulated reads can be either 150 or 250')
    parser.add_argument('--spades_thread_num',type=int, default=params.spades_thread_num,
            help = 'the number of threads used for MetaSPAdes')
    parser.add_argument('--seq_length', '-L', type = int, default=params.seq_length,
            help = 'the length of AMR gene\'s neighbourhood to be extracted')
    parser.add_argument('--gfa_file', type = str, default = params.gfa_file,
            help = 'the address of the file for assembly graph')
    parser.add_argument('--contig_file', type = str, default = params.contig_file,
            help = 'the address of the file containing contigs after assembly')
    parser.add_argument('--reads', type = check_reads, default = params.reads,
            help = 'the address of the files containing paired-end reads')
    parser.add_argument('--find_amr_genes', type = str2bool, default = params.find_amr_genes,
            help = 'Whether to assume the AMR genes (in metagenome) are known or to look for them in assembly graph')
    parser.add_argument('--amr_identity_threshold', type = int, default = params.amr_identity_threshold,
            help = 'the threshold used for amr alignment: a hit is returned if identity/coverage >= threshold')
    parser.add_argument('--ref_genomes_available', type = str2bool, default = params.ref_genomes_available,
            help = 'Whether we have access to reference genome(s)')
    parser.add_argument('--multi_processor', type = str2bool, default = params.multi_processor,
            help = 'Whether to use multi processors for parallel programming')
    parser.add_argument('--core_num', type = int, default = params.core_num,
            help = 'the number of cores used in case of parallel programming')
    parser.add_argument('--coverage_thr', nargs="+", default = params.coverage_thr,
            help = 'coverage threshold to check if an annotated gene is truly AMR neighbor or just a false positive')
    # parser.add_argument('--ref_ng_annotations_file', type = str, default = params.ref_ng_annotations_file,
   # parser.add_argument('--ref_ng_annotations_file', type = str, default = params.ref_ng_annotations_file,
    #       help = 'the file containing the annotation of all neighborhoods extracted from ref genomes.')
    # parser.add_argument('--prokka_prefix', type = str, default = params.PROKKA_COMMAND_PREFIX,
    #       help = 'Set only if prokka is run through docker')
    # parser.add_argument('--amr_files','-A', type=str, default = params.amr_files,
    #       help = 'the path of the file(s) containing the AMR gene sequence(s)')
    # parser.add_argument('--output_dir', '-O', type = str, default=params.output_dir,
    #       help = 'the output dir to store the results')
    # parser.add_argument('--metagenome_file', type = str, default= params.metagenome_file,
    #       help = 'the address of metagenome file')
    # parser.add_argument('--assembler_output_dir',type=str, default=params.assembler_output_dir,
    #       help = 'the output dir to store MetaSPAdes results')
    # parser.add_argument('--graph_distance', '-D', type = int, default=params.graph_distance,
    #       help = 'the maximum distance of neighborhood nodes to be extracted from the AMR gene')
    # parser.add_argument('--ng_seq_files', nargs="+", default = params.ng_seq_files,
    #       help = 'the address of the files containing all extracted neighborhood sequences in assembly graph')
    # parser.add_argument('--ng_path_info_files', nargs="+", default = params.ng_path_info_files,
    #       help = 'the address of the files containing all path information for extracted neighborhood sequences in assembly graph')
    # parser.add_argument('--spades_error_correction', type = str2bool, default = params.spades_error_correction,
    #       help = 'Whether to turn on or off error correction in MetaSPAdes')
    # parser.add_argument('--use_RGI', type = str2bool, default = params.use_RGI,
    #       help = 'Whether to contribute RGI annotation in Prokka result')
    # parser.add_argument('--RGI_include_loose', type = str2bool, default = params.RGI_include_loose,
    #       help = 'Whether to include loose cases in RGI result')
    # parser.add_argument('--path_node_threshold', type = int, default = params.path_node_threshold,
    #       help = 'the threshold used for recursive pre_path and post_path search as long as the length of the path is less that this threshold')
    # parser.add_argument('--path_seq_len_percent_threshold', type = int, default = params.path_seq_len_percent_threshold,
    #       help = 'the threshold used for recursive pre_seq and post_seq until we have this percentage of the required length\
    #        after which we just extract from the longest neighbor')
    # parser.add_argument('--number_of_insertions', type = int, default=params.number_of_insertions,
    #       help = 'the number of genomes generated by inserting AMR in different locations of reference genome')
    # parser.add_argument('--insertion_type', type =Insertion_type , default=params.insertion_type,
    #       help = 'Should insertion locations be selected randomly (1) or using some defined values (2)')
    # parser.add_argument('--insertion_locations', nargs="+", default = params.insertion_locations,
    #       help = 'list of predefined insertion locations to insert AMR in reference genome\
    #               if you chose multiple reference genomes first type the list of insertion locations\
    #               for the first genome followed by the list for the second and ...\
    #               please make sure the number of locations for each genome matches the\
    #               corresponding one specified in --num_insertion argument')
    # parser.add_argument('--artificial_amr_insertion', type = str2bool, default = params.artificial_amr_insertion,
    #       help = 'Whether to insert the AMR gene in genomes artificially')
    # parser.add_argument('--genome_amr_files', nargs="+", default = params.genome_amr_files,
    #       help = 'the address of the files containing genome after AMR insertion')

    return parser

if __name__=="__main__":
	import params
	text = 'This code is used to find the context of a given AMR gene'
	parser = argparse.ArgumentParser(description=text)
	parser.add_argument('--config_file', '-C', type = str, default='',
		help = 'the config file to set parameters for full_pipeline()')
	args = parser.parse_args()
    # Read config file into a dictionery
	print("Reading the config file '"+args.config_file+"' ...")
	with open(args.config_file, 'r') as yamlfile:
		data = yaml.load(yamlfile, Loader=yaml.FullLoader)
	params = update_full_pipeline_params(params, data)
	log_name = 'logger_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.log'
	initialize_logger(params.main_dir, log_name)
	validate_print_parameters_tools(params)
	#logging.info(str(params.__dict__))
	#create the output directory; if it exists, delete it and create a new one
	if not os.path.exists(params.output_dir):
		os.makedirs(params.output_dir)
	# else:
	# 	try:
	# 		shutil.rmtree(params.output_dir)
	# 	except OSError as e:
	# 		logging.error("Error: %s - %s." % (e.filename, e.strerror))
	# 	os.makedirs(params.output_dir)
	full_pipeline_main(params)
