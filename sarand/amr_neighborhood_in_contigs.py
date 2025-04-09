
"""
File:		amr_neighborhood_in_contigs.py
Aythor:		Somayeh Kafaie
Date:		March 2021
Purpose:	To find the neighborhood of AMRs in a contig file, compare them with
			that of the ref genomes and calculate the sentivity and precision

To run:
$ conda activate rgi
$ python amr_neighborhood_in_contigs.py
NOTE: It reads required parameters from params.py and the most important parameters need
to be set correctly there are:
params.seq_length, params.contig_file, params.amr_identity_threshold, params.amr_files,
params.ref_ng_annotations_file, params.main_dir, params.output_dir,
params.PROKKA_COMMAND_PREFIX, params.use_RGI,params.RGI_include_loose,
params.ref_genomes_available
NOTE: The result are available in the following directory:
params.output_dir+'contigs_output_'+str(params.seq_length)

"""

################################################################################

import os
import re
import glob
import argparse
import datetime
import csv
import collections
import subprocess
from csv import DictReader
import pandas as pd
from Bio import SeqIO
from gfapy.sequence import rc

from sarand.util.logger import LOG
from sarand.config import AMR_SEQ_DIR
from sarand.utils import amr_name_from_comment, split_up_down_info,\
		annotate_sequence,restricted_amr_name_from_modified_name,\
		create_fasta_file

NOT_FOUND_FILE = 'not_found_amrs_in_contigs.txt'


def annotate_sequence_bundle(contig_ng_info, out_dir, no_RGI,
								RGI_include_loose):
	"""
	"""
	annotate_dir = os.path.join(out_dir,'contig_annotations')
	if not os.path.exists(annotate_dir):
		os.makedirs(annotate_dir)
	for amr_info_list in contig_ng_info:
		amr_name = amr_info_list[0]
		restricted_amr_name = restricted_amr_name_from_modified_name(amr_name)
		# mydir = annotate_dir+restricted_amr_name
		# if not os.path.exists(mydir):
		# 	os.makedirs(mydir)
		annotation_file_name =os.path.join(annotate_dir, 'contig_annotation_'+restricted_amr_name+'.csv')
		with open(annotation_file_name, 'a') as fd:
			writer = csv.writer(fd)
			writer.writerow(['seq_name', 'seq_value', 'seq_length', 'gene', 
							'product', 'length', 'start_pos', 'end_pos', 'RGI_prediction_type',
							 'family', 'target_amr'])
		for index, amr_info in enumerate(amr_info_list[1]):
			contig_name = amr_info['contig']
			seq = amr_info['seq']
			seq_description = 'contig_'+amr_name+'_'+contig_name.replace(' ','_').replace('.','').replace(',','')
			annotation_prefix = 'contig_'+ restricted_amr_name +'__'+str(index)
			seq_info = annotate_sequence(seq+"\n", annotation_prefix, annotate_dir,
                    no_RGI, RGI_include_loose)
			found, _ , _, _, _ = split_up_down_info(seq, seq_info)
			if not found:
				LOG.error("no target amr was found in this contig sequence: "+contig_name)
			with open(annotation_file_name, 'a') as fd:
				writer = csv.writer(fd)
				for gene_info in seq_info:
					writer.writerow([seq_description, gene_info['seq_value'],
										len(gene_info['seq_value']),
										gene_info['gene'],
										gene_info['product'], gene_info['length'],
										gene_info['start_pos'], gene_info['end_pos'],
										gene_info['RGI_prediction_type'],
										gene_info['family'], gene_info['target_amr']])
		LOG.info("NOTE: The annotation of neighborhood sequences in contigs for "+\
			amr_name+"has been stroed in " + annotation_file_name)

def extract_amr_length(target_genes_file):
    """
    """
    amr_objects = []
    with open(target_genes_file) as fp:
        for line in fp:
            if line.startswith('>'):
                amr_comment = line[1:-1]
                continue
            amr_name = amr_name_from_comment(amr_comment)
            amr_object={'name':amr_name, 'length':len(line)-1, 'seq': line[:-1], 'title':amr_comment}
            amr_objects.append(amr_object)
    return amr_objects

def find_all_amrs_and_neighborhood(target_genes_file, genome_file, out_dir,
									neighborhood_len = 1000, threshold = 95):
	"""
	"""
	if genome_file == '':
		LOG.error('Please enter the address of contig file!')
		sys.exit()
	blast_file_name = os.path.join(out_dir, 'blast_out_contig.csv')
	if not os.path.isfile(blast_file_name):
		# Find the length of each AMR sequence
		amr_objects = extract_amr_length(target_genes_file)
		#creat blast database from the (meta)genome file
		db_command = subprocess.run(["makeblastdb","-in", genome_file, "-parse_seqids",
								"-dbtype", "nucl"], stdout=subprocess.PIPE, check= True)
		LOG.info(db_command.stdout.decode('utf-8'))
		blast_file = open(blast_file_name, "w")
		blast_command = subprocess.run(["blastn", "-query", target_genes_file, "-db", genome_file,
						"-outfmt", "10", "-evalue", "0.5", "-perc_identity", str(threshold-1),
						"-num_threads", "4"], stdout=blast_file, check= True)
		blast_file.close()
	AMR_dir = os.path.join(out_dir, AMR_SEQ_DIR)
	ng_file = os.path.join(out_dir, 'AMR_contig_neighborhood.fasta')
	if not os.path.exists(AMR_dir) or not os.path.isfile(ng_file):
		#Read the blast result
		amr_list = []
		amr = collections.namedtuple('amr', 'amr_name seq_name identity matched_length q_start q_end c_start c_end')
		with open(blast_file_name, 'r') as file1:
			myfile = csv.reader(file1)
			for row in myfile:
				myamr = amr(amr_name=row[0], seq_name=row[1], identity=row[2],
                    matched_length=row[3], q_start=row[6], q_end=row[7],
                    c_start=row[8], c_end=row[9])
				amr_list.append(myamr)
        #Find the list of detected AMRs
		amr_start_list = []
		amr_end_list = []
		record_name_list = []
		amr_name_list = []
		if not os.path.exists(AMR_dir):
			os.makedirs(AMR_dir)
        # find the start and end of AMR for all found cases above the threshold
		for record in amr_list:
			if int(float(record.identity))>=threshold:
				target_amr = next((amr_obj for amr_obj in amr_objects if amr_obj['title'].split(' ')[0] == record.amr_name), None)
				if target_amr and int(float(record.matched_length)/target_amr['length']*100)>=threshold:
					if target_amr['name'] not in amr_name_list:
						amr_file_name = restricted_amr_name_from_modified_name(target_amr['name'])
						amr_file = create_fasta_file(target_amr['seq']+'\n', AMR_dir, '>'+target_amr['title']+'\n', amr_file_name)
					amr_name_list.append(target_amr['name'])
					amr_start_list.append(int(record.c_start))
					amr_end_list.append(int(record.c_end))
					record_name_list.append(record.seq_name)
				elif not target_amr:
					LOG.error("ERROR: couldn't find the length of AMR: "+str(record.amr_name))
					#import pdb; pdb.set_trace()
					sys.exit()
        # extract neighborhood sequence(s)
		ng_lists = []
		amr_list = []
		with open(ng_file, 'w') as myfile:
			seq_list = []
			contig_name_list = []
			for i, amr_start in enumerate(amr_start_list):
				amr_end = amr_end_list[i]
				record_name = record_name_list[i]
				seq = ''
				reverse_complement = False
                # extract sequence from both sides
				if amr_start > amr_end:
					amr_start, amr_end = amr_end, amr_start
					reverse_complement = True
				for record in SeqIO.parse(open(genome_file,'r'),'fasta'):
                    #if record.id == record_name:
					if record_name == record.id:
						amr_start = amr_start -1
						if amr_start-neighborhood_len >= 0:
							seq = str(record.seq[amr_start-neighborhood_len:amr_start]).upper()+str(record.seq[amr_start:amr_end]).lower()
							amr_start = neighborhood_len
						else:
							seq = str(record.seq[:amr_start]).upper()+str(record.seq[amr_start:amr_end]).lower()
						if (amr_end-1+neighborhood_len)<=len(record):
							seq+=str(record.seq[amr_end:amr_end+neighborhood_len])
						else:
							seq+=str(record.seq[amr_end:])
						if reverse_complement:
							seq = rc(seq)
						ng_item ={'contig':record.description, 'seq':seq}
						if amr_name_list[i] not in amr_list:
							amr_list.append(amr_name_list[i])
							ng_lists.append([ng_item])
						else:
							amr_index = amr_list.index(amr_name_list[i])
							ng_lists[amr_index].append(ng_item)
						myfile.write('>'+amr_name_list[i]+'::'+record.description+'\n')
						myfile.write(seq+'\n')
						break
	else:
		amr_list = []
		ng_lists = []
		with open(ng_file, 'r') as myfile:
			for line in myfile:
				if line.startswith('>'):
					items = line[1:-1].split('::')
					amr_name = items[0]
					contig_name = items[1]
				else:
					ng_item ={'contig':contig_name, 'seq':line[:-1]}
					if amr_name not in amr_list:
						amr_list.append(amr_name)
						ng_lists.append([ng_item])
					else:
						amr_index = amr_list.index(amr_name)
						ng_lists[amr_index].append(ng_item)

	#delete tempoarary blastdb created for genome_file
	directory = os.path.dirname(genome_file)
	base = os.path.basename(genome_file)
	pattern = os.path.join(directory, f"{base}.*")
	for filepath in glob.glob(pattern):
    		if filepath != genome_file and os.path.isfile(filepath):
        		os.remove(filepath)
        		
	return zip(amr_list, ng_lists)

def find_contig_amrs_main(params):
	"""
	"""
	#creating a directory for results
	contig_dir = os.path.join(params.output_dir, 'contigs_output_'+str(params.neighbourhood_length))
	if not os.path.exists(contig_dir):
		os.makedirs(contig_dir)
	found_file = os.path.join(contig_dir , NOT_FOUND_FILE)

	#read card sequences and do neighborhood extraction and annotations for blast hits
	# here we assume params.input_gfa is representing the contig file from which
	# the neighborhoods will be extracted
	contig_ng_info = find_all_amrs_and_neighborhood(params.target_genes, params.input_gfa,
								contig_dir, params.neighbourhood_length,
								params.min_target_identity)
	annotate_sequence_bundle(contig_ng_info, contig_dir,
						params.no_rgi, params.rgi_include_loose,)
	LOG.info("neighborhood extraction and annotation from cotigs is done!")

if __name__=="__main__":
	main()
