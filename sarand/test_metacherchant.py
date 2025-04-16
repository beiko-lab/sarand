
"""
There should be a main directory <main_dir> that contains:
    1- AMR_seqs_full.fasta
    2- AMR_info/sequences folder that contains all the gene sequeces in separate files
    3- output/<amr_name>/graph.gfa for all genes
To run the code make sure to use:
 * -a/--assembler metacherchant
 * --meta_main_dir <main_dir>
"""
import os
import sys
import csv
import pandas as pd
import datetime
import gfapy
import subprocess
from pathlib import Path

from sarand.util.logger import LOG
from sarand.external.blastn import Blastn
from sarand.utils import (
    restricted_amr_name_from_modified_name,
    amr_name_from_comment,
    read_path_info_from_align_file,
    create_fasta_file,
    retrieve_AMR,
    extract_name_from_file_name,
)
#from sarand.extract_neighborhood import neighborhood_sequence_extraction
#from sarand.extract_neighborhood import sequence_neighborhood_main
from sarand.new_extract_neighborhood import sequence_neighborhood_main
from sarand.full_pipeline import (
    neighborhood_annotation,
    check_coverage_consistency_remove_rest_seq,
    find_corrsponding_seq_path_file,
)
from sarand.config import SEQ_NAME_PREFIX, ANNOTATION_DIR, EVAL_DIR, AMR_DIR_NAME, \
    AMR_SEQ_DIR, SEQ_DIR_NAME, CONDA_BANDAGE_NAME, CONDA_EXE_NAME, AMR_ALIGN_DIR, \
    CONDA_BLAST_NAME
from sarand.external.bandage import Bandage, BandageParams
RUN_BANDAGE = False

def extract_files(gfiles, message):
    """
    To extract file(s) address from an object
    # if gfiles is a list and the first item is a file address (it would be more
    # accurate to check this for all items) return gfiles
    # else if gfiles is a file address return [gfiles] as a list with one item
    # else if gfiles is a directory address return the list of files in it
    Parameters:
        gfiles:		a string or list of strings
        message:	an error message in case that no file was extracted successfully
    Return:
        the list of file(s) address
    """
    if isinstance(gfiles, list):
        # check if these are files (not directories)
        if os.path.isfile(gfiles[0]):
            return gfiles
        else:
            LOG.error(message)
            sys.exit(1)
        # elif os.path.isdir(gfiles[0])
    elif os.path.isfile(gfiles):
        return [gfiles]
    elif os.path.isdir(gfiles):
        myfiles = [
            os.path.join(gfiles, f)
            for f in os.listdir(gfiles)
            if os.path.isfile(os.path.join(gfiles, f))
        ]
        return myfiles
    elif message != "":
        LOG.error(message)
        sys.exit(1)

def find_annotation_file(annotation_file, out_dir,amr_name, seq_len):
    """
    """
    if annotation_file!='':
        return annotation_file
    annotate_dir = os.path.join(out_dir, ANNOTATION_DIR, ANNOTATION_DIR+'_'+str(seq_len),
                                    'annotation_'+amr_name+'_'+str(seq_len))
    for f in os.listdir(annotate_dir):
        if f=='trimmed_annotation_info_'+amr_name+'.csv':
            return os.path.join(annotate_dir, f)
    return ''

def find_sequence_file(seq_file, sequence_dir, amr_name, seq_len):
    """
    """
    #seq_file has already been found
    if seq_file!='':
        return seq_file
    seq_dir = os.path.join(sequence_dir, 'sequences')
    for f in os.listdir(seq_dir):
        if f.startswith(SEQ_NAME_PREFIX+amr_name+'_'+str(seq_len)) and\
            os.path.isfile(os.path.join(seq_dir, f)) and\
            os.path.splitext(f)[1]=='.txt':
            return os.path.join(seq_dir, f)
    return ''

def extract_amr_info_from_file(amr_info_file):
    """
    """
    amr_names = []
    amr_seqs = []
    with open(amr_info_file) as fd:
        for line in fd:
            if line.startswith('>'):
                amr_comment = line[1:-1]
                continue
            amr_name = amr_name_from_comment(amr_comment)
            amr_names.append(amr_name)
            amr_seqs.append(line[:-1])
    return amr_names, amr_seqs

def find_corresponding_amr_file(amr_name, amr_files):
    """
    """
    for amr_file in amr_files:
        if os.path.splitext(os.path.basename(amr_file))[0]==amr_name:
            return amr_file
    LOG.error("Error: No amr_file for amr "+amr_name+" was found!")
    sys.exit()

def find_amr_related_nodes(amr_file, gfa_file, output_dir,
							threshold =  95, output_pre = ''):
	"""
	Run bandage+blast to find the sequence in amr_file (as the query) in the assembly
	graph (gfa_file), and extract the path(s) representing the AMR sequence (i.e.,
	AMR path) in the assembly graph.
	e.g., (121) A+, B-, C+ (143) --> paths_info =	[
										{nodes:['A', 'B', 'C'],
										orientations: ['+', '-', '+'],
										start_pos: 121,
										end_pos: 143}
													]
	Parameters:
		amr_file: the FASTA file containing the AMR sequence (our query)
		gfa_file: the GFA file containing the assembly graph (our DB)
		output_dir: the output directory to store files
		bandage_path: the address of bandage executation file
		threshold: threshold for identity and coverage
		output_pre: used for naming output file
	Return:
		A boolean value which is True if any path was found and
		A list of dictionaries each denoting an AMR path
	"""
	#Run bandage+blast
	output_name=os.path.join(output_dir, output_pre+'_align_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
	if os.path.isfile(output_name+'.tsv'):
		os.remove(output_name+'.tsv')
	bandage_command = subprocess.run([CONDA_EXE_NAME, 'run', '-n', CONDA_BANDAGE_NAME, 'Bandage',
                    "querypaths", gfa_file, amr_file,
					output_name, "--pathnodes", "50", "--minpatcov",
					str((threshold-1)/100.0), "--minmeanid", str((threshold-1)/100.0),
					"--minhitcov", str((threshold-1)/100.0)],
					stdout=subprocess.PIPE, check= True )
	LOG.info(bandage_command.stdout.decode('utf-8'))
	# command = bandage_path +' querypaths '+gfa_file+' '+amr_file+' '+output_name + ' --pathnodes 50'
	# os.system(command)
	#Process the output tsv file
	found, paths_info = read_path_info_from_align_file(output_name+".tsv", threshold)
	if (found):
		return found, [paths_info[0]]
	else:
		return found, []

def order_path_nodes(path_nodes, amr_file, out_dir, threshold = 90):
	"""
	Given that we have a list of nodes that are supposed to represent a given AMR
	(i.e., AMR path), this method returns their right order to represent the AMR
	sequence.
	Curretly, this method is only used for Metacherchant
	Parameters:
		path_nodes: the list of nodes representing AMR
		amr_file: the file containing AMR sequence
		out_dir: the dir to store some temporary files
		threshold: used for identity and coverage in bandage+blast
	Returns:
		the lists of sorted nodes and their corresponding orientation
	"""
	path_nodes_info = []
	no_path_nodes = []
	for i, node in enumerate(path_nodes):
		node_info_list = []
		#write the sequence into a fasta file
		query_file = create_fasta_file(node.sequence, out_dir, comment = ">"+node.name+"\n", file_name = 'query')
		#run blast query for alignement
		blast_file_name = os.path.join(out_dir,'blast.csv')
		blast_file = open(blast_file_name, "w")
		#cmd = [CONDA_EXE_NAME,'run', '-n', CONDA_BLAST_NAME, 'blastn',
		cmd = ['blastn',
                        '-query', query_file, "-subject", amr_file,
			"-task", "blastn", "-outfmt", "10", "-max_target_seqs", "10",
			"-evalue", "0.5", "-perc_identity", str(threshold-1)]
		LOG.debug(' '.join(map(str, cmd)))
		proc = subprocess.Popen(
			cmd,
			#stdout=subprocess.PIPE,
			stdout=blast_file,
			stderr=subprocess.PIPE,
			encoding='utf-8'
		)
		stdout, stderr = proc.communicate()
		if proc.returncode != 0:
			raise RuntimeError(f'blastn failed: {stderr}')
		#blast_command = subprocess.run([CONDA_EXE_NAME,'run', '-n',CONDA_BLAST_NAME,
                #        "-query", query_file, "-subject", amr_file,
		#				"-task", "blastn", "-outfmt", "10", "-max_target_seqs", "10",
		#				"-evalue", "0.5", "-perc_identity", str(threshold-1)],
		#				stdout=blast_file, check= True)
		blast_file.close()
		# command = 'blastn -query '+query_file+' -subject '+amr_file+\
		# 	' -task blastn -outfmt 10 -max_target_seqs 10 -evalue 0.5 -perc_identity '+\
		# 	str(threshold)+' > '+ blast_file_name
		# os.system(command)
		with open(blast_file_name, 'r') as file1:
			myfile = csv.reader(file1)
			for row in myfile:
				identity=int(float(row[2]))
				coverage = int(float(row[3])/len(node.sequence)*100)
				if identity>=threshold and coverage>=threshold:
					node_info = {'name':node.name, 'c_start':int(row[8]), 'c_end':int(row[9])}
					node_info_list.append(node_info)
		if not node_info_list:
			no_path_nodes.append(i)
			LOG.error(node.name+' was not found in the sequence of '+amr_file)
		path_nodes_info.append(node_info_list)
	#order nodes
	start_list = []
	orientations = []
	for path_info in path_nodes_info:
		if len(path_info)>0:
			if path_info[0]['c_start'] < path_info[0]['c_end']:
				start_list.append(path_info[0]['c_start'])
				orientations.append('+')
			else:
				start_list.append(path_info[0]['c_end'])
				orientations.append('-')
		else:
			start_list.append(-1)
			orientations.append('/')
	# start_list =[e[0]['c_start'] if len(e)>0 else -1 for e in path_nodes_info ]
	sorted_indeces = sorted(range(len(start_list)), key=lambda k: start_list[k])
	sorted_path_nodes = []
	sorted_orientations = []
	for index in sorted_indeces:
		if index not in no_path_nodes:
			sorted_path_nodes.append(path_nodes[index])
			sorted_orientations.append(orientations[index])

	return sorted_path_nodes, sorted_orientations

def extract_amr_align_from_file(gfa_file):
	"""
	Retrieve the list of segments that represent AMR path
	This method is use for Metacherchant in which such nodes ends with '_start'
	Parameters:
	 	gfa_file: the assembly graph file
	Return:
		the list of graph segments representing the AMR path
	"""
	myGraph = gfapy.Gfa.from_file(gfa_file)
	path_nodes = []
	for segment in myGraph.segments:
		if segment.name.endswith('_start'):
			path_nodes.append(segment)
	return path_nodes

def gene_alignment_extraction_metacherchant(gfa_file, output_dir,
									threshold = 90, amr_file = ''):
	"""
	to extract the AMR alignment to the local graph for metacherchant
	Parameters:
		gfa_file:	the GFA file containing the assembly graph
		length:		the length of all sequences around the AMR gene to be extracted
		output_dir: the output directory to store files
		threshold: 	threshold for identity and coverage
		amr_file:	the FASTA file containing the AMR sequence
	Return:
		the name of file containing the list of extracted sequences/paths
	"""
	with open(os.path.join(output_dir, "metacherchant_no_path.txt"), 'a') as no_path_file:
		no_path_file.write(amr_file+'\n')
	# try:
	# 	myGraph = gfapy.Gfa.from_file(gfa_file)
	# except Exception as e:
	# 	LOG.error(e)
	# 	return ''
	#find nodes ending at _start as the nodes that create amr_path
	path_nodes = extract_amr_align_from_file(gfa_file)
	#find the order of these nodes in the amr path
	ordered_path_nodes, orientations = order_path_nodes(path_nodes, amr_file,
									os.path.join(output_dir,AMR_ALIGN_DIR), threshold)
	node_list = [e.name for e in ordered_path_nodes]

	last_segment = ordered_path_nodes[-1]
	start_pos = 1
	end_pos = len(str(last_segment.sequence))

	path_info = {
        	"nodes": node_list,
        	"orientations": orientations,
        	"start_pos": start_pos,
        	"end_pos": end_pos,
    	}
	LOG.debug('last_segment = '+last_segment.name+' start_pos = '+str(start_pos)+' end_pos= '+str(end_pos))

	return [path_info]

def test_metacherchant_main(params):
    """
    Main runner function for testing metacherchant
    """
    LOG.info("Starting testing Metacherchant...")

    amr_info_file = os.path.join(params.meta_main_dir, "AMR_seqs_full.fasta")
    amr_names, amr_seqs = extract_amr_info_from_file(amr_info_file)
    amr_files = os.path.join(params.meta_main_dir ,AMR_DIR_NAME, AMR_SEQ_DIR)
    ref_amr_files = extract_files(amr_files, 'please provide the address of the AMR gene(s)')
    #if params.ref_genomes_available:
    #    df = pd.read_csv(params.ref_ng_annotations_file, skipinitialspace=True,  keep_default_na=False)
    #    amr_groups = df.groupby('target_amr')
    sequence_dir = os.path.join(params.output_dir, SEQ_DIR_NAME)
    if not os.path.exists(sequence_dir):
        os.makedirs(sequence_dir)
    ##to make use of previously generated alignment files
    #align_files = []
    #align_dir = os.path.join(sequence_dir, 'alignment_files')
    #if os.path.exists(align_dir):
    #    align_files = [os.path.join(align_dir, f) for f in os.listdir(align_dir) \
	#			if os.path.isfile(os.path.join(align_dir, f)) and os.path.splitext(f)[1]=='.tsv']
    #evaluation_dir = os.path.join(params.output_dir, EVAL_DIR, EVAL_DIR+'_'+str(params.seq_length))
    #if not os.path.exists(evaluation_dir):
    #    os.makedirs(evaluation_dir)
    #summary_file = os.path.join(evaluation_dir, 'summaryMetrics_up_down_metacherchant_'+\
    #    datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.csv')
    #with open(summary_file,'a') as fd:
    #    writer = csv.writer(fd)
    #    writer.writerow(['AMR', 'Unique_TP#', 'FP#', 'Unique_True#', 'found#','sensitivity', 'precision', 'not_found'])

    #average_precision = 0
    #average_sensitivity = 0
    #no_align_file = open(os.path.join(params.output_dir, 'no_align_found.txt'), 'w')
    align_dir = os.path.join(params.output_dir, AMR_ALIGN_DIR)
    if not os.path.exists(align_dir):
        os.makedirs(align_dir)
    amr_seq_align_info = []
    seq_files = []
    path_info_files = []
    LOG.info('amr_names: ' + str(amr_names))
    for amr_name in amr_names:
        LOG.info('Processing '+amr_name+' ...')
        restricted_amr_name = restricted_amr_name_from_modified_name(amr_name)
        #amr_name1 = amr_name.replace(";",'SS')
        #restricted_amr_name1 = ''.join(e for e in amr_name1 if e.isalpha() or e.isnumeric() or e=='_' or e=='-')
        seq_file = ''
        # sequence extraction
        #if Pipeline_tasks.sequence_neighborhood.value in task_list:
        amr_file = find_corresponding_amr_file(restricted_amr_name, ref_amr_files)
        #found = False
        #for myfile in align_files:
        #    if os.path.basename(myfile).startswith('ng_sequences_'+restricted_amr_name+'_align') or\
        #        os.path.basename(myfile).startswith('ng_sequences_'+restricted_amr_name+restricted_amr_name+'_align'):
        #            align_file = myfile
        #            found, paths_info = read_path_info_from_align_file(align_file, 95)
        #            break
        #    if found:
        #        amr_pair = (amr_file,paths_info)
        #    else:
        #gfa_file = params.output_dir+'output/'+str(i+1)+'/graph.gfa'
        gfa_file = os.path.join(params.meta_main_dir, 'output', restricted_amr_name, 'graph.gfa')
        command = "sed -e 's/\tCL:z:GREEN//g' -i "+gfa_file
        print(command)
        os.system(command)
        LOG.info('gfa_file: '+gfa_file)
        output_dir = params.output_dir
        found = False
        if RUN_BANDAGE:
                LOG.debug('calling Bandage to find ' + amr_name)
                found, amr_paths_info = find_amr_related_nodes(amr_file, gfa_file,
						align_dir,
						params.min_target_identity, restricted_amr_name)
                LOG.info("found: " + str(found))
                LOG.info("amr_paths_info:" + str(amr_paths_info))
        if not found:
            LOG.debug("Calling gene_alignment_extraction_metacherchant for "+amr_name)
            amr_paths_info = gene_alignment_extraction_metacherchant(gfa_file, params.output_dir,
	        params.min_target_identity, amr_file)

        #amr_pair = (amr_file, amr_paths_info)
        amr_seq_align_info.append((amr_file, amr_paths_info))
        amr_seq_align_info_fake = []
        amr_seq_align_info_fake.append((amr_file, amr_paths_info))
        LOG.debug("amr_paths_info:" + str(amr_paths_info))
        LOG.info('Neighborhood Extraction for ' + amr_file + '...')
        seq_files_fake, path_info_files_fake = sequence_neighborhood_main(
            params,
            Path(gfa_file),
            amr_seq_align_info_fake,
            params.debug
        )
        seq_files.extend(seq_files_fake)
        path_info_files.extend(path_info_files_fake)
    #no_align_file.close()
    #LOG.info('Neighborhood Extraction ...')
    #seq_files, path_info_files = sequence_neighborhood_main(
    #    params,
    #    Path(gfa_file),
    #    amr_seq_align_info,
    #    params.debug
    #)

    #Annotation
    LOG.info('Neighborhood Annoation ...')
    if seq_files:
        neighborhood_files = seq_files
    else:
        LOG.error("No file containing the extracted neighborhood sequences is available!")
        import pdb
        pdb.set_trace()

    if path_info_files:
        nodes_info_files = path_info_files
    else:
        LOG.error("No file containing path info for neighborhood sequences is available!")
        import pdb
        pdb.set_trace()
    all_seq_info_lists = []
    annotation_files = []
    for amr_name in amr_names:
        LOG.info('Annnnnotation for '+amr_name+' ...')
        restricted_amr_name = restricted_amr_name_from_modified_name(amr_name)
        amr_file = find_corresponding_amr_file(restricted_amr_name, ref_amr_files)
        #restricted_amr_name = extract_name_from_file_name(amr_file)
        #_, amr_name = retrieve_AMR(amr_file)
        neighborhood_file, nodes_info_file = find_corrsponding_seq_path_file(
            restricted_amr_name,
            neighborhood_files,
            nodes_info_files,
            params.neighbourhood_length,
        )
        if neighborhood_file == -1:
            LOG.error(
                "no sequence file for "
                + amr_file
                + " was found! We looked for a file like "
                + restricted_amr_name
            )
            import pdb
            pdb.set_trace()
            sys.exit(1)
        all_seq_info_list, annotation_file = neighborhood_annotation(
            amr_name,
            neighborhood_file,
            nodes_info_file,
            params.neighbourhood_length,
            params.output_dir,
            params.no_rgi,
            params.rgi_include_loose,
            "_" + restricted_amr_name,
            params.num_cores,
        )
        all_seq_info_lists.append(all_seq_info_list)
        annotation_files.append(annotation_file)

if __name__=="__main__":
    #text = 'This code is used to test metacherchant'
    #log_name = 'logger_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.log'
    #initialize_logger(params.main_dir, log_name)
    # parser = argparse.ArgumentParser(description=text)
    # parser.add_argument('--amr', type=str, default='',
    #     help='the path of the directory containing the amr files')
    # parser.add_argument('--seq', type=str, default = '',
    #     help = 'the path of the fasta file containing all AMR sequences')
    # args = parser.parse_args()
    main()
