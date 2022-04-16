"""
File:		params.py
Aythor:		Somayeh Kafaie
Date:		September 2020

Purpose:	To set the parameters used in the project

NOTE: all parameters can be set in params.py
NOTE: if use_RGI = TRUE, make sure either RGI has been installed system-wide or
	you already are in the environment RGI installed in!
"""
import enum
import os
import pkg_resources

"""
# The list of valid tasks
# to use them set task in params.py:
# for a single task: insert its value as the single item in task list
# for a range of tasks: insert the value of start and end as the first and the second (last) item in the list
# for all tasks: task = [0]
"""
class Pipeline_tasks(enum.Enum):
	all = 0
	sequence_neighborhood = 1
	neighborhood_annotation = 2

class Assembler_name(enum.Enum):
    metaspades = 1
    megahit = 2
    bcalm = 3
    metacherchant = 4
    spacegraphcats = 5

BANDAGE_PATH = 'Bandage'
#BANDAGE_PATH = '/media/Data/tools/Bandage_Ubuntu_dynamic_v0_8_1/Bandage'
#cwd = os.getcwd()
#PROKKA_COMMAND_PREFIX = 'docker run -v '+cwd+':/data staphb/prokka:latest '
PROKKA_COMMAND_PREFIX = ''

amr_db = pkg_resources.resource_filename(__name__, 'data/CARD_AMR_seq.fasta')
main_dir = 'test'
output_dir = os.path.join(main_dir , 'output_dir')

multi_processor = True
core_num = 4
#-1 when no gene-coverage threshold is used
coverage_thr = 30
task = [0]

amr_files = os.path.join(output_dir ,'AMR_info/sequences/')
find_amr_genes = True
amr_identity_threshold =  95

assembler = Assembler_name.metaspades
#Setting for assembler
if assembler == Assembler_name.metaspades:
	gfa_file = os.path.join(os.path.join(main_dir , 'spade_output') , 'assembly_graph_with_scaffolds.gfa')
	#contig_file = os.path.join(assembler_output_dir , 'contigs.fasta')
	max_kmer_size = 55
elif assembler == Assembler_name.megahit:
	gfa_file = os.path.join(os.path.join(main_dir , 'megahit_output') , 'k59.gfa')
	#contig_file = os.path.join(assembler_output_dir , 'intermediate_contigs/k59.contigs.fa')
	max_kmer_size = 59
	#max_kmer_size = 119
elif assembler == Assembler_name.bcalm:
	gfa_file = os.path.join(os.path.join(main_dir , 'bcalm_output') , 'bcalm_graph_55.gfa')
	#doesn't produce any contig file; so just sue the one from meta-pades for evaluation
	#contig_file = os.path.join(main_dir , 'spade_output/contigs.fasta')
	max_kmer_size = 54
elif assembler == Assembler_name.metacherchant:
	max_kmer_size = 30
elif assembler == Assembler_name.spacegraphcats:
	gfa_file = 'bcalm_graph.gfa'
	max_kmer_size = 30

seq_length = 1000

# to be used in sequence extraction and after having path_node_threshold number of
#nodes in our path or already extracting path_seq_len_percent_threshold percent of
#neighborhood sequence, we stop or only traverse the path from the longest neighbor of the current node
path_node_threshold =  1000
#time-out in seconds to stop pre_seq or post_seq extraction
#default is -1, meaning that there is no time out
ng_extraction_time_out = -1

use_RGI =  True
RGI_include_loose = False

ng_seq_files = os.path.join(output_dir , 'sequences_info/sequences_info_'+str(seq_length)+'/sequences/')
ng_path_info_files = os.path.join(output_dir , 'sequences_info/sequences_info_'+str(seq_length)+'/paths_info/')
