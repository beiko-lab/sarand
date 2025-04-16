"""
File:		extract_neighborhood.py
Author:		Somayeh Kafaie
Date:		July 2020
Purpose:	To extract the neighborhood of an AMR gene from an assembly graph

To run:
- The Bandage+BLAST implementation:
	python extract_neighborhood.py --amr/-A <AMR gene file path in FASTA format>
	--gfa/-G <GFA assembly graph>
	--length/-L <length of the linear sequence around AMR gene to be extracted (default = 1000)>

"""

################################################################################

import sys
import os
import tempfile

import gfapy
import datetime
import threading
import time

from functools import partial
from multiprocessing.pool import Pool
from sarand.util.file import try_dump_to_disk
import csv
import subprocess
from gfapy.sequence import rc
import shutil
import multiprocessing
from csv import DictReader
from pathlib import Path
from sarand.util.logger import LOG


from sarand.config import AMR_DIR_NAME, AMR_SEQ_DIR, AMR_ALIGN_DIR, AMR_OVERLAP_FILE, SEQ_DIR_NAME, SEQ_NAME_PREFIX, \
    ANNOTATION_DIR

from sarand.utils import (
    reverse_sign,
    find_node_name,
    exist_in_path,
    compare_two_sequences,
    read_path_info_from_align_file,
    retrieve_AMR,
    create_fasta_file,
    split_up_down_info,
    seqs_annotation_are_identical,
    similar_seq_annotation_already_exist,
    amr_name_from_comment,
    extract_name_from_file_name,
    restricted_amr_name_from_modified_name,
    delete_lines_started_with,
    read_path_info_from_align_file_with_multiple_amrs,
    annotate_sequence,
    extract_amr_sequences,
)

OUT_DIR = "output"
TEMP_DIR = "temp"


import networkx as nx



import gfapy
import gzip

import os
import subprocess
import shutil
import multiprocessing



def calculate_coverage(node, max_kmer_size, node_name, assembler):
    """
    To calculate node coverage based on assembler type
    Parameters:
            segment/node: the node whose covearge is supposed to eb calculated
            max_kmer_size: the max_kmer_size used by assembler
    Return:
            coverage

    """
    if assembler == "metaspades" or assembler == "bcalm":
        coverage = node['KC'] / (len(node['sequence']) - max_kmer_size)
    elif assembler == "megahit":
        coverage = float(node_name.split("cov_")[1].split("_")[0])
    elif assembler == "metacherchant" or assembler == "spacegraphcats":
        coverage = -1
    else:
        LOG.error(
            "no way of calculating node coverage has been defined for this assembler!"
        )
        import pdb

        pdb.set_trace()
        sys.exit()
    return coverage

def get_sequnce_path(directed_graph, path, threshold, len_amr, up_down):
    sequence = ""
    if (up_down == "down"):
        sequence = directed_graph.nodes[path[0]]['sequence'][-len_amr:]
        for i in range(1, len(path)):
            #print(i)
            sequence = sequence + \
                directed_graph.nodes[path[i]
                                     ]['sequence'][directed_graph[path[i-1]][path[i]]['overlap']:]
        sequence = sequence[0:threshold]
    elif (up_down == "up"):
        for i in range(len(path)-1, 0, -1):
            sequence = sequence + directed_graph.nodes[path[i]]['sequence'][:-(directed_graph[path[i-1]][path[i]]['overlap'])]
        sequence = sequence + \
            directed_graph.nodes[path[0]]['sequence'][0:len_amr]
        sequence = sequence[-threshold:]
    return sequence


def reverse_complement(sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_sequence = sequence[::-1]
    reverse_complement_sequence = ''.join(
        [complement[base] for base in reverse_sequence])
    return reverse_complement_sequence


def create_directed_graph_nx(gfa_graph):
    directed_graph = nx.DiGraph()

    # create nodes
    for segment in gfa_graph.segments:
        directed_graph.add_node(
            segment.name + "+", name=segment.name + "+", sequence=segment.sequence, KC=segment.KC)
        # when we have - , that the reverse complement of the sequence!
        directed_graph.add_node(
            segment.name + "-", name=segment.name + "-", sequence=reverse_complement(segment.sequence), KC=segment.KC)

    # create edges
    edges = [line for line in gfa_graph.lines if isinstance(
        line, gfapy.Line) and line.record_type == "L"]
    for edge in edges:
        src = edge.from_segment.name
        dist = edge.to_segment.name

        if edge.from_orient == "-":
            src = src + "-"
        else:
            src = src + "+"
        if edge.to_orient == "-":
            dist = dist + "-"
        else:
            dist = dist + "+"

        overlap = int(str(edge.overlap).replace("M", ""))
        weight = len(directed_graph.nodes[dist]['sequence']) - overlap

        directed_graph.add_edge(src, dist, overlap=overlap, weight=weight)

        if "-" in src:
            src = src[:-1]+"+"
        else:
            src = src[:-1] + "-"

        if "-" in dist:
            dist = dist[:-1] + "+"
        else:
            dist = dist[:-1] + "-"

        weight = len(
            directed_graph.nodes[src]['sequence']) - overlap
        directed_graph.add_edge(dist, src, overlap=overlap, weight=weight)

    return directed_graph

##### compare downstream and upstreams seperately for similarity (new job)
#### final method for get path
def get_paths_from_big_nx_graph_4(directed_graph, amr_gene_node, len_amr, up_down, params, stop_flag, file_name):

    stop_flag.clear()
    print("flag: ", stop_flag.is_set())
    threshold = params.neighbourhood_length
    cluster_threshold = params.max_down_up_paths

    source = f"{params.output_dir}/clustered_{amr_gene_node}_{len_amr}_{up_down}.fasta"
    destination = f"{params.output_dir}/final_down_up/clustered_{amr_gene_node}_{len_amr}_{up_down}.fasta"

    if os.path.exists(destination):
        print(f"yes {amr_gene_node} in {up_down} found ")
        file_name["file_name"] = destination
        return

    max_number_seq_for_cdhit = params.max_number_seq_for_cdhit
    radius = threshold - len_amr
    selected_paths = {}
    ego_nx_graph = nx.ego_graph(
        directed_graph, amr_gene_node, radius=radius, distance="weight")


    nodes_one_step_more = set(ego_nx_graph.nodes)


    for node in ego_nx_graph.nodes:
        adjacent_nodes = set(directed_graph.successors(node))
        nodes_one_step_more.update(adjacent_nodes)
    ego_nx_graph = directed_graph.subgraph(nodes_one_step_more)
    if len(nodes_one_step_more)==1 :
        file_name["file_name"] = ""
        return

    print("yes len of ego_nx_graph : ", len(ego_nx_graph.nodes))
    paths = {tuple([amr_gene_node]):len_amr}

    counter_for_paths = 0
    print("flag: ", stop_flag.is_set())
    while not stop_flag.is_set() and len(paths)>0:
        #while(len(paths)>0):
            p, w = paths.popitem()
            #print(p, w)
            print("paths : ", len(paths))

            #print(ego_nx_graph.out_edges(p[-1], data = True))
            neighbor = set(ego_nx_graph.successors(p[-1]))
            #print("neighbor: ", neighbor)
            if(neighbor):
                check_loop = all(item in set(p) for item in neighbor)
                #print("check loop: ", check_loop)
                if(not check_loop):
                    for edge in ego_nx_graph.out_edges(p[-1], data = True):
                        #print("n : ", edge[1])
                        #print("yes or not in a loop: ", edge[1] in p)
                        if(edge[1] not in p):
                            new_w = w + edge[2]['weight']
                            new_p = p + tuple([edge[1]])
                            if( new_w < threshold):
                                paths[new_p] = new_w
                            else:
                                selected_paths[new_p] = get_sequnce_path(ego_nx_graph, new_p, threshold, len_amr, up_down)
                else:
                    selected_paths[p] = get_sequnce_path(ego_nx_graph, p, threshold, len_amr, up_down)
            else:
                selected_paths[p] = get_sequnce_path(ego_nx_graph, p, threshold, len_amr, up_down)
            #print("len selected_paths: ", len(selected_paths))
            if(len(selected_paths) > max_number_seq_for_cdhit):
                counter_for_paths = counter_for_paths + 1
                #print("counter_for_paths : " , counter_for_paths)
                number_cluster = check_similarity_down_up_streams(selected_paths, f"{params.output_dir}/{amr_gene_node}.fasta", source, params.similarity)
                selected_paths = {}

                if(cluster_threshold != -1):
                    if(number_cluster > cluster_threshold):
                        paths = {}


    #print("len selected_paths: end ", len(selected_paths))
    if(len(selected_paths)):
        check_similarity_down_up_streams(selected_paths, f"{params.output_dir}/{amr_gene_node}.fasta", source, params.similarity)

    ##### move to final paths

    shutil.move(source, destination)

    file_name["file_name"] =  destination
    #print(file_name)
    #return destination
def merge_upstream_amr_downstream_metacherchant(upstream_paths_file, downstream_paths_file, amr, amr_name, amr_seq, len_before_amr, len_after_amr, params, seq_file):
    threshold = params.neighbourhood_length

    mergepaths = {}
    if (upstream_paths_file != "" and downstream_paths_file == ""):
        # Initialize an empty dictionary to store the data
        upstream_paths = {}

        # Open the file and read it line by line
        with open(upstream_paths_file, 'r') as file:
            while True:
                # Read the first line (key)
                key_line = file.readline().strip()

                # If the line is empty, break the loop (EOF)
                if not key_line:
                    break

                # Remove the "> " at the start and convert the string to a tuple
                key_tuple = eval(key_line.replace("> ", ""))

                # Read the second line (value)
                value_line = file.readline().strip()

                # Add the key-value pair to the dictionary
                upstream_paths[key_tuple] = value_line


        for up_path, seq in upstream_paths.items():
            #print("len temp:" , len(temp_mergepaths))
            new_path = up_path[::-1][:-1]
            new_path = new_path + (tuple(amr),)
            new_seq = seq + amr_seq.lower()
            print("new path :", type(new_path), new_path)
            mergepaths[new_path] = new_seq

        check_for_similarity(mergepaths, f"{params.output_dir}/final_result/{amr_name}.fasta", seq_file, params.similarity)

    elif (upstream_paths_file == "" and downstream_paths_file != ""):

        # Initialize an empty dictionary to store the data
        downstream_paths = {}

        # Open the file and read it line by line
        with open(downstream_paths_file, 'r') as file:
            while True:
                # Read the first line (key)
                key_line = file.readline().strip()

                # If the line is empty, break the loop (EOF)
                if not key_line:
                    break

                # Remove the "> " at the start and convert the string to a tuple
                key_tuple = eval(key_line.replace("> ", ""))

                # Read the second line (value)
                value_line = file.readline().strip()

                # Add the key-value pair to the dictionary
                downstream_paths[key_tuple] = value_line
        for down_path, seq in downstream_paths.items():
            new_path = (tuple(amr),)
            new_path = new_path + down_path[1:]
            new_seq = amr_seq.lower() + seq
            mergepaths[new_path] = new_seq

        check_for_similarity(mergepaths, f"{params.output_dir}/final_result/{amr_name}.fasta", seq_file, params.similarity)

    elif (upstream_paths_file != "" and downstream_paths_file != ""):
        # Initialize an empty dictionary to store the data
        upstream_paths = {}

        # Open the file and read it line by line
        with open(upstream_paths_file, 'r') as file:
            while True:
                # Read the first line (key)
                key_line = file.readline().strip()

                # If the line is empty, break the loop (EOF)
                if not key_line:
                    break

                # Remove the "> " at the start and convert the string to a tuple
                key_tuple = eval(key_line.replace("> ", ""))

                # Read the second line (value)
                value_line = file.readline().strip()

                # Add the key-value pair to the dictionary
                upstream_paths[key_tuple] = value_line

        # Initialize an empty dictionary to store the data
        downstream_paths = {}

        # Open the file and read it line by line
        with open(downstream_paths_file, 'r') as file:
            while True:
                # Read the first line (key)
                key_line = file.readline().strip()

                # If the line is empty, break the loop (EOF)
                if not key_line:
                    break

                # Remove the "> " at the start and convert the string to a tuple
                key_tuple = eval(key_line.replace("> ", ""))

                # Read the second line (value)
                value_line = file.readline().strip()

                # Add the key-value pair to the dictionary
                downstream_paths[key_tuple] = value_line



            for up_path, up_seq in upstream_paths.items():
                for down_path, down_seq in downstream_paths.items():
                    new_path = up_path[::-1][:-1]
                    new_path = new_path + (tuple(amr),)
                    new_path = new_path + down_path[1:]
                    new_seq = up_seq + amr_seq.lower() + down_seq
                    mergepaths[new_path] = new_seq

            check_for_similarity(mergepaths, f"{params.output_dir}/final_result/{amr_name}.fasta", seq_file, params.similarity)

    else:
        #import pdb; pdb.set_trace()
        new_path = (tuple(amr),)
        new_seq = amr_seq.lower()
        mergepaths[new_path] = new_seq
        check_for_similarity(mergepaths, f"{params.output_dir}/final_result/{amr_name}.fasta", seq_file, params.similarity)

    #return seq_file

def merge_upstream_amr_downstream_5(upstream_paths_file, downstream_paths_file, amr, amr_name, amr_seq, len_before_amr, len_after_amr, params, seq_file):
    threshold = params.neighbourhood_length

    mergepaths = {}
    if (upstream_paths_file != "" and downstream_paths_file == ""):
        # Initialize an empty dictionary to store the data
        upstream_paths = {}

        # Open the file and read it line by line
        with open(upstream_paths_file, 'r') as file:
            while True:
                # Read the first line (key)
                key_line = file.readline().strip()

                # If the line is empty, break the loop (EOF)
                if not key_line:
                    break

                # Remove the "> " at the start and convert the string to a tuple
                key_tuple = eval(key_line.replace("> ", ""))

                # Read the second line (value)
                value_line = file.readline().strip()

                # Add the key-value pair to the dictionary
                upstream_paths[key_tuple] = value_line


        for up_path, seq in upstream_paths.items():
            #print("len temp:" , len(temp_mergepaths))
            new_path = up_path[::-1][:-1]
            new_path = new_path + (tuple(amr),)
            if(len_after_amr == 0):
                new_seq = seq + amr_seq[len_before_amr:].lower()
            else:
                new_seq = seq + amr_seq[len_before_amr:-len_after_amr].lower() + amr_seq[-len_after_amr:]

            if(len_after_amr > threshold):
                new_seq = new_seq[0:-(len_after_amr-threshold)]
            print("new path :", type(new_path), new_path)
            mergepaths[new_path] = new_seq




        check_for_similarity(mergepaths, f"{params.output_dir}/final_result/{amr_name}.fasta", seq_file, params.similarity)

    elif (upstream_paths_file == "" and downstream_paths_file != ""):

        # Initialize an empty dictionary to store the data
        downstream_paths = {}

        # Open the file and read it line by line
        with open(downstream_paths_file, 'r') as file:
            while True:
                # Read the first line (key)
                key_line = file.readline().strip()

                # If the line is empty, break the loop (EOF)
                if not key_line:
                    break

                # Remove the "> " at the start and convert the string to a tuple
                key_tuple = eval(key_line.replace("> ", ""))

                # Read the second line (value)
                value_line = file.readline().strip()

                # Add the key-value pair to the dictionary
                downstream_paths[key_tuple] = value_line
        for down_path, seq in downstream_paths.items():
            new_path = (tuple(amr),)
            new_path = new_path + down_path[1:]

            if(len_before_amr > threshold):
                new_seq = amr_seq[:len_before_amr][-threshold:]
            else:
                new_seq = amr_seq[:len_before_amr]

            if(len_after_amr == 0):
                new_seq = new_seq + amr_seq[len_before_amr:].lower() + seq
            else:
                new_seq = new_seq + amr_seq[len_before_amr:-len_after_amr].lower() + seq

            mergepaths[new_path] = new_seq


        check_for_similarity(mergepaths, f"{params.output_dir}/final_result/{amr_name}.fasta", seq_file, params.similarity)

    elif (upstream_paths_file != "" and downstream_paths_file != ""):
        # Initialize an empty dictionary to store the data
        upstream_paths = {}

        # Open the file and read it line by line
        with open(upstream_paths_file, 'r') as file:
            while True:
                # Read the first line (key)
                key_line = file.readline().strip()

                # If the line is empty, break the loop (EOF)
                if not key_line:
                    break

                # Remove the "> " at the start and convert the string to a tuple
                key_tuple = eval(key_line.replace("> ", ""))

                # Read the second line (value)
                value_line = file.readline().strip()

                # Add the key-value pair to the dictionary
                upstream_paths[key_tuple] = value_line

        # Initialize an empty dictionary to store the data
        downstream_paths = {}

        # Open the file and read it line by line
        with open(downstream_paths_file, 'r') as file:
            while True:
                # Read the first line (key)
                key_line = file.readline().strip()

                # If the line is empty, break the loop (EOF)
                if not key_line:
                    break

                # Remove the "> " at the start and convert the string to a tuple
                key_tuple = eval(key_line.replace("> ", ""))

                # Read the second line (value)
                value_line = file.readline().strip()

                # Add the key-value pair to the dictionary
                downstream_paths[key_tuple] = value_line



            for up_path, up_seq in upstream_paths.items():
                for down_path, down_seq in downstream_paths.items():
                    new_path = up_path[::-1][:-1]
                    new_path = new_path + (tuple(amr),)
                    new_path = new_path + down_path[1:]
                    if(len_after_amr==0):
                        new_seq = up_seq + \
                            amr_seq[len_before_amr:].lower() + \
                        	down_seq
                    else:
                        new_seq = up_seq + \
                            amr_seq[len_before_amr:-len_after_amr].lower() + \
                        	down_seq
                    mergepaths[new_path] = new_seq

            check_for_similarity(mergepaths, f"{params.output_dir}/final_result/{amr_name}.fasta", seq_file, params.similarity)

    else:
        #import pdb; pdb.set_trace()
        new_path = (tuple(amr),)

        if(len_before_amr > threshold):
            new_seq = amr_seq[:len_before_amr][-threshold:]
        else:
            new_seq = amr_seq[:len_before_amr]

        if(len_after_amr == 0):
            new_seq = new_seq + amr_seq[len_before_amr:].lower()
        else:
            new_seq = new_seq + amr_seq[len_before_amr:-len_after_amr].lower() + amr_seq[-len_after_amr:]

        if(len_after_amr > threshold):
            new_seq = new_seq[0:-(len_after_amr-threshold)]

        mergepaths[new_path] = new_seq
        check_for_similarity(mergepaths, f"{params.output_dir}/final_result/{amr_name}.fasta", seq_file, params.similarity)

    #return seq_file



def check_add_temp_for_similarity(temp_mergepaths, input_file, output_file, similarity):

    with open(input_file, 'w') as file:
        for path, seq in temp_mergepaths.items():
            file.write("> " + str(path) + "\n")
            file.write(str(seq)+"\n")

    #os.remove(output_file)

    cdhit_command = [
        'cd-hit',                   # Command to execute
        '-i', input_file,           # Input file
        '-o', output_file,          # Output file
        '-c', str(similarity),                # Sequence identity threshold (e.g., 90%)
        '-n', '5',                   # Word length (default is 5)
        '-T', '14'
    ]
    # Execute the command
    try:
        subprocess.run(cdhit_command, check=True)
        print(f"CD-HIT completed for {input_file}.")
    except subprocess.CalledProcessError as e:
        print(f"Error running CD-HIT on {input_file}: {e}")




    #os.remove(input_file)
    shutil.copyfile(output_file, input_file)

def check_similarity_down_up_streams(paths, input_file, output_file, similarity):

    if(os.path.exists(input_file)):
        with open(input_file, 'w') as file:
            for path, seq in paths.items():
                file.write("> " + str(path) + "\n")
                file.write(str(seq)+"\n")

        cdhit_command = [
            'cd-hit',                   # Command to execute
            '-i', input_file,           # Input file
            '-o', "temp_file.fasta",          # Output file
            '-c', str(similarity),                # Sequence identity threshold (e.g., 90%)
            '-n', '5',                   # Word length (default is 5)
            '-T', '14'
        ]
        # Execute the command
        try:
            subprocess.run(cdhit_command, check=True)
            print(f"CD-HIT completed for {input_file}.")
        except subprocess.CalledProcessError as e:
            print(f"Error running CD-HIT on {input_file}: {e}")

        with open('temp_file.fasta', 'r') as file:
           content_to_append = file.read()


        with open(output_file, 'a') as file:
           file.write(content_to_append)

        # Copy the file
        shutil.copy(output_file, input_file)

        cdhit_command = [
            'cd-hit',                   # Command to execute
            '-i', input_file,           # Input file
            '-o', output_file,          # Output file
            '-c', str(similarity),                # Sequence identity threshold (e.g., 90%)
            '-n', '5',                   # Word length (default is 5)
            '-T', '14'
        ]
        # Execute the command
        try:
            subprocess.run(cdhit_command, check=True)
            print(f"CD-HIT completed for {input_file}.")
        except subprocess.CalledProcessError as e:
            print(f"Error running CD-HIT on {input_file}: {e}")

    else:
        with open(input_file, 'w') as file:
            for path, seq in paths.items():
                file.write("> " + str(path) + "\n")
                file.write(str(seq)+"\n")

        cdhit_command = [
            'cd-hit',                   # Command to execute
            '-i', input_file,           # Input file
            '-o', output_file,          # Output file
            '-c', str(similarity),                # Sequence identity threshold (e.g., 90%)
            '-n', '5',                   # Word length (default is 5)
            '-T', '14'
        ]
        # Execute the command
        try:
            subprocess.run(cdhit_command, check=True)
            print(f"CD-HIT completed for {input_file}.")
        except subprocess.CalledProcessError as e:
            print(f"Error running CD-HIT on {input_file}: {e}")


    number_cluster = 0
    with open(output_file, 'r') as file:
        number_cluster = sum(1 for line in file) / 2

    #print(f"number of cluster : {number_cluster}")
    return number_cluster


def check_for_similarity(mergepaths, input_file, output_file, similarity):
    with open(input_file, 'a') as file:
        for path, seq in mergepaths.items():
            file.write("> " + str(path) + "\n")
            file.write(str(seq)+"\n")
    # Define the CD-HIT command
    cdhit_command = [
        'cd-hit',                   # Command to execute
        '-i', input_file,           # Input file
        '-o', output_file,          # Output file
        '-c', str(similarity),                # Sequence identity threshold (e.g., 90%)
        '-n', '5'                   # Word length (default is 5)
    ]

    # Execute the command
    try:
        subprocess.run(cdhit_command, check=True)
        print(f"CD-HIT completed for {input_file}.")
    except subprocess.CalledProcessError as e:
        print(f"Error running CD-HIT on {input_file}: {e}")

def write_paths_info_to_file(paths_info_list, paths_info_file):
    """ """

    with open(paths_info_file, "a") as fd:
        writer = csv.writer(fd)
        for path_info in paths_info_list:
            #print(path_info)
            writer.writerow(
                [
                    path_info["sequence"],
                    path_info["node"],
                    path_info["coverage"],
                    path_info["start"],
                    path_info["end"],
                ]
            )



def neighborhood_sequence_extraction(
    directed_graph, reverse_directed_graph, params, amr_info
):
    #print(params.output_dir)
    # Define a flag for stopping the loop
    stop_flag_downstream = threading.Event()
    stop_flag_upstream = threading.Event()

    duration = params.extraction_timeout * 60

    file_path, amr_paths = amr_info
    #LOG.debug("amr_file = " + file_path)
    amr_name = file_path.split("/")[-1].split(".")[0]
    LOG.info(f"exteracting downstream and upstream paths for : {amr_name}")
    #print("amr_name :", amr_name)
    #print("len amr paths :", len(amr_paths))
    # csv file for path_info
    os.makedirs(os.path.join(params.output_dir, SEQ_DIR_NAME,
        SEQ_DIR_NAME + "_" + str(params.neighbourhood_length), "paths_info"), exist_ok=True)
    paths_info_file = os.path.join(
        params.output_dir,
        SEQ_DIR_NAME,
        SEQ_DIR_NAME + "_" + str(params.neighbourhood_length),
        "paths_info",
        SEQ_NAME_PREFIX + amr_name + "_" + str(params.neighbourhood_length)
        + "_"
        + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
        + ".csv",
    )

    with open(paths_info_file, "a") as fd:
        writer = csv.writer(fd)
        writer.writerow(["sequence", "node", "coverage", "start", "end"])

    # Extract the sequenc of AMR neighborhood
    LOG.debug(f"Calling extract_neighborhood_sequence for {os.path.basename(amr_name)}...")
    seq_counter = 0
    output_name = SEQ_NAME_PREFIX + amr_name
    seq_output_dir = os.path.join(params.output_dir, SEQ_DIR_NAME,
        SEQ_DIR_NAME + "_" + str(params.neighbourhood_length),AMR_SEQ_DIR)
    os.makedirs(seq_output_dir, exist_ok=True)
    threshold = params.neighbourhood_length
    seq_file = os.path.join(
        seq_output_dir,
        output_name
        + "_"
        + str(threshold)
        + "_"
        + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
        + ".txt",
    )
    for amr_path in amr_paths:
        find_downstream = True
        find_upstream = True
        #print(amr_path)
        # amr_gene
        amr_gene = [pair[0]+pair[1]
                    for pair in zip(amr_path['nodes'], amr_path['orientations'])]
        print("amr gene : ", amr_gene)
        #len_after_amr = len(
        #    directed_graph.nodes[amr_gene[-1]]['sequence']) - amr_path['end_pos']
        #len_before_amr = amr_path['start_pos'] - 1

        if(amr_path['end_pos'] == 0):
            len_after_amr = 0
        else:
            len_after_amr = len(
                directed_graph.nodes[amr_gene[-1]]['sequence']) - amr_path['end_pos']

        if(amr_path['start_pos'] == 0):
            len_before_amr = 0
        else:
            len_before_amr = amr_path['start_pos'] - 1


        print("len after: ", len_after_amr)
        print("len before: ", len_before_amr)


        if len_after_amr >= params.neighbourhood_length:
            find_downstream = False

        if len_before_amr >= params.neighbourhood_length:
            find_upstream = False

        #print(find_downstream)
        #print(find_upstream)



        downstream_paths_file = {"file_name" : "" }
        upstream_paths_file = {"file_name" : "" }

        if params.assembler == "metacherchant":
        	#read amr_seq directly from the amr file
        	amr_seq, _ = retrieve_AMR(file_path)
        	amr_seq = amr_seq.rstrip('\n')
        else:
        	amr_seq = directed_graph.nodes[amr_gene[0]]['sequence']
        	if(len(amr_gene)>1):
        		for gene_index in range(1, len(amr_gene)):
        			amr_seq = amr_seq + directed_graph.nodes[amr_gene[gene_index]]['sequence'][directed_graph[amr_gene[gene_index-1]][amr_gene[gene_index]]['overlap']:]

        # amr_seq = amr_seq[len_before_amr:-len_after_amr]
        # downstream

        if find_downstream:
            print("into down")
            amr_gene_node = amr_gene[-1]

            loop_thread = threading.Thread(target=get_paths_from_big_nx_graph_4, args=(directed_graph, amr_gene_node, len_after_amr, "down", params, stop_flag_downstream, downstream_paths_file))
            loop_thread.start()
            # Wait for the thread to finish or until the time limit expires
            start_time = time.time()
            print(f"duration : {duration}")
            while time.time() - start_time < duration:
                if not loop_thread.is_alive():  # If the thread has finished naturally, break
                    #print("Loop finished early.")
                    break
                time.sleep(0.001)  # Check periodically

            # If time expires, set the stop flag to stop the thread
            if loop_thread.is_alive():
                #print("Time limit reached. Stopping the loop.")
                stop_flag_downstream.set()
            loop_thread.join()

            #downstream_paths_file = get_paths_from_big_nx_graph_4(
            #    directed_graph, amr_gene_node, len_after_amr, "down", params, stop_flag)
            print("done downstream")

        # upstream
        if find_upstream:
            print("into up")
            amr_gene_node = amr_gene[0]
            #stop_flag.clear()
            print("flag: ", stop_flag_upstream.is_set())
            loop_thread = threading.Thread(target=get_paths_from_big_nx_graph_4, args=(reverse_directed_graph, amr_gene_node, len_before_amr, "up", params, stop_flag_upstream, upstream_paths_file))
            loop_thread.start()
            # Wait for the thread to finish or until the time limit expires
            start_time = time.time()
            while time.time() - start_time < duration:
                if not loop_thread.is_alive():  # If the thread has finished naturally, break
                    #print("Loop finished early.")
                    break
                time.sleep(0.001)  # Check periodically

            # If time expires, set the stop flag to stop the thread
            if loop_thread.is_alive():
                #print("Time limit reached. Stopping the loop.")
                stop_flag_upstream.set()
            loop_thread.join()

            #upstream_paths_file = get_paths_from_big_nx_graph_4(
            #    reverse_directed_graph, amr_gene_node, len_before_amr, "up", params, stop_flag)
            print("done upstream")

        print(f"downstream file :{downstream_paths_file}")
        print(f"upstream file : {upstream_paths_file}")
        if params.assembler == "metacherchant":
        	merge_upstream_amr_downstream_metacherchant(
        		upstream_paths_file=upstream_paths_file["file_name"], downstream_paths_file=downstream_paths_file["file_name"], amr=amr_gene, amr_name=amr_name,
        		amr_seq=amr_seq, len_before_amr=len_before_amr, len_after_amr=len_after_amr, params=params, seq_file=seq_file)
        else:
        	merge_upstream_amr_downstream_5(
            		upstream_paths_file=upstream_paths_file["file_name"], downstream_paths_file=downstream_paths_file["file_name"], amr=amr_gene, amr_name=amr_name,
            		amr_seq=amr_seq, len_before_amr=len_before_amr, len_after_amr=len_after_amr, params=params, seq_file=seq_file)


    if params.assembler != "metacherchant":
        path_info_list = create_paths_info_list(seq_file, directed_graph, amr_paths, params.neighbourhood_length, params.max_kmer_size, params.assembler)

        write_paths_info_to_file(path_info_list, paths_info_file)

    return seq_file, paths_info_file



def create_paths_info_list(seq_file, ego_graph, amr_paths, threshold, max_kmer_size, assembler):
    paths_info_list = []
    paths = {}

    # Open the file and read it line by line
    with open(seq_file, 'r') as file:
        while True:
            key_value = file.readline().strip()
            if not key_value:
                break
            path = eval(key_value.replace("> ", ""))
            seq = file.readline().strip()
            paths[path] = seq
            #paths.append(path)
    counter = 0
    for path, seq in paths.items():
        seq = seq.upper()
        end = -1
        #print(f"path : {path} , {type(path)}, amr_gene, {amr_gene}, {type(amr_gene)}")

        for index, element in enumerate(path):
            if isinstance(element, tuple):
                #print(f"The first tuple is at index: {index}")
                index_amr_in_path = index
                for amr_path in amr_paths:
                    #print(amr_path)
                    # amr_gene
                    amr_gene = [pair[0]+pair[1]
                        for pair in zip(amr_path['nodes'], amr_path['orientations'])]
                    if(tuple(amr_gene)==element):
                        if(amr_path['end_pos'] == 0):
                            len_after_amr = 0
                        else:
                            len_after_amr = len(
                                ego_graph.nodes[amr_gene[-1]]['sequence']) - amr_path['end_pos']

                        if(amr_path['start_pos'] == 0):
                            len_before_amr = 0
                        else:
                            len_before_amr = amr_path['start_pos'] - 1
                break

        print(f"len_befor : {len_before_amr}")
        print(f"len after : {len_after_amr}")
        print(f"amr : {amr_gene}")
        #index_amr_in_path = path.index(amr_gene)
        #print("index amr : ", index_amr_in_path)

        is_upstreams = (index_amr_in_path > 0)
        is_downstreams = (len(path)-1) > index_amr_in_path

        #print(is_upstreams)
        #print(is_downstreams)

        #### upstreams
        if(is_upstreams):
            for up_path_index in range(index_amr_in_path):
                node = path[up_path_index]
                coverage = calculate_coverage(ego_graph.nodes[path[up_path_index]], max_kmer_size, path[up_path_index], assembler)
                start = end + 1
                if(type(path[up_path_index+1]) == tuple):
                    #print(path[up_path_index + 1])
                    #print("path[up_path_index + 1][0] : ", path[up_path_index + 1][0])
                    end = seq.find(ego_graph.nodes[path[up_path_index + 1][0]]['sequence'][0:len_before_amr]) - 1
                    #print("end : ", end)
                else:
                    #print(path[up_path_index + 1])
                    end = seq.find(ego_graph.nodes[path[up_path_index + 1]]['sequence']) - 1
                    #print(end)
                info = {"sequence" : counter,
                    "node" : node,
                    "coverage" : coverage,
                    "start" : start,
                    "end" : end}
                paths_info_list.append(info)

        #print(len(path[index_amr_in_path]))
        ##### amr
        if(len(path[index_amr_in_path]) == 1):
            print(f"len_befor : {len_before_amr}")
            print(f"len after : {len_after_amr}")
            node = path[index_amr_in_path]
            coverage = calculate_coverage(ego_graph.nodes[path[index_amr_in_path][0]], max_kmer_size, path[index_amr_in_path][0], assembler)
            ####upstream
            if(len_before_amr > 0):
                start = end + 1
                end = (end + len_before_amr) if(len_before_amr <= threshold) else (end + threshold)
                #print("info : ", info)
                info = {"sequence" : counter,
                        "node" : node,
                        "coverage" : coverage,
                        "start" : start,
                        "end" : end}
                paths_info_list.append(info)

            ####amr
            start = end + 1
            end = end + len(ego_graph.nodes[path[index_amr_in_path][0]]['sequence'][len_before_amr:]) if (len_after_amr==0) else end + len(ego_graph.nodes[path[index_amr_in_path][0]]['sequence'][len_before_amr:- len_after_amr])
            info = {"sequence" : counter,
                    "node" : node,
                    "coverage" : coverage,
                    "start" : start,
                    "end" : end}
            #print("info : ", info)
            paths_info_list.append(info)

            ###downstream
            if(len_after_amr>0):
                start = end + 1
                end = end + threshold  if (len_after_amr >= threshold) else end + len_after_amr
                info = {"sequence" : counter,
                        "node" : node,
                        "coverage" : coverage,
                        "start" : start,
                        "end" : end}
                #print("info : ", info)
                paths_info_list.append(info)

        elif len(path[index_amr_in_path]) == 2 :
            ##### upstream
            node = path[index_amr_in_path][0]
            coverage = calculate_coverage(ego_graph.nodes[path[index_amr_in_path][0]], max_kmer_size, path[index_amr_in_path][0], assembler)
            if(len_before_amr > 0):
                start = end + 1
                end = (end + len_before_amr) if(len_before_amr <= threshold) else (end + threshold)
                info = {"sequence" : counter,
                        "node" : node,
                        "coverage" : coverage,
                        "start" : start,
                        "end" : end}
                paths_info_list.append(info)

            #####amr
            start = end + 1
            end = end + len(ego_graph.nodes[path[index_amr_in_path][0]]['sequence'][len_before_amr:])
            info = {"sequence" : counter,
                    "node" : node,
                    "coverage" : coverage,
                    "start" : start,
                    "end" : end}
            paths_info_list.append(info)

            node = path[index_amr_in_path][-1]
            coverage = calculate_coverage(ego_graph.nodes[path[index_amr_in_path][-1]], max_kmer_size, path[index_amr_in_path][-1], assembler)


            start = end + 1
            #import pdb; pdb.set_trace()
            end = end + len(ego_graph.nodes[path[index_amr_in_path][-1]]['sequence']) - len_after_amr - int(ego_graph[path[index_amr_in_path][-2]][path[index_amr_in_path][-1]]['overlap'])
            info = {"sequence" : counter,
                    "node" : node,
                    "coverage" : coverage,
                    "start" : start,
                    "end" : end}
            paths_info_list.append(info)

            #### downstream
            if(len_after_amr>0):
                start = end + 1
                if(index_amr_in_path < len(path) - 1):
                    end = (end + len_after_amr) if(len_after_amr <= threshold) else (end + threshold)
                else:
                    end = len(seq) - 1
                info = {"sequence" : counter,
                        "node" : node,
                        "coverage" : coverage,
                        "start" : start,
                        "end" : end}
                paths_info_list.append(info)

        else:
            ##### upstream
            node = path[index_amr_in_path][0]
            coverage = calculate_coverage(ego_graph.nodes[path[index_amr_in_path][0]], max_kmer_size, path[index_amr_in_path][0], assembler)
            if(len_before_amr > 0):
                start = end + 1
                end = (end + len_before_amr) if(len_before_amr <= threshold) else (end + threshold)
                info = {"sequence" : counter,
                        "node" : node,
                        "coverage" : coverage,
                        "start" : start,
                        "end" : end}
                paths_info_list.append(info)

            #####amr
            start = end + 1
            end = end + len(ego_graph.nodes[node]['sequence'][len_before_amr:])
            info = {"sequence" : counter,
                    "node" : node,
                    "coverage" : coverage,
                    "start" : start,
                    "end" : end}
            paths_info_list.append(info)


            for amr_index in range(1, len(path[index_amr_in_path])-1) :
                node = path[index_amr_in_path][amr_index]
                coverage = calculate_coverage(ego_graph.nodes[path[index_amr_in_path][amr_index]], max_kmer_size, path[index_amr_in_path][amr_index], assembler)
                start = end + 1
                end = end + ego_graph[path[index_amr_in_path][amr_index-1]][path[index_amr_in_path][amr_index]]['weight']
                info = {"sequence" : counter,
                    "node" : node,
                    "coverage" : coverage,
                    "start" : start,
                    "end" : end}
                paths_info_list.append(info)

            node = path[index_amr_in_path][-1]
            coverage = calculate_coverage(ego_graph.nodes[path[index_amr_in_path][-1]], max_kmer_size, path[index_amr_in_path][-1], assembler)

            start = end + 1
            end = end + len(ego_graph.nodes[path[index_amr_in_path][-1]]['sequence']) - len_after_amr - int(ego_graph[path[index_amr_in_path][-2]][path[index_amr_in_path][-1]]['overlap'])
            info = {"sequence" : counter,
                    "node" : node,
                    "coverage" : coverage,
                    "start" : start,
                    "end" : end}
            paths_info_list.append(info)

            #### downstream
            if(len_after_amr>0):
                start = end + 1
                if(index_amr_in_path < len(path) - 1):
                    end = (end + len_after_amr) if(len_after_amr <= threshold) else (end + threshold)
                else:
                    end = len(seq) - 1
                info = {"sequence" : counter,
                        "node" : node,
                        "coverage" : coverage,
                        "start" : start,
                        "end" : end}
                paths_info_list.append(info)

        #import pdb

        #pdb.set_trace()
        #### downstreams
        #print("is_downstrams : ", is_downstreams)
        if(is_downstreams):
            #print("index_amr_in_path : ", index_amr_in_path)
            #print("len path :", len(path))

            for down_path_index in range(index_amr_in_path+1, len(path)-1):
                #print(down_path_index)
                node = path[down_path_index]
                coverage = calculate_coverage(ego_graph.nodes[path[down_path_index]], max_kmer_size, path[down_path_index], assembler)
                start = end + 1
                if(down_path_index==index_amr_in_path+1):
                    end = end + ego_graph[path[down_path_index-1][-1]][path[down_path_index]]['weight']
                    #print("end in the if : ", end)
                else:
                    end = end + ego_graph[path[down_path_index-1]][path[down_path_index]]['weight']
                    #print("end out of the if : ", end)
                info = {"sequence" : counter,
                    "node" : node,
                    "coverage" : coverage,
                    "start" : start,
                    "end" : end}
                paths_info_list.append(info)

            node = path[-1]
            coverage = calculate_coverage(ego_graph.nodes[path[-1]], max_kmer_size, path[-1], assembler)
            start = end + 1
            end = len(seq) - 1
            info = {"sequence" : counter,
                    "node" : node,
                    "coverage" : coverage,
                    "start" : start,
                    "end" : end}
            paths_info_list.append(info)
        counter = counter + 1

    #print("paths info list : ", paths_info_list)
    return paths_info_list





if __name__ == "__main__":
    pass


def sequence_neighborhood_main(
        params,
        gfa_file,
        amr_seq_align_info,
        debug: bool
):
    sequence_dir = os.path.join(
        params.output_dir,
        SEQ_DIR_NAME,
        f'{SEQ_DIR_NAME}_{params.neighbourhood_length}',
    )
    os.makedirs(sequence_dir, exist_ok=True)

    os.makedirs(f"{params.output_dir}/final_result", exist_ok=True)
    os.makedirs(f"{params.output_dir}/final_down_up", exist_ok=True)
    # Read GFA file
    gfa_graph = gfapy.Gfa.from_file(gfa_file)
    directed_graph = create_directed_graph_nx(gfa_graph)

    reverse_directed_graph = directed_graph.reverse()
    for edge in reverse_directed_graph.edges:
        reverse_directed_graph[edge[0]][edge[1]]['weight'] = len(
            reverse_directed_graph.nodes[edge[1]]['sequence']) - reverse_directed_graph[edge[0]][edge[1]]['overlap']

    #import pdb; pdb.set_trace()
    # If running single threaded do not add any overhead using multiprocessing pool
    if params.num_cores == 1:
        lists = list()
        for x in amr_seq_align_info:
            file_path, amr_paths = x
            #LOG.debug("amr_file = " + file_path)
            amr_name = file_path.split("/")[-1].split(".")[0]
            #if(amr_name == "TEM-181") :
            lists.append(neighborhood_sequence_extraction(directed_graph, reverse_directed_graph, params, x))
    else:
        p_extraction = partial(
            neighborhood_sequence_extraction,
            directed_graph, reverse_directed_graph, params
        )
        with Pool(params.num_cores) as p:
            lists = list(p.imap(p_extraction, amr_seq_align_info))
    seq_files, path_info_files = zip(*lists)
    LOG.info("delete files")
    for item in os.listdir(params.output_dir):
        item_path = os.path.join(params.output_dir, item)

        # Check if it is a file, and if so, delete it
        if os.path.isfile(item_path):
            if not (item_path.endswith('.log') or item_path.endswith('.txt')):
                os.remove(item_path)



    shutil.rmtree(f"{params.output_dir}/final_result")
    shutil.rmtree(f"{params.output_dir}/final_down_up")

    return seq_files, path_info_files
