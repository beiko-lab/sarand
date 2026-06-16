import sys
import datetime
import gc
import threading
import time
import csv
import shutil
import multiprocessing
from pathlib import Path

import gfapy
import networkx as nx

from sarand.util.logger import LOG
from sarand.config import SEQ_DIR_NAME, SEQ_NAME_PREFIX, NEIGHBORHOOD_SEQ_DIR
from sarand.external.cdhit import Cdhit


def read_paths_file(file_path):
    """Read a paths file of alternating "> <path-tuple>" / sequence lines into a dict."""
    paths = {}
    with open(file_path, 'r') as fd:
        while True:
            key_line = fd.readline().strip()
            # An empty line marks EOF
            if not key_line:
                break
            path = eval(key_line.replace("> ", ""))
            paths[path] = fd.readline().strip()
    return paths


def write_paths_file(file_path, paths, mode='w'):
    """Write a {path: sequence} dict as alternating "> <path>" / sequence lines."""
    with open(file_path, mode) as fd:
        for path, seq in paths.items():
            fd.write("> " + str(path) + "\n")
            fd.write(str(seq) + "\n")

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
    else:
        LOG.error(
            "no way of calculating node coverage has been defined for this assembler!"
        )
        sys.exit()
    return coverage


def get_sequence_path(directed_graph, path, threshold, len_target, up_down):
    """
    Reconstruct the nucleotide sequence spelled out by a node path, trimmed to
    ``threshold`` bases on the requested side of the target gene.
    Parameters:
        directed_graph: the assembly graph
        path: the ordered list of node names making up the path
        threshold: maximum number of bases to keep
        len_target: length of the target-gene portion contributed by the first node
        up_down: "down" for a downstream path, "up" for an upstream path
    Return:
        the (trimmed) path sequence as a string
    """
    sequence = ""
    if (up_down == "down"):
        sequence = directed_graph.nodes[path[0]]['sequence'][-len_target:]
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
            directed_graph.nodes[path[0]]['sequence'][0:len_target]
        sequence = sequence[-threshold:]
    return sequence


def reverse_complement(sequence):
    """Return the reverse complement of a DNA ``sequence``."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_sequence = sequence[::-1]
    reverse_complement_sequence = ''.join(
        [complement[base] for base in reverse_sequence])
    return reverse_complement_sequence


def create_directed_graph_nx(gfa_graph):
    """
    Build a networkx directed graph from a gfapy GFA graph.

    Each segment becomes two nodes ("+" for the segment sequence and "-" for its
    reverse complement) and every GFA link becomes a pair of directed edges
    (forward and the reverse-complement direction) weighted by the non-overlapping
    sequence length.
    """
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
def get_paths_from_big_nx_graph_4(directed_graph, target_gene_node, len_target, up_down, params, stop_flag, file_name):
    """
    Enumerate the paths radiating from a target-gene node out to the neighborhood
    length, cluster their sequences with cd-hit, and record the resulting file.

    The search is bounded by ``stop_flag`` (set by the caller on timeout). The
    clustered paths are written to a fasta under ``final_down_up`` and the path of
    that file is stored in ``file_name["file_name"]`` (empty string if no paths).
    Parameters:
        directed_graph: the (forward or reversed) assembly graph to traverse
        target_gene_node: the node to start the traversal from
        len_target: bases of the target gene contributed by the start node
        up_down: "down" for downstream traversal, "up" for upstream
        params: the parsed CLI parameters
        stop_flag: threading.Event used to abort the traversal on timeout
        file_name: dict whose "file_name" key receives the output path
    """
    stop_flag.clear()
    print("flag: ", stop_flag.is_set())
    threshold = params.neighborhood_length

    source = f"{params.output_dir}/clustered_{target_gene_node}_{len_target}_{up_down}.fasta"
    destination = f"{params.output_dir}/final_down_up/clustered_{target_gene_node}_{len_target}_{up_down}.fasta"

    if Path(destination).exists():
        print(f"yes {target_gene_node} in {up_down} found ")
        file_name["file_name"] = destination
        return

    radius = threshold - len_target
    selected_paths = {}
    ego_nx_graph = nx.ego_graph(
        directed_graph, target_gene_node, radius=radius, distance="weight")


    nodes_one_step_more = set(ego_nx_graph.nodes)


    for node in ego_nx_graph.nodes:
        adjacent_nodes = set(directed_graph.successors(node))
        nodes_one_step_more.update(adjacent_nodes)
    ego_nx_graph = directed_graph.subgraph(nodes_one_step_more)
    if len(nodes_one_step_more)==1 :
        file_name["file_name"] = ""
        return

    print("yes len of ego_nx_graph : ", len(ego_nx_graph.nodes))
    paths = {tuple([target_gene_node]):len_target}

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
                                selected_paths[new_p] = get_sequence_path(ego_nx_graph, new_p, threshold, len_target, up_down)
                else:
                    selected_paths[p] = get_sequence_path(ego_nx_graph, p, threshold, len_target, up_down)
            else:
                selected_paths[p] = get_sequence_path(ego_nx_graph, p, threshold, len_target, up_down)

    if(len(selected_paths)):
        check_similarity_down_up_streams(selected_paths, f"{params.output_dir}/{target_gene_node}.fasta", source, params.deduplication_identity)

    ##### move to final paths

    shutil.move(source, destination)

    file_name["file_name"] =  destination
    #print(file_name)
    #return destination

def merge_upstream_target_downstream_5(upstream_paths_file, downstream_paths_file, target, target_name, target_seq, len_before_target, len_after_target, params, seq_file):
    """
    Join the upstream paths, the target gene, and the downstream paths into full
    neighborhood sequences and cluster them into the per-target result fasta.

    The target gene is written in lower case (with any portion already covered by
    the gene's own node trimmed via ``len_before_target``/``len_after_target``) so
    it can be located again during annotation.
    Parameters:
        upstream_paths_file: file of clustered upstream paths ("" if none)
        downstream_paths_file: file of clustered downstream paths ("" if none)
        target: the list of node names making up the target gene
        target_name: the target gene name (used for the output file)
        target_seq: the target gene nucleotide sequence
        len_before_target: bases of the start node preceding the target gene
        len_after_target: bases of the end node following the target gene
        params: the parsed CLI parameters
        seq_file: the cumulative per-target neighborhood sequence file
    """
    threshold = params.neighborhood_length

    mergepaths = {}
    if (upstream_paths_file != "" and downstream_paths_file == ""):
        upstream_paths = read_paths_file(upstream_paths_file)
        for up_path, seq in upstream_paths.items():
            new_path = up_path[::-1][:-1]
            new_path = new_path + (tuple(target),)
            if(len_after_target == 0):
                new_seq = seq + target_seq[len_before_target:].lower()
            else:
                new_seq = seq + target_seq[len_before_target:-len_after_target].lower() + target_seq[-len_after_target:]

            if(len_after_target > threshold):
                new_seq = new_seq[0:-(len_after_target-threshold)]
            print("new path :", type(new_path), new_path)
            mergepaths[new_path] = new_seq

    elif (upstream_paths_file == "" and downstream_paths_file != ""):
        downstream_paths = read_paths_file(downstream_paths_file)
        for down_path, seq in downstream_paths.items():
            new_path = (tuple(target),)
            new_path = new_path + down_path[1:]

            if(len_before_target > threshold):
                new_seq = target_seq[:len_before_target][-threshold:]
            else:
                new_seq = target_seq[:len_before_target]

            if(len_after_target == 0):
                new_seq = new_seq + target_seq[len_before_target:].lower() + seq
            else:
                new_seq = new_seq + target_seq[len_before_target:-len_after_target].lower() + seq

            mergepaths[new_path] = new_seq

    elif (upstream_paths_file != "" and downstream_paths_file != ""):
        upstream_paths = read_paths_file(upstream_paths_file)
        downstream_paths = read_paths_file(downstream_paths_file)
        for up_path, up_seq in upstream_paths.items():
            for down_path, down_seq in downstream_paths.items():
                new_path = up_path[::-1][:-1]
                new_path = new_path + (tuple(target),)
                new_path = new_path + down_path[1:]
                if(len_after_target==0):
                    new_seq = up_seq + \
                        target_seq[len_before_target:].lower() + \
                        down_seq
                else:
                    new_seq = up_seq + \
                        target_seq[len_before_target:-len_after_target].lower() + \
                        down_seq
                mergepaths[new_path] = new_seq

    else:
        new_path = (tuple(target),)

        if(len_before_target > threshold):
            new_seq = target_seq[:len_before_target][-threshold:]
        else:
            new_seq = target_seq[:len_before_target]

        if(len_after_target == 0):
            new_seq = new_seq + target_seq[len_before_target:].lower()
        else:
            new_seq = new_seq + target_seq[len_before_target:-len_after_target].lower() + target_seq[-len_after_target:]

        if(len_after_target > threshold):
            new_seq = new_seq[0:-(len_after_target-threshold)]

        mergepaths[new_path] = new_seq

    check_for_similarity(mergepaths, f"{params.output_dir}/final_result/{target_name}.fasta", seq_file, params.deduplication_identity)

    #return seq_file



def check_similarity_down_up_streams(paths, input_file, output_file, similarity):
    """
    Cluster the given up/downstream ``paths`` with cd-hit, accumulating the unique
    sequences into ``output_file``, and return the number of clusters.
    Parameters:
        paths: a {path: sequence} dict of candidate paths
        input_file: scratch fasta written for cd-hit input
        output_file: fasta accumulating the clustered (unique) sequences
        similarity: cd-hit identity threshold
    Return:
        the number of clusters currently in ``output_file``
    """
    if Path(input_file).exists():
        write_paths_file(input_file, paths)
        Cdhit.cluster(input_file, "temp_file.fasta", similarity)

        with open('temp_file.fasta', 'r') as file:
           content_to_append = file.read()


        with open(output_file, 'a') as file:
           file.write(content_to_append)

        # Copy the file
        shutil.copy(output_file, input_file)

        Cdhit.cluster(input_file, output_file, similarity)

    else:
        write_paths_file(input_file, paths)
        Cdhit.cluster(input_file, output_file, similarity)


    number_cluster = 0
    with open(output_file, 'r') as file:
        number_cluster = sum(1 for line in file) / 2

    #print(f"number of cluster : {number_cluster}")
    return number_cluster


def check_for_similarity(mergepaths, input_file, output_file, similarity):
    """
    Append the merged neighborhood sequences to ``input_file`` and cluster them
    with cd-hit into ``output_file`` to drop near-duplicates.
    Parameters:
        mergepaths: a {path: sequence} dict of merged neighborhood sequences
        input_file: fasta the sequences are appended to (cd-hit input)
        output_file: fasta receiving the clustered result
        similarity: cd-hit identity threshold
    """
    write_paths_file(input_file, mergepaths, mode='a')
    Cdhit.cluster(input_file, output_file, similarity)

def write_paths_info_to_file(paths_info_list, paths_info_file):
    """
    Compute per-path coverage and append the per-node path info rows to the
    path-info CSV.
    Parameters:
        paths_info_list: the list of per-node info dicts produced for every path
        paths_info_file: the CSV file to append the rows to
    """
    # calculate path coverage
    coverage_list = []
    sequence_num= -1
    first_seq_num = 0
    coverage = 0
    for node_info in paths_info_list:
        if (node_info["sequence"] != sequence_num) and (coverage!=0):
            coverage_list.append(coverage/path_length)
            coverage = 0
        elif (node_info["sequence"] != sequence_num) and (coverage==0):
            first_seq_num = node_info["sequence"]
        #it is updated by each node but when it's used it has the "end" of the last node in the path with that sequence #
        path_length=node_info["end"]
        sequence_num = node_info["sequence"]
        coverage += node_info["coverage"] * (node_info["end"] - node_info["start"] + 1)
    coverage_list.append(coverage/path_length)
    # write into paths_info_file
    sequence_num = first_seq_num
    coverage_index = 0
    with open(paths_info_file, "a") as fd:
        writer = csv.writer(fd)
        for node_info in paths_info_list:
            if (node_info["sequence"] != sequence_num):
                coverage_index = coverage_index + 1
            sequence_num =  node_info["sequence"] 
            #print(path_info)
            writer.writerow(
                [
                    node_info["sequence"],
                    node_info["node"],
                    node_info["coverage"],
                    node_info["start"],
                    node_info["end"],
                    coverage_list[coverage_index]
                ]
            )



def neighborhood_sequence_extraction(
    directed_graph, reverse_directed_graph, params, target_info
):
    """
    Extract the up/downstream neighborhood sequences for a single target gene.

    For each alignment hit of the target the upstream and downstream paths are
    enumerated (each in a timed-out worker thread), merged with the target gene
    into full neighborhood sequences, and the per-node path/coverage info is
    written out.
    Parameters:
        directed_graph: the forward assembly graph
        reverse_directed_graph: the reversed assembly graph (for upstream search)
        params: the parsed CLI parameters
        target_info: a (target_file_path, target_hits) pair
    Return:
        (seq_file, paths_info_file): the neighborhood sequence and path-info files
    """
    # Define a flag for stopping the loop
    stop_flag_downstream = threading.Event()
    stop_flag_upstream = threading.Event()

    duration = params.extraction_timeout * 60

    file_path, target_hits = target_info
    #LOG.debug("target_file = " + file_path)
    target_name = file_path.split("/")[-1].split(".")[0]
    LOG.info(f"Extracting downstream and upstream paths for : {target_name}")
    #print("target_name :", target_name)
    #print("len target paths :", len(target_hits))
    # csv file for path_info
    length_dir = (
        Path(params.output_dir)
        / SEQ_DIR_NAME
    )
    paths_info_dir = length_dir / "neighborhood_paths"
    paths_info_dir.mkdir(parents=True, exist_ok=True)
    # kept as str: the returned file path is substring-matched downstream
    paths_info_file = str(
        paths_info_dir
        / (
            SEQ_NAME_PREFIX + target_name + "_" + str(params.neighborhood_length)
            + "_"
            + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
            + ".csv"
        )
    )

    with open(paths_info_file, "a") as fd:
        writer = csv.writer(fd)
        writer.writerow(["sequence", "node", "coverage", "start", "end", "path_coverage"])

    # Extract the sequence of target neighborhood
    LOG.debug(f"Calling extract_neighborhood_sequences for {Path(target_name).name}...")
    output_name = SEQ_NAME_PREFIX + target_name
    seq_output_dir = length_dir / NEIGHBORHOOD_SEQ_DIR
    seq_output_dir.mkdir(parents=True, exist_ok=True)
    threshold = params.neighborhood_length
    seq_file = str(
        seq_output_dir
        / (
            output_name
            + "_"
            + str(threshold)
            + "_"
            + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
            + ".txt"
        )
    )
    for target_hit in target_hits:
        find_downstream = True
        find_upstream = True
        #print(target_hit)
        # target_gene
        target_gene = [pair[0]+pair[1]
                    for pair in zip(target_hit['nodes'], target_hit['orientations'])]
        print("target gene : ", target_gene)
        #len_after_target = len(
        #    directed_graph.nodes[target_gene[-1]]['sequence']) - target_hit['end_pos']
        #len_before_target = target_hit['start_pos'] - 1

        if(target_hit['end_pos'] == 0):
            len_after_target = 0
        else:
            len_after_target = len(
                directed_graph.nodes[target_gene[-1]]['sequence']) - target_hit['end_pos']

        if(target_hit['start_pos'] == 0):
            len_before_target = 0
        else:
            len_before_target = target_hit['start_pos'] - 1


        print("len after: ", len_after_target)
        print("len before: ", len_before_target)


        if len_after_target >= params.neighborhood_length:
            find_downstream = False

        if len_before_target >= params.neighborhood_length:
            find_upstream = False

        #print(find_downstream)
        #print(find_upstream)



        downstream_paths_file = {"file_name" : "" }
        upstream_paths_file = {"file_name" : "" }
       	target_seq = directed_graph.nodes[target_gene[0]]['sequence']
        if(len(target_gene)>1):
        	for gene_index in range(1, len(target_gene)):
        		target_seq = target_seq + directed_graph.nodes[target_gene[gene_index]]['sequence'][directed_graph[target_gene[gene_index-1]][target_gene[gene_index]]['overlap']:]

        # target_seq = target_seq[len_before_target:-len_after_target]
        # downstream

        if find_downstream:
            print("into down")
            target_gene_node = target_gene[-1]

            loop_thread = threading.Thread(target=get_paths_from_big_nx_graph_4, args=(directed_graph, target_gene_node, len_after_target, "down", params, stop_flag_downstream, downstream_paths_file))
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
            #    directed_graph, target_gene_node, len_after_target, "down", params, stop_flag)
            print("done downstream")

        # upstream
        if find_upstream:
            print("into up")
            target_gene_node = target_gene[0]
            #stop_flag.clear()
            print("flag: ", stop_flag_upstream.is_set())
            loop_thread = threading.Thread(target=get_paths_from_big_nx_graph_4, args=(reverse_directed_graph, target_gene_node, len_before_target, "up", params, stop_flag_upstream, upstream_paths_file))
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
            #    reverse_directed_graph, target_gene_node, len_before_target, "up", params, stop_flag)
            print("done upstream")

        print(f"downstream file :{downstream_paths_file}")
        print(f"upstream file : {upstream_paths_file}")
       	merge_upstream_target_downstream_5(
            		upstream_paths_file=upstream_paths_file["file_name"], downstream_paths_file=downstream_paths_file["file_name"], target=target_gene, target_name=target_name,
            		target_seq=target_seq, len_before_target=len_before_target, len_after_target=len_after_target, params=params, seq_file=seq_file)


    path_info_list = create_paths_info_list(seq_file, directed_graph, target_hits, params.neighborhood_length, params.max_kmer_size, params.assembler)

    write_paths_info_to_file(path_info_list, paths_info_file)

    return seq_file, paths_info_file



def create_paths_info_list(seq_file, ego_graph, target_hits, threshold, max_kmer_size, assembler):
    """
    Build per-node path/coverage info for every extracted neighborhood sequence.

    For each path it walks the upstream nodes, the target-gene node(s), and the
    downstream nodes, recording each node's coverage and its start/end offset
    within the assembled sequence.
    Parameters:
        seq_file: the file of extracted neighborhood sequences (path -> sequence)
        ego_graph: the assembly graph the nodes belong to
        target_hits: the alignment hits used to locate the target gene in a path
        threshold: the neighborhood length
        max_kmer_size: the assembler's max k-mer size (for coverage calculation)
        assembler: the assembler name (for coverage calculation)
    Return:
        a flat list of per-node info dicts (sequence, node, coverage, start, end)
    """
    paths_info_list = []
    paths = read_paths_file(seq_file)
    counter = 0
    for path, seq in paths.items():
        seq = seq.upper()
        end = -1
        #print(f"path : {path} , {type(path)}, target_gene, {target_gene}, {type(target_gene)}")

        for index, element in enumerate(path):
            if isinstance(element, tuple):
                #print(f"The first tuple is at index: {index}")
                index_target_in_path = index
                for target_hit in target_hits:
                    #print(target_hit)
                    # target_gene
                    target_gene = [pair[0]+pair[1]
                        for pair in zip(target_hit['nodes'], target_hit['orientations'])]
                    if(tuple(target_gene)==element):
                        if(target_hit['end_pos'] == 0):
                            len_after_target = 0
                        else:
                            len_after_target = len(
                                ego_graph.nodes[target_gene[-1]]['sequence']) - target_hit['end_pos']

                        if(target_hit['start_pos'] == 0):
                            len_before_target = 0
                        else:
                            len_before_target = target_hit['start_pos'] - 1
                break

        print(f"len_befor : {len_before_target}")
        print(f"len after : {len_after_target}")
        print(f"target : {target_gene}")
        #index_target_in_path = path.index(target_gene)
        #print("index target : ", index_target_in_path)

        is_upstreams = (index_target_in_path > 0)
        is_downstreams = (len(path)-1) > index_target_in_path

        #print(is_upstreams)
        #print(is_downstreams)

        #### upstreams
        if(is_upstreams):
            for up_path_index in range(index_target_in_path):
                node = path[up_path_index]
                coverage = calculate_coverage(ego_graph.nodes[path[up_path_index]], max_kmer_size, path[up_path_index], assembler)
                start = end + 1
                if(type(path[up_path_index+1]) == tuple):
                    #print(path[up_path_index + 1])
                    #print("path[up_path_index + 1][0] : ", path[up_path_index + 1][0])
                    end = seq.find(ego_graph.nodes[path[up_path_index + 1][0]]['sequence'][0:len_before_target]) - 1
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

        #print(len(path[index_target_in_path]))
        ##### target
        if(len(path[index_target_in_path]) == 1):
            print(f"len_befor : {len_before_target}")
            print(f"len after : {len_after_target}")
            node = path[index_target_in_path]
            coverage = calculate_coverage(ego_graph.nodes[path[index_target_in_path][0]], max_kmer_size, path[index_target_in_path][0], assembler)
            ####upstream
            if(len_before_target > 0):
                start = end + 1
                end = (end + len_before_target) if(len_before_target <= threshold) else (end + threshold)
                #print("info : ", info)
                info = {"sequence" : counter,
                        "node" : node,
                        "coverage" : coverage,
                        "start" : start,
                        "end" : end}
                paths_info_list.append(info)

            ####target
            start = end + 1
            end = end + len(ego_graph.nodes[path[index_target_in_path][0]]['sequence'][len_before_target:]) if (len_after_target==0) else end + len(ego_graph.nodes[path[index_target_in_path][0]]['sequence'][len_before_target:- len_after_target])
            info = {"sequence" : counter,
                    "node" : node,
                    "coverage" : coverage,
                    "start" : start,
                    "end" : end}
            #print("info : ", info)
            paths_info_list.append(info)

            ###downstream
            if(len_after_target>0):
                start = end + 1
                end = end + threshold  if (len_after_target >= threshold) else end + len_after_target
                info = {"sequence" : counter,
                        "node" : node,
                        "coverage" : coverage,
                        "start" : start,
                        "end" : end}
                #print("info : ", info)
                paths_info_list.append(info)

        elif len(path[index_target_in_path]) == 2 :
            ##### upstream
            node = path[index_target_in_path][0]
            coverage = calculate_coverage(ego_graph.nodes[path[index_target_in_path][0]], max_kmer_size, path[index_target_in_path][0], assembler)
            if(len_before_target > 0):
                start = end + 1
                end = (end + len_before_target) if(len_before_target <= threshold) else (end + threshold)
                info = {"sequence" : counter,
                        "node" : node,
                        "coverage" : coverage,
                        "start" : start,
                        "end" : end}
                paths_info_list.append(info)

            #####target
            start = end + 1
            end = end + len(ego_graph.nodes[path[index_target_in_path][0]]['sequence'][len_before_target:])
            info = {"sequence" : counter,
                    "node" : node,
                    "coverage" : coverage,
                    "start" : start,
                    "end" : end}
            paths_info_list.append(info)

            node = path[index_target_in_path][-1]
            coverage = calculate_coverage(ego_graph.nodes[path[index_target_in_path][-1]], max_kmer_size, path[index_target_in_path][-1], assembler)


            start = end + 1
            end = end + len(ego_graph.nodes[path[index_target_in_path][-1]]['sequence']) - len_after_target - int(ego_graph[path[index_target_in_path][-2]][path[index_target_in_path][-1]]['overlap'])
            info = {"sequence" : counter,
                    "node" : node,
                    "coverage" : coverage,
                    "start" : start,
                    "end" : end}
            paths_info_list.append(info)

            #### downstream
            if(len_after_target>0):
                start = end + 1
                if(index_target_in_path < len(path) - 1):
                    end = (end + len_after_target) if(len_after_target <= threshold) else (end + threshold)
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
            node = path[index_target_in_path][0]
            coverage = calculate_coverage(ego_graph.nodes[path[index_target_in_path][0]], max_kmer_size, path[index_target_in_path][0], assembler)
            if(len_before_target > 0):
                start = end + 1
                end = (end + len_before_target) if(len_before_target <= threshold) else (end + threshold)
                info = {"sequence" : counter,
                        "node" : node,
                        "coverage" : coverage,
                        "start" : start,
                        "end" : end}
                paths_info_list.append(info)

            #####target
            start = end + 1
            end = end + len(ego_graph.nodes[node]['sequence'][len_before_target:])
            info = {"sequence" : counter,
                    "node" : node,
                    "coverage" : coverage,
                    "start" : start,
                    "end" : end}
            paths_info_list.append(info)


            for target_index in range(1, len(path[index_target_in_path])-1) :
                node = path[index_target_in_path][target_index]
                coverage = calculate_coverage(ego_graph.nodes[path[index_target_in_path][target_index]], max_kmer_size, path[index_target_in_path][target_index], assembler)
                start = end + 1
                end = end + ego_graph[path[index_target_in_path][target_index-1]][path[index_target_in_path][target_index]]['weight']
                info = {"sequence" : counter,
                    "node" : node,
                    "coverage" : coverage,
                    "start" : start,
                    "end" : end}
                paths_info_list.append(info)

            node = path[index_target_in_path][-1]
            coverage = calculate_coverage(ego_graph.nodes[path[index_target_in_path][-1]], max_kmer_size, path[index_target_in_path][-1], assembler)

            start = end + 1
            end = end + len(ego_graph.nodes[path[index_target_in_path][-1]]['sequence']) - len_after_target - int(ego_graph[path[index_target_in_path][-2]][path[index_target_in_path][-1]]['overlap'])
            info = {"sequence" : counter,
                    "node" : node,
                    "coverage" : coverage,
                    "start" : start,
                    "end" : end}
            paths_info_list.append(info)

            #### downstream
            if(len_after_target>0):
                start = end + 1
                if(index_target_in_path < len(path) - 1):
                    end = (end + len_after_target) if(len_after_target <= threshold) else (end + threshold)
                else:
                    end = len(seq) - 1
                info = {"sequence" : counter,
                        "node" : node,
                        "coverage" : coverage,
                        "start" : start,
                        "end" : end}
                paths_info_list.append(info)

        #### downstreams
        if(is_downstreams):
            #print("index_target_in_path : ", index_target_in_path)
            #print("len path :", len(path))

            for down_path_index in range(index_target_in_path+1, len(path)-1):
                #print(down_path_index)
                node = path[down_path_index]
                coverage = calculate_coverage(ego_graph.nodes[path[down_path_index]], max_kmer_size, path[down_path_index], assembler)
                start = end + 1
                if(down_path_index==index_target_in_path+1):
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

# Worker-process globals for neighborhood extraction. The assembly graph is
# large, so rather than binding it into a per-task callable (which re-pickles the
# whole graph for every target), it is handed to each worker once via the pool
# initializer and read from these globals during extraction.
_WORKER_DIRECTED_GRAPH = None
_WORKER_REVERSE_GRAPH = None
_WORKER_PARAMS = None


def _init_extraction_worker(directed_graph, reverse_directed_graph, params):
    """Pool initializer: store the graph in each worker process once."""
    global _WORKER_DIRECTED_GRAPH, _WORKER_REVERSE_GRAPH, _WORKER_PARAMS
    _WORKER_DIRECTED_GRAPH = directed_graph
    _WORKER_REVERSE_GRAPH = reverse_directed_graph
    _WORKER_PARAMS = params


def _extract_in_worker(target_info):
    """Per-task worker: only the small per-target payload is pickled per call."""
    return neighborhood_sequence_extraction(
        _WORKER_DIRECTED_GRAPH, _WORKER_REVERSE_GRAPH, _WORKER_PARAMS, target_info
    )


def extract_target_neighborhoods(
        params,
        gfa_file,
        target_seq_align_info,
        debug: bool
):
    """
    Extract the neighborhood sequences for every target gene found in the graph.

    Builds the directed (and reversed) assembly graph from the GFA file once, then
    runs neighborhood extraction for each target (in a multiprocessing pool unless
    single-threaded), and cleans up the scratch files afterwards.
    Parameters:
        params: the parsed CLI parameters
        gfa_file: the assembly graph in GFA format
        target_seq_align_info: list of (target_file_path, target_hits) pairs
        debug: kept for interface compatibility
    Return:
        (seq_files, path_info_files): the per-target sequence and path-info files
    """
    output_dir = Path(params.output_dir)
    sequence_dir = output_dir / SEQ_DIR_NAME
    sequence_dir.mkdir(parents=True, exist_ok=True)

    (output_dir / "final_result").mkdir(parents=True, exist_ok=True)
    (output_dir / "final_down_up").mkdir(parents=True, exist_ok=True)
    # Read GFA file
    gfa_graph = gfapy.Gfa.from_file(gfa_file)
    directed_graph = create_directed_graph_nx(gfa_graph)
    del gfa_graph  # no longer needed; free the gfapy structure before extraction

    reverse_directed_graph = directed_graph.reverse()
    for edge in reverse_directed_graph.edges:
        reverse_directed_graph[edge[0]][edge[1]]['weight'] = len(
            reverse_directed_graph.nodes[edge[1]]['sequence']) - reverse_directed_graph[edge[0]][edge[1]]['overlap']

    # If running single threaded do not add any overhead using multiprocessing pool
    if params.num_cores == 1:
        lists = list()
        for x in target_seq_align_info:
            lists.append(neighborhood_sequence_extraction(directed_graph, reverse_directed_graph, params, x))
    else:
        # The workers only traverse (never mutate) the graph, so share a single
        # physical copy via copy-on-write fork inheritance rather than giving each
        # worker its own. Under "fork" the initargs are inherited, not pickled;
        # gc.freeze() stops the cyclic collector from touching (and thus copying)
        # the shared pages in the children. Fall back to the default context
        # (which pickles one copy per worker) where fork is unavailable.
        try:
            ctx = multiprocessing.get_context("fork")
        except ValueError:
            ctx = multiprocessing.get_context()
        gc.freeze()
        try:
            with ctx.Pool(
                params.num_cores,
                initializer=_init_extraction_worker,
                initargs=(directed_graph, reverse_directed_graph, params),
            ) as p:
                lists = list(p.imap(_extract_in_worker, target_seq_align_info))
        finally:
            gc.unfreeze()
    seq_files, path_info_files = zip(*lists)
    LOG.info("delete files")
    for item_path in Path(params.output_dir).iterdir():
        # Check if it is a file, and if so, delete it
        if item_path.is_file():
            if item_path.suffix not in ('.log', '.txt'):
                item_path.unlink()



    shutil.rmtree(f"{params.output_dir}/final_result")
    shutil.rmtree(f"{params.output_dir}/final_down_up")

    return seq_files, path_info_files
