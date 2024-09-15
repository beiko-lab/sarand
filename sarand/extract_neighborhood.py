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
import csv
import subprocess
from gfapy.sequence import rc
import shutil
import multiprocessing
from csv import DictReader
from pathlib import Path

from sarand.util.logger import LOG
from sarand.utils import (
    reverse_sign,
    find_node_name,
    exist_in_path,
    compare_two_sequences,
    read_path_info_from_align_file,
    delete_lines_started_with
)

from sarand.config import SEQ_DIR_NAME, SEQ_NAME_PREFIX

OUT_DIR = "output"
TEMP_DIR = "temp"


"""
AM: This is never used.
"""
# def retrieve_AMR(file_path):
#     """
#     To read the AMR gene from the fasta file.
#     Parameters:
#             file_path:	the address of the file containing the AMR gene
#     Return:
#             the sequence of the AMR gene in lower case
#     """
#     with open(file_path) as fp:
#         for line in fp:
#             if line.startswith(">"):
#                 continue
#             return line.lower()


def similar_sequence_exits(seq_list, query_seq, output_dir):
    """
    To check if any similar sequence to query_seq exists in seq_list
    if it finds any similar one with a shorter length, replaces that with query_Seq
    Parameters:
            seq_list: 	the list of sequences
            query_seq:	the query sequence
            output_dir:	the output directory
    Return:
            the index of the similar seq if it has a shorter length than query_seq (to be replace)
            and
            True if a seq similar to query_seq was found in seq_list

    """
    for i, seq in enumerate(seq_list):
        if compare_two_sequences(seq, query_seq, output_dir, blast_ext=str(i)):
            # if the new sequence has a longer length replace the shorter similar oe with it
            if len(seq) < len(query_seq):
                return i, True
            else:
                return -1, True
    return -1, False


def write_sequences_to_file(sequence_list, path_list, file_name):
    """
    To write the extracted sequences and paths in neighborhood of the AMR gene in a file
    Parameters:
            sequence_list: the list of sequences extracted in AMR gene neighborhood
            path_list: the list of all paths (a list of nodes) in AMR gene neighborhood
            file_name: the name of file to be writtn in
    """
    file = open(file_name, "a+")
    for seq, path in zip(sequence_list, path_list):
        myLine = "> " + ", ".join(path) + ":"
        file.write(myLine + "\n")
        file.write(seq + "\n")
    file.close()


def sequence_on_orientation(seq, orient):
    """
    To	return the sequence based on the orientation used in the graph
    Parameters:
            orient:	+ / -
            seq:	the original sequence
    Return: reverse complement for -; otherwise seq itself
            for orient = '-' and seq = 'AACTG' -> 'CAGTT'
    """
    if orient == "-":
        return rc(seq)
    else:
        return seq


def find_overlap(
    from_node, to_node, from_orient, to_orient, validate_overlap_size=True
):
    """
    To return the size of overlap beween  from_node and to_node
    e.g.; from_node = 2, to_node = 3, from_orient = +, to_orient= +, edge = [L 2 + 3 + 55M]
                            ---->  return 55
    Parameters:
            from_node:	the first node of the edge
            to_node:	the second node of the edge
            from_orient:the orientation of the first node of the edge
            to_orient:	the orientation of the second node of the edge
            validate_overlap_size: to verify if the extracted size for overlap is correct and matches
    Return:
            the size of overlap
    """
    size = 0
    for edge in from_node.dovetails_L + from_node.dovetails_R:
        # 3561+ --> 15082+
        if (
            edge.from_segment == from_node
            and edge.to_segment == to_node
            and edge.from_orient == from_orient
            and edge.to_orient == to_orient
        ):
            if edge.overlap[0].code == "M":
                size = edge.overlap[0].length
                if size > 0 and validate_overlap_size:
                    seq1 = sequence_on_orientation(
                        edge.from_segment.sequence, edge.from_orient
                    )
                    seq2 = sequence_on_orientation(
                        edge.to_segment.sequence, edge.to_orient
                    )
                    if seq1[len(seq1) - size :] != seq2[:size]:
                        LOG.info("WARNING: the overlap doesn't match: " + str(edge))
                        LOG.info("first sequence: " + seq1)
                        LOG.info("second sequence: " + seq2)
            else:
                LOG.info(
                    "WARNING: the sequence from node "
                    + from_node.name
                    + " to node "
                    + to_node.name
                    + "has "
                    + edge.overlap[0]
                    + " instead of overlap!"
                )
            return size
        # sometimes we have the reverse complement of the edge and not itself
        # 15082- --> 3561-
        elif (
            edge.from_segment == to_node
            and edge.to_segment == from_node
            and edge.from_orient == reverse_sign(to_orient)
            and edge.to_orient == reverse_sign(from_orient)
        ):
            if edge.overlap[0].code == "M":
                size = edge.overlap[0].length
                if validate_overlap_size:
                    seq1 = sequence_on_orientation(
                        edge.to_segment.sequence, reverse_sign(edge.to_orient)
                    )
                    seq2 = sequence_on_orientation(
                        edge.from_segment.sequence, reverse_sign(edge.from_orient)
                    )
                    if seq1[len(seq1) - size :] != seq2[:size]:
                        LOG.info("WARNING: the overlap doesn't match: " + str(edge))
                        LOG.info("first sequence: " + seq1)
                        LOG.info("second sequence: " + seq2)
            else:
                LOG.info(
                    "WARNING: the sequence from node "
                    + from_node.name
                    + " to node "
                    + to_node.name
                    + "has "
                    + edge.overlap[0]
                    + " instead of overlap!"
                )
            return size
    return size


"""
AM: This is never used.
"""
# def remove_identical_sub_sequences(sequence_list, path_list):
#     """
#     To keep only one instance in case of multiple identical extracted sequences or
#     when a sequence is a subsequence of another one
#     Parameters:
#             sequence_list:	list of extracted sequences
#             path_list:		list of extracted paths
#     Return sequence_list and path_list in which every sequence/path is unique and
#             none of them is subssequence/subpath of another one
#     """
#     # In case two sequences are identical keep the one with the lower index in the list
#     identical_sub_set = set()
#     for i, seq1 in enumerate(sequence_list):
#         for j, seq2 in enumerate(sequence_list):
#             if i < j:
#                 if seq1 == seq2 or seq2 in seq1:
#                     identical_sub_set.add(j)
#                 elif seq1 in seq2:
#                     identical_sub_set.add(i)
#
#     identical_sub_list = sorted(identical_sub_set)
#     for index in reversed(identical_sub_list):
#         del sequence_list[index]
#         del path_list[index]
#
#     return sequence_list, path_list


def reverse_path(path):
    """
    To generate the reverse of a path;
            e.g., [1+, [8-, 12+], 9-] --> [9+, [12-, 8+,] 1-]
    """
    mypath = []
    items = ["[", "]", "{", "}"]
    reversed_items = ["]", "[", "}", "{"]
    for node in reversed(path):
        mynode = ""
        num = ""
        for ch in node:
            if ch in items:
                mynode = reversed_items[items.index(ch)] + mynode
            elif ch == "-" or ch == "+":
                num += reverse_sign(ch)
                mynode = num + mynode
            else:
                num += ch
        mypath.append(mynode)
    return mypath


def extract_found_amr(myGraph, node_list, orientation_list, start_pos, end_pos):
    """
    To extract the sequence representing AMR gene from the path found by
    Bandage+BLAST presented in node_list; I could have used the original AMR gene
    sequence but I thought there might be cases that the assembler can't assemble
    the entire AMR gene or there are some errors in assembled one.
    Parameters:
            myGraph: The original graph
            node_list:	the list of nodes containing the AMR gene (i.e., path)
            orientation_list: list of orientation of nodes in the path
            start_pos:	the position of AMR in the first node from which the AMR gene has started
            end_pos: the position of AMR in the last node in which the AMR gene ended
    Return:
                    the sequence of the AMR gene extracted from the path presented by node_list
    """
    if not node_list:
        LOG.error("ERROR: There is no node in node_list representing the AMR gene!")
        import pdb

        pdb.set_trace()
        sys.exit()

    if len(node_list) == 1:
        return sequence_on_orientation(
            myGraph.segment(node_list[0]).sequence, orientation_list[0]
        )[start_pos - 1 : end_pos]

    # Add first node
    found_amr_seq = sequence_on_orientation(
        myGraph.segment(node_list[0]).sequence, orientation_list[0]
    )[start_pos - 1 :]
    # Add middle nodes
    if len(node_list) > 2:
        for i in range(len(node_list) - 2):
            seq = sequence_on_orientation(
                myGraph.segment(node_list[i + 1]).sequence, orientation_list[i + 1]
            )
            overlap_size = find_overlap(
                myGraph.segment(node_list[i]),
                myGraph.segment(node_list[i + 1]),
                orientation_list[i],
                orientation_list[i + 1],
            )
            seq = seq[overlap_size:]
            found_amr_seq += seq
    # Add last node
    seq = sequence_on_orientation(
        myGraph.segment(node_list[-1]).sequence, orientation_list[-1]
    )[:end_pos]
    overlap_size = find_overlap(
        myGraph.segment(node_list[-2]),
        myGraph.segment(node_list[-1]),
        orientation_list[-2],
        orientation_list[-1],
    )
    # if only a part of overlap is in the AMR sequence:
    # e.g., with overlap_size = 55, (93) 26719919-, 5498944+ (32)
    if end_pos < overlap_size:
        len_not_in_amr = overlap_size - end_pos
        found_amr_seq = found_amr_seq[:-len_not_in_amr]
    else:
        seq = seq[overlap_size:]
        found_amr_seq += seq

    return found_amr_seq


def generate_sequence_path(
    myGraph,
    node_list,
    orientation_list,
    start_pos,
    end_pos,
    pre_sequence_list,
    post_sequence_list,
    pre_path_list,
    post_path_list,
    path_length_list,
    output_name,
    path_dir="same",
):
    """
    To concatenate  the sequences and paths before and after the AMR sequence.
    also, for any new sequence we add it if no similar enough sequence exists in the
    list. If such a sequence already exits in the list and its length is shorter,
    we replace it with the new sequence.
    e.g., pre_path_list = [['1-','2+']], post_path_list=[['7+'],['8+','9+','10-']],
    node_dir_list = ['[13+','14-]']
    - if path_dir = 'same' --> 		[
                                            ['1-','2+','[13+','14-]','7+'],
                                            ['1-','2+','[13+','14-]','8+','9+','10-']
                                                                    ]
    - if path_dir = 'reverse' --> 	[
                                            ['7-', '[14+', '13-]', '2-', '1+'],
                                            ['10+', '9-', '8-', '[14+', '13-]', '2-', '1+']
                                                                    ]
    Parameters:
            myGraph: The original graph
            node_list:	the list of nodes containing the AMR gene (i.e., path)
            orientation_list: list of orientation of nodes in the path
            start_pos:	the position of AMR in the first node from which the AMR gene has started
            end_pos: the position of AMR in the last node in which the AMR gene ended
            pre_sequence_list: the list of sequences preceding the AMR sequence in the graph
            post_sequence_list: the list of sequences following the AMR sequence in the graph
            pre_path_list: the list of nodes (and their orientation) preceding the AMR node
            post_path_list: the list of nodes (and their orientation) following the AMR node
            path_dir: 'same' if the AMR path direction (node_dir_list) is the same as what
                    blast returned (presenting actual AMR sequence) and 'reverse' if the AMR path
                    direction (node_dir_list) is the reverse of blast result (presenting the
                    reverse complement of AMR sequence); in the latter, we need to return the
                    reverse complement of the concatenated sequence to have the actual AMR
                    sequence and not its reverse.

    """
    # extract found AMR sequence
    found_amr_seq = extract_found_amr(
        myGraph, node_list, orientation_list, start_pos, end_pos
    )
    # create a temporaty directory
    with tempfile.TemporaryDirectory() as temp_dir:

        # to be able to go over both 'for' loops make lists non-empty
        if not pre_path_list:
            pre_path_list = [[]]
        if not post_path_list:
            post_path_list = [[]]
        if not pre_sequence_list:
            pre_sequence_list = [""]
        if not post_sequence_list:
            post_sequence_list = [""]
        # add direction to nodes representing AMR gene
        node_orient_list = [
            node + orient for node, orient in zip(node_list, orientation_list)
        ]
        # annotate the path representing AMR gene
        node_orient_list[0] = "[" + node_orient_list[0]
        node_orient_list[-1] = node_orient_list[-1] + "]"
        # Generating all found sequences	and paths
        path_list = []
        sequence_list = []
        counter = 0
        path_info_list = []
        for pre_seq, pre_path in zip(pre_sequence_list, pre_path_list):
            for post_seq, post_path in zip(post_sequence_list, post_path_list):
                if path_dir == "same":
                    sequence = pre_seq.upper() + found_amr_seq.lower() + post_seq.upper()
                    path = pre_path + node_orient_list + post_path
                elif path_dir == "reverse":
                    sequence = rc(
                        pre_seq.upper() + found_amr_seq.lower() + post_seq.upper()
                    )
                    path = reverse_path(pre_path + node_orient_list + post_path)
                index, found = similar_sequence_exits(sequence_list, sequence, temp_dir)
                if not found:
                    sequence_list.append(sequence)
                    path_list.append(path)
                    path_info_list.append(path_length_list[counter])
                elif index >= 0:
                    sequence_list[index] = sequence
                    path_list[index] = path
                    path_info_list[index] = path_length_list[counter]
                counter += 1

    return sequence_list, path_list, path_info_list


def append_path_sequence(
    sequence,
    path,
    sequence_list,
    path_list,
    output_dir,
    path_length,
    path_length_list,
    outfile,
):
    """ """
    index, found = similar_sequence_exits(sequence_list, sequence, output_dir)
    if not found:
        sequence_list.append(sequence)
        path_list.append(path)
        path_length_list.append(path_length)
        # write outputs to the temp file
        if sequence != "":
            with open(outfile, "a") as fd:
                writer = csv.writer(fd)
                writer.writerow([len(sequence_list) - 1, sequence, path, path_length])

    elif index >= 0:
        sequence_list[index] = sequence
        path_list[index] = path
        path_length_list[index] = path_length
        # change corresponding entry in temp file
        if sequence != "":
            r = csv.reader(open(outfile))
            lines = list(r)
            # lines[0] contains the headers
            lines[index + 1] = [str(index), str(sequence), str(path), str(path_length)]
            writer = csv.writer(open(outfile, "w"))
            writer.writerows(lines)

    return sequence_list, path_list, path_length_list


def extract_post_sequence_recursively_both_dir(
    node,
    node_orient,
    current_seq,
    current_path,
    length,
    sequence_list,
    path_list,
    node_list,
    compare_dir,
    path_thr,
    current_path_length,
    path_length_list,
    out_file,
):
    """
    To extract recursively the sequences following the AMR gene sequence and their paths
    In this function, we do this:
    - Find node's immediate neighbors like B
    - if there is an edge node --> B and the orient of node in that is the same as node_orient,
            B is eligible to be considered in the post_path for AMR sequence.
    - if there is an edge B --> node and the orient of node in that is reverse of
            node_orient, B with thereverse orient is eligible to be considered in the
            post_path for AMR sequence.
    Parameters:
            node:		the staring node in next paths
            node_orient:the orientation of the node in the edge this function instance was called from
            current_seq:the seq found so far
            current_path: the path found so far
            length: the remained length of the sequences around the AMR gene to be extracted
            sequence_list: the list of sequences following the AMR gene sequence
            path_list:	the list of nodes (and their orientation) following the nodes presenting the AMR gene
            node_list:	the list of nodes presenting the AMR gene
            path_thr: the threshold used for recursive pre_path and post_path
                    search as long as the length of the path is less that this threshold
    Return:
            modified sequence_list and path_list
    """
    seq = ""
    found_any_edge = False
    if length > 0 and len(current_path) < path_thr:
        for edge in node.dovetails_L + node.dovetails_R:
            # the second part of condition ensures that we are finding only nodes that are not presenting the AMR gene itself
            # Logically that shouldn't be the case unless there is a loop in the network:
            # path A -> B -> C represents the AMR; however we have an edge C -> A
            # so when we are looking for the nodes following AMR , A is selected too!!!!!
            to_segment = None
            to_orient = ""
            if (
                edge.from_segment.name == node.name
                and edge.to_segment.name != node.name
                and edge.to_segment.name not in node_list
                and edge.from_orient == node_orient
            ):
                to_segment = edge.to_segment
                to_orient = edge.to_orient
            elif (
                edge.to_segment.name == node.name
                and edge.from_segment.name != node.name
                and edge.from_segment.name not in node_list
                and edge.to_orient == reverse_sign(node_orient)
            ):
                to_segment = edge.from_segment
                to_orient = reverse_sign(edge.from_orient)
            if to_segment and to_orient != "":
                # taking care of loops in the path
                index = exist_in_path(current_path, to_segment.name + to_orient)
                if index == -1:
                    found_any_edge = True
                    new_seq = sequence_on_orientation(to_segment.sequence, to_orient)
                    # Remove the overlap between nodes' sequences
                    overlap_size = find_overlap(
                        edge.from_segment,
                        edge.to_segment,
                        edge.from_orient,
                        edge.to_orient,
                    )
                    new_seq = new_seq[overlap_size:]
                    seq = current_seq + new_seq
                    path = list.copy(current_path)
                    path.append(str(to_segment.name) + str(to_orient))
                    path_length = list.copy(current_path_length)
                    path_length.append(min(len(new_seq), length))
                    if len(new_seq) >= length:
                        (
                            sequence_list,
                            path_list,
                            path_length_list,
                        ) = append_path_sequence(
                            seq[: len(current_seq) + length],
                            path,
                            sequence_list,
                            path_list,
                            compare_dir,
                            path_length,
                            path_length_list,
                            out_file,
                        )
                    else:
                        (
                            sequence_list,
                            path_list,
                            path_length_list,
                        ) = extract_post_sequence_recursively_both_dir(
                            to_segment,
                            to_orient,
                            seq,
                            path,
                            length - len(new_seq),
                            sequence_list,
                            path_list,
                            node_list,
                            compare_dir,
                            path_thr,
                            path_length,
                            path_length_list,
                            out_file,
                        )
                # graph loops are specified in { }
                elif index > -1:
                    current_path[index] = "{" + current_path[index]
                    current_path[-1] += "}"
    if not found_any_edge:
        sequence_list, path_list, path_length_list = append_path_sequence(
            current_seq,
            current_path,
            sequence_list,
            path_list,
            compare_dir,
            current_path_length,
            path_length_list,
            out_file,
        )
    return sequence_list, path_list, path_length_list


def extract_post_sequence(
    node,
    node_orient,
    node_list,
    length,
    end_pos,
    output_dir,
    output_name,
    path_thr,
    time_out_counter,
):
    """
    The initial function to extract the post_sequence of AMR.
    it adds the sequence following the AMR sequence (if there is any) in the last
    node in AMR path and if more is required, calls a recursive function to go over neighbors.
    Parameters:
            node:		the staring node in next paths
            node_orient:the orientation of the node in AMR path
            node_list:	the list of nodes presenting the AMR gene
            length: 	the length of the sequences around the AMR gene to be extracted
            end_pos: 	the position of AMR in the last node in which the AMR gene ended
            path_thr: the threshold used for recursive pre_path and post_path
                    search as long as the length of the path is less that this threshold
    Return:
            the lists of sequences and paths following AMR
    """
    post_sequence_list = []
    post_path_list = []
    path_length = []
    path_length_list = []
    post_sequence = ""
    postfix_length = length
    if end_pos > 0:
        # attach the end of the AMR node (not included in AMR seq) to post_sequence
        postfix_length -= len(node.sequence) - end_pos
        post_sequence = sequence_on_orientation(node.sequence, node_orient)[end_pos:]
    # Find the sequence after the AMR gene
    if postfix_length <= 0:
        post_sequence = post_sequence[:length]
        post_sequence_list.append(post_sequence)
        path_length.append(len(post_sequence))
        path_length_list.append(path_length)
    # check all edges started from last_segment
    else:
        # create a temporaty directory
        compare_dir = os.path.join(output_dir, "temp_comparison_" + output_name)
        os.makedirs(compare_dir, exist_ok=True)
        # File to store pre_sequence results in case we need to stop the task by time-out before it completes
        temp_result = os.path.join(
            output_dir,
            "post_temp_result_"
            + output_name
            + "_"
            + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            + ".csv",
        )
        with open(temp_result, "w") as fd:
            writer = csv.writer(fd)
            writer.writerow(["index", "post_seq", "post_path", "post_path_len"])
        if len(post_sequence) > 0:
            path_length.append(len(post_sequence))
        p = ""
        if time_out_counter <= 0:
            extract_post_sequence_recursively_both_dir(
                node=node,
                node_orient=node_orient,
                current_seq=post_sequence,
                current_path=[],
                length=postfix_length,
                sequence_list=[],
                path_list=[],
                node_list=node_list,
                compare_dir=compare_dir,
                path_thr=path_thr,
                current_path_length=path_length,
                path_length_list=path_length_list,
                out_file=temp_result,
            )
        else:
            p = multiprocessing.Process(
                target=extract_post_sequence_recursively_both_dir,
                args=(
                    node,
                    node_orient,
                    post_sequence,
                    [],
                    postfix_length,
                    [],
                    [],
                    node_list,
                    compare_dir,
                    path_thr,
                    path_length,
                    path_length_list,
                    temp_result,
                ),
            )
            p.start()
            # Wait for time_out_counter seconds or until process finishes
            p.join(time_out_counter)

        if p != "" and p.is_alive():
            LOG.info("extract_post_sequence is still running... let's kill it...")
            # Terminate - may not work if process is stuck for good
            p.terminate()
            # OR Kill - will work for sure, no chance for process to finish nicely however
            # p.kill()
            p.join()
        # Read the result from temp file
        with open(temp_result) as fd:
            myreader = DictReader(fd)
            for row in myreader:
                post_sequence_list.append(row["post_seq"])
                if row["post_path"] != "[]":
                    path = [a.strip("'") for a in row["post_path"][1:-1].split(", ")]
                else:
                    path = []
                post_path_list.append(path)
                path_length = [
                    int(l.strip())
                    for l in row["post_path_len"][1:-1].split(",")
                    if l.strip() != ""
                ]
                path_length_list.append(path_length)
        # delete temp file
        if os.path.isfile(temp_result):
            os.remove(temp_result)
        # delete the temp folder
        if os.path.exists(compare_dir):
            try:
                shutil.rmtree(compare_dir)
            except OSError as e:
                LOG.error("Error: %s - %s." % (e.filename, e.strerror))

    return post_sequence_list, post_path_list, path_length_list


def extract_pre_sequence_recursively_both_dir(
    node,
    node_orient,
    current_seq,
    current_path,
    length,
    sequence_list,
    path_list,
    node_list,
    compare_dir,
    path_thr,
    current_path_length,
    path_length_list,
    out_file,
):
    """
    To extract recursively the sequences preceding the AMR gene sequence and their paths
    In this function, we do this:
    - Find node's immediate neighbors like B
    - if there is an edge B --> node and the orient of node in that is the same as node_orient,
            B is eligible to be considered in the pre_path for AMR sequence.
    - if there is an edge node --> B and the orient of node in that is reverse of
            node_orient, B with the reverse orient is eligible to be considered in the
            pre_path for AMR sequence.
    Parameters:
            node:		the staring node in next paths
            node_orient:the orientation of the node in the edge this function instance was called from
            current_seq:the seq found so far
            current_path: the path found so far
            length: the remained length of the sequences around the AMR gene to be extracted
            sequence_list: the list of sequences preceding the AMR gene sequence
            path_list:	the list of nodes (and their orientation) preceding the nodes presenting the AMR gene
            node_list:	the list of nodes presenting the AMR gene
            path_thr: the threshold used for recursive pre_path and post_path
                    search as long as the length of the path is less that this threshold
    Return:
            modified sequence_list and path_list
    """
    seq = ""
    found_any_edge = False
    if length > 0 and len(current_path) < path_thr:
        for edge in node.dovetails_R + node.dovetails_L:
            # the second part ofcondition ensures that we are finding only nodes that are not presenting the AMR gene itself
            # Logically that shouldn't be the case unless there is a loop in the network:
            # path A -> B -> C represents the AMR; however we have an edge C -> A
            # so when we are looking for the nodes preceding AMR , C is selected too!!!!!
            from_segment = None
            from_orient = ""
            if (
                edge.to_segment.name == node.name
                and edge.from_segment.name != node.name
                and edge.from_segment.name not in node_list
                and edge.to_orient == node_orient
            ):
                from_segment = edge.from_segment
                from_orient = edge.from_orient
            elif (
                edge.from_segment.name == node.name
                and edge.to_segment.name != node.name
                and edge.to_segment.name not in node_list
                and edge.from_orient == reverse_sign(node_orient)
            ):
                from_segment = edge.to_segment
                from_orient = reverse_sign(edge.to_orient)
            if from_segment and from_orient != "":
                # taking care of loops in the path
                index = exist_in_path(current_path, from_segment.name + from_orient)
                if index == -1:
                    found_any_edge = True
                    new_seq = sequence_on_orientation(
                        from_segment.sequence, from_orient
                    )
                    # Remove the overlap between nodes' sequences
                    overlap_size = find_overlap(
                        edge.from_segment,
                        edge.to_segment,
                        edge.from_orient,
                        edge.to_orient,
                    )
                    new_seq = new_seq[: len(new_seq) - overlap_size]
                    seq = new_seq + current_seq
                    path = list.copy(current_path)
                    path.insert(0, str(from_segment.name) + str(from_orient))
                    path_length = list.copy(current_path_length)
                    path_length.insert(0, min(len(new_seq), length))
                    if len(new_seq) >= length:
                        (
                            sequence_list,
                            path_list,
                            path_length_list,
                        ) = append_path_sequence(
                            seq[len(new_seq) - length :],
                            path,
                            sequence_list,
                            path_list,
                            compare_dir,
                            path_length,
                            path_length_list,
                            out_file,
                        )
                    else:
                        (
                            sequence_list,
                            path_list,
                            path_length_list,
                        ) = extract_pre_sequence_recursively_both_dir(
                            from_segment,
                            from_orient,
                            seq,
                            path,
                            length - len(new_seq),
                            sequence_list,
                            path_list,
                            node_list,
                            compare_dir,
                            path_thr,
                            path_length,
                            path_length_list,
                            out_file,
                        )
                # graph loops are specified in { }
                elif index > -1:
                    current_path[index] = current_path[index] + "}"
                    current_path[0] = "{" + current_path[0]
    if not found_any_edge:
        sequence_list, path_list, path_length_list = append_path_sequence(
            current_seq,
            current_path,
            sequence_list,
            path_list,
            compare_dir,
            current_path_length,
            path_length_list,
            out_file,
        )
    return sequence_list, path_list, path_length_list


def extract_pre_sequence(
    node,
    node_orient,
    node_list,
    length,
    start_pos,
    output_dir,
    output_name,
    path_thr,
    time_out_counter,
):
    """
    The initial function to extract the pre_sequence of AMR.
    it adds the sequence preceding the AMR sequence (if there is any) in the first
    node in AMR path and if more sequence is required, calls a recursive function to go over neighbors.
    Parameters:
            node:		the staring node in next paths
            node_orient:the orientation of the node in AMR path
            node_list:	the list of nodes presenting the AMR gene
            length: 	the length of the sequences around the AMR gene to be extracted
            end_pos: 	the position of AMR in the last node in which the AMR gene ended
            path_thr: the threshold used for recursive pre_path and post_path
                    search as long as the length of the path is less that this threshold
    Return:
            the lists of sequences and paths preceding AMR

    """
    pre_sequence_list = []
    pre_path_list = []
    pre_sequence = ""
    path_length = []
    path_length_list = []
    prefix_length = length
    if start_pos > 0:
        # attach the beginning of the AMR node (not included in AMR seq) to pre_sequence
        prefix_length -= start_pos - 1
        pre_sequence = sequence_on_orientation(node.sequence, node_orient)[
            : start_pos - 1
        ]
    # Find the sequence before the AMR gene
    if prefix_length <= 0:
        pre_sequence = pre_sequence[-length:]
        pre_sequence_list.append(pre_sequence)
        path_length.append(len(pre_sequence))
        path_length_list.append(path_length)
    # check all edges ended at first node
    else:
        # create a temporaty directory
        compare_dir = os.path.join(output_dir, "temp_comparison_" + output_name)
        os.makedirs(compare_dir, exist_ok=True)

        # File to store pre_sequence results in case we need to stop the task by time-out before it completes
        temp_result = os.path.join(
            output_dir,
            "pre_temp_result_"
            + output_name
            + "_"
            + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            + ".csv",
        )
        with open(temp_result, "w") as fd:
            writer = csv.writer(fd)
            writer.writerow(["index", "pre_seq", "pre_path", "pre_path_len"])
        if len(pre_sequence) > 0:
            path_length.append(len(pre_sequence))
        p = ""
        if time_out_counter <= 0:
            extract_pre_sequence_recursively_both_dir(
                node=node,
                node_orient=node_orient,
                current_seq=pre_sequence,
                current_path=[],
                length=prefix_length,
                sequence_list=[],
                path_list=[],
                node_list=node_list,
                compare_dir=compare_dir,
                path_thr=path_thr,
                current_path_length=path_length,
                path_length_list=path_length_list,
                out_file=temp_result,
            )
        else:
            p = multiprocessing.Process(
                target=extract_pre_sequence_recursively_both_dir,
                args=(
                    node,
                    node_orient,
                    pre_sequence,
                    [],
                    prefix_length,
                    [],
                    [],
                    node_list,
                    compare_dir,
                    path_thr,
                    path_length,
                    path_length_list,
                    temp_result,
                ),
            )
            p.start()
            # Wait for time_out_counter seconds or until process finishes
            p.join(time_out_counter)
        if p != "" and p.is_alive():
            LOG.info("extract_pre_sequence is still running... let's kill it...")
            # Terminate - may not work if process is stuck for good
            p.terminate()
            # OR Kill - will work for sure, no chance for process to finish nicely however
            # p.kill()
            p.join()
        # Read the result from temp file
        with open(temp_result) as fd:
            myreader = DictReader(fd)
            for row in myreader:
                pre_sequence_list.append(row["pre_seq"])
                if row["pre_path"] != "[]":
                    path = [a.strip("'") for a in row["pre_path"][1:-1].split(", ")]
                else:
                    path = []
                pre_path_list.append(path)
                path_length = [
                    int(l.strip())
                    for l in row["pre_path_len"][1:-1].split(",")
                    if l.strip() != ""
                ]
                path_length_list.append(path_length)
        # delete temp file
        if os.path.isfile(temp_result):
            os.remove(temp_result)
        # delete the temp folder
        if os.path.exists(compare_dir):
            try:
                shutil.rmtree(compare_dir)
            except OSError as e:
                LOG.error("Error: %s - %s." % (e.filename, e.strerror))

    return pre_sequence_list, pre_path_list, path_length_list


def write_paths_info_to_file(paths_info_list, paths_info_file, seq_counter):
    """ """
    counter = seq_counter
    with open(paths_info_file, "a") as fd:
        writer = csv.writer(fd)
        for i, path_info in enumerate(paths_info_list):
            for node_info in path_info:
                counter = i + seq_counter + 1
                writer.writerow(
                    [
                        counter,
                        node_info["node"],
                        node_info["coverage"],
                        node_info["start"],
                        node_info["end"],
                    ]
                )
    return counter


def calculate_coverage(segment, max_kmer_size, node, assembler):
    """
    To calculate node coverage based on assembler type
    Parameters:
            segment/node: the node whose covearge is supposed to eb calculated
            max_kmer_size: the max_kmer_size used by assembler
    Return:
            coverage

    """
    if assembler == "metaspades" or assembler == "bcalm":
        coverage = segment.KC / (len(segment.sequence) - max_kmer_size)
    elif assembler == "megahit":
        coverage = float(node.split("cov_")[1].split("_")[0])
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


def generate_node_range_coverage(
    myGraph,
    node_list,
    orientation_list,
    start_pos,
    end_pos,
    pre_path_list,
    pre_path_length_list,
    post_path_list,
    post_path_length_list,
    max_kmer_size,
    assembler="metaspades",
):
    """
    To generate a list of all nodes representing a sequence their start and end
    in the sequence and their coverage
    Parameters:
            myGraph: the assembly graph
            node_list: the list of nodes reprsenting AMR
            orientation_list: the orientation of nodes representing AMR
            start_pos: the start position of sequence in the first node in node_list
            end_pos: the end position of sequence in the last node in node_list
            pre_path_list: the node details for nodes representing upstream
            pre_path_length_list: the node length list for nodes representing upstream
            post_path_list: the node details for nodes representing downstream
            post_path_length_list: the node length list for nodes representing downstream
            max_kmer_size: the maximum kmer-size used by assembler
            assembler: the assembler used to generate the graph
    Results:
            the list of nodes and their start and end in the sequence for every extracted sequence
    """
    if not pre_path_list:
        pre_path_list = [[]]
    if not post_path_list:
        post_path_list = [[]]
    if not pre_path_length_list:
        pre_path_length_list = [[]]
    if not post_path_length_list:
        post_path_length_list = [[]]
    # find amr_path_length_info
    start = 0
    amr_path_length_info = []
    for i, node in enumerate(node_list):
        segment = myGraph.segment(node)
        coverage = calculate_coverage(segment, max_kmer_size, node, assembler)
        # the length of first node in amr
        if i == 0 and len(node_list) == 1:
            length = end_pos - start_pos + 1
        elif i == 0 and len(node_list) > 1:
            length = len(segment.sequence) - start_pos + 1
        # the length of intermediate nodes in amr
        elif i < (len(node_list) - 1):
            overlap_size = find_overlap(
                myGraph.segment(node_list[i - 1]),
                segment,
                orientation_list[i - 1],
                orientation_list[i],
            )
            length = len(segment.sequence) - overlap_size
        # the length of last node in amr
        else:
            overlap_size = find_overlap(
                myGraph.segment(node_list[i - 1]),
                segment,
                orientation_list[i - 1],
                orientation_list[i],
            )
            length = end_pos - overlap_size
        node_info = {
            "node": node,
            "coverage": coverage,
            "start": start,
            "end": start + length - 1,
        }
        amr_path_length_info.append(node_info)
        start = start + length

    # pre_path_length_info
    pre_paths_info = []
    for path, path_length in zip(pre_path_list, pre_path_length_list):
        start = 0
        pre_path_info = []
        assert (
            len(path) == len(path_length) or len(path) == len(path_length) - 1
        ), "inconsistent length of arrays: path vs path_length"
        difference = len(path_length) - len(path)
        end_common = len(path_length) - difference
        for node, length in zip(path, path_length[:end_common]):
            pure_node = find_node_name(node)
            segment = myGraph.segment(pure_node)
            coverage = calculate_coverage(segment, max_kmer_size, pure_node, assembler)
            node_info = {
                "node": pure_node,
                "coverage": coverage,
                "start": start,
                "end": start + length - 1,
            }
            start = start + length
            pre_path_info.append(node_info)
        if len(path) == len(path_length) - 1:
            segment = myGraph.segment(node_list[0])
            coverage = calculate_coverage(
                segment, max_kmer_size, node_list[0], assembler
            )
            node_info_last = {
                "node": node_list[0],
                "coverage": coverage,
                "start": start,
                "end": start + path_length[-1] - 1,
            }
            pre_path_info.append(node_info_last)
        pre_paths_info.append(pre_path_info)

    # attach amr info to pre_path
    for i in range(len(pre_paths_info)):
        if len(pre_paths_info[i]) > 0:
            lag = pre_paths_info[i][-1]["end"] + 1
        else:
            lag = 0
        for amr_info in amr_path_length_info:
            tmp = amr_info.copy()
            tmp["start"] += lag
            tmp["end"] += lag
            pre_paths_info[i].append(tmp)

    # post_path_length_info
    post_paths_info = []
    for path, path_length in zip(post_path_list, post_path_length_list):
        start = 0
        post_path_info = []
        assert (
            len(path) == len(path_length) or len(path) == len(path_length) - 1
        ), "inconsistent length of arrays: path vs path_length"
        start_index = len(path_length) - len(path)
        if len(path) == len(path_length) - 1:
            segment = myGraph.segment(node_list[-1])
            coverage = calculate_coverage(
                segment, max_kmer_size, node_list[-1], assembler
            )
            node_info = {
                "node": node_list[-1],
                "coverage": coverage,
                "start": start,
                "end": start + path_length[0] - 1,
            }
            post_path_info.append(node_info)
            start = start + path_length[0]
        for node, length in zip(path, path_length[start_index:]):
            pure_node = find_node_name(node)
            segment = myGraph.segment(pure_node)
            coverage = calculate_coverage(segment, max_kmer_size, pure_node, assembler)
            node_info = {
                "node": pure_node,
                "coverage": coverage,
                "start": start,
                "end": start + length - 1,
            }
            start = start + length
            post_path_info.append(node_info)
        post_paths_info.append(post_path_info)

    # generate the entire path info by adding post_path
    paths_info = []
    for pre_path in pre_paths_info:
        if len(pre_path) > 0:
            lag = pre_path[-1]["end"] + 1
        else:
            lag = 0
        for post_path in post_paths_info:
            path_info = list.copy(pre_path)
            for node_info in post_path:
                tmp = node_info.copy()
                tmp["start"] += lag
                tmp["end"] += lag
                path_info.append(tmp)
            paths_info.append(path_info)

    return paths_info


def extract_neighborhood_sequence(
    gfa_file,
    length,
    amr_path_info,
    path_thr,
    output_name,
    max_kmer_size,
    seq_info,
    output_dir,
    time_out_counter,
    assembler="metaspades",
):
    """
    To extract all linear sequences with length = length from the start and end of the AMR gene
    In some cases, we might need to extract only the upstream or downstream
    Parameters:
            myGraph: The assembly graph
            length: the length of all sequences around the AMR gene to be extracted
            amr_path_info: contains node_list (the list of nodes containing the AMR gene (i.e., path)),
                    their orientation, start_pos (the position of AMR in the first node from which
                    the AMR gene has started) and end_pos (the position of AMR in the last node
                    in which the AMR gene ended)
            path_thr: the threshold used for recursive pre_path and post_path
                    search as long as the length of the path is less that this threshold
            output_name: used for naming files
            max_kmer_size: the maximum kmer size used by the assembler
            seq_info: a dictionary to store upstream (presequence), downstream
                    (postsequence) and their length
            assembler: the assembler used to assemble the graph
    Results:
            the list of extracted sequences and paths (the nodes representing the sequence)
                    and their info as well as modified seq_info
    """
    try:
        LOG.debug(f"Loading the graph from {gfa_file}...")
        myGraph = gfapy.Gfa.from_file(gfa_file)
    except Exception as e:
        LOG.error("Graph not loaded successfully: " + str(e))
        if assembler == "metacherchant":
            return [], [], [], []
        else:
            import pdb

            pdb.set_trace()
    node_list = amr_path_info["nodes"]
    orientation_list = amr_path_info["orientations"]
    start_pos = amr_path_info["start_pos"]
    end_pos = amr_path_info["end_pos"]

    last_segment = myGraph.segment(node_list[-1])
    if start_pos == 0:
        start_pos = 1
    if end_pos == 0:
        end_pos = len(str(last_segment.sequence))

    LOG.debug(
        "last_segment = "
        + last_segment.name
        + " start_pos = "
        + str(start_pos)
        + " end_pos= "
        + str(end_pos)
    )

    # Find the sequences and paths for the path provided by bandage+blast for AMR gene
    if not seq_info["pre_seq"]:
        # Find the sequence before the AMR gene
        LOG.debug("Running extract_pre_sequence " + output_name + " ...")
        pre_sequence_list, pre_path_list, pre_path_length_list = extract_pre_sequence(
            myGraph.segment(node_list[0]),
            orientation_list[0],
            node_list,
            length,
            start_pos,
            output_dir,
            output_name,
            path_thr,
            time_out_counter,
        )
        seq_info["pre_seq"] = pre_sequence_list
        seq_info["pre_path"] = pre_path_list
        seq_info["pre_len"] = pre_path_length_list
    else:
        # up_stream has already found for another path with the same start node and position
        pre_sequence_list = seq_info["pre_seq"]
        pre_path_list = seq_info["pre_path"]
        pre_path_length_list = seq_info["pre_len"]

    if not seq_info["post_seq"]:
        # Find the sequence after the AMR gene
        LOG.debug("Running extract_post_sequence " + output_name + " ...")
        (
            post_sequence_list,
            post_path_list,
            post_path_length_list,
        ) = extract_post_sequence(
            last_segment,
            orientation_list[-1],
            node_list,
            length,
            end_pos,
            output_dir,
            output_name,
            path_thr,
            time_out_counter,
        )
        seq_info["post_seq"] = post_sequence_list
        seq_info["post_path"] = post_path_list
        seq_info["post_len"] = post_path_length_list
    else:
        # down_stream has already found for another path with the same end node and position
        post_sequence_list = seq_info["post_seq"]
        post_path_list = seq_info["post_path"]
        post_path_length_list = seq_info["post_len"]

    # combine path_info from pre and post
    LOG.debug("Running generate_node_range_coverage " + output_name + " ...")
    path_length_list = generate_node_range_coverage(
        myGraph,
        node_list,
        orientation_list,
        start_pos,
        end_pos,
        pre_path_list,
        pre_path_length_list,
        post_path_list,
        post_path_length_list,
        max_kmer_size,
        assembler,
    )
    # Combine pre_ and post_ sequences and paths
    LOG.debug("Running generate_sequence_path " + output_name + " ...")
    sequence_list, path_list, path_info_list = generate_sequence_path(
        myGraph,
        node_list,
        orientation_list,
        start_pos,
        end_pos,
        pre_sequence_list,
        post_sequence_list,
        pre_path_list,
        post_path_list,
        path_length_list,
        output_name,
        "same",
    )

    return sequence_list, path_list, path_info_list, seq_info


def find_amr_related_nodes(
    amr_file, gfa_file, output_dir, threshold=95, output_pre="", align_file=""
):
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
            threshold: threshold for identity and coverage
            output_pre: used for naming output file
    Return:
            A boolean value which is True if any path was found and
            A list of dictionaries each denoting an AMR path
    """
    if align_file == "":
        # Run bandage+blast
        output_name = os.path.join(
            output_dir,
            output_pre
            + "_align_"
            + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"),
        )
        if os.path.isfile(output_name + ".tsv"):
            os.remove(output_name + ".tsv")
        bandage_command = subprocess.run(
            [
				"Bandage",
                "querypaths",
                gfa_file,
                amr_file,
                output_name,
                "--pathnodes",
                "50",
                "--minpatcov",
                str((threshold - 1) / 100.0),
                "--minmeanid",
                str((threshold - 1) / 100.0),
                "--minhitcov",
                str((threshold - 1) / 100.0),
            ],
            stdout=subprocess.PIPE,
            check=True,
        )
        LOG.info(bandage_command.stdout.decode("utf-8"))

        align_file = output_name + ".tsv"
    # Process the output tsv file
    found, paths_info = read_path_info_from_align_file(output_name + ".tsv", threshold)

    return found, paths_info


def check_if_similar_ng_extractions_exist(
    amr_path_info, amr_paths_info, assembler="metaspades"
):
    """
    To check the new AMR path found by alignment with the ones already processed!
    If there is a similar AMR path in amr_paths_info, we don't need to process
    amr_path_info!
    Parameters:
            amr_path_info: the new AMR path
            amr_paths_info: the list of AMR paths processed so far
    Results:
            a dictionary containing the index of an entry from amr_paths_info with similar
            start node-position or end node-position
    """
    found_up_stream = False
    found_down_stream = False
    similar_path = {"up_stream": -1, "down_stream": -1}
    new_first_node = amr_path_info["nodes"][0] + amr_path_info["orientations"][0]
    new_last_node = amr_path_info["nodes"][-1] + amr_path_info["orientations"][-1]
    new_start_pos = amr_path_info["start_pos"]
    new_end_pos = amr_path_info["end_pos"]
    # first look for a path whose both first and last nodes are the same as the new path
    # If not found a path with similar first node and similar last node
    # then look for a path that its first_node=new_first node
    # and a path that its last_node=new_last_node
    for i, path_info in enumerate(amr_paths_info):
        start_pos = path_info["start_pos"]
        end_pos = path_info["end_pos"]
        if start_pos != new_start_pos and end_pos != new_end_pos:
            continue
        first_node = path_info["nodes"][0] + path_info["orientations"][0]
        last_node = path_info["nodes"][-1] + path_info["orientations"][-1]
        if (
            len(path_info) == 1
            and len(amr_path_info) == 1
            and new_first_node == first_node
        ):
            LOG.error(
                "there shouldnt be two separate amr paths with the same single node!"
            )
            import pdb

            pdb.set_trace()
        if (
            new_first_node == first_node
            and new_start_pos == start_pos
            and new_last_node == last_node
            and new_end_pos == end_pos
        ):
            return {"up_stream": i, "down_stream": i}
        if (
            not found_up_stream
            and new_first_node == first_node
            and new_start_pos == start_pos
        ):
            similar_path["up_stream"] = i
            found_up_stream = True
        if (
            not found_down_stream
            and new_last_node == last_node
            and new_end_pos == end_pos
        ):
            similar_path["down_stream"] = i
            found_down_stream = True

    if found_up_stream and found_down_stream:
        return similar_path
    # there are some cases that a path starts with one more node than another path
    # E.g. P1 = {a1, a2} P2 = {a0, a1, a2}
    # Or a path ends at one more node than another path
    # E.g., P1= {a1, a2} P2 = {a1, a2, a3}
    new_second_node = ""
    new_second_last_node = ""
    if len(amr_path_info["nodes"]) > 2:
        new_second_node = amr_path_info["nodes"][1] + amr_path_info["orientations"][1]
        new_second_last_node = (
            amr_path_info["nodes"][-2] + amr_path_info["orientations"][-2]
        )
        for i, path_info in enumerate(amr_paths_info):
            first_node = path_info["nodes"][0] + path_info["orientations"][0]
            last_node = path_info["nodes"][-1] + path_info["orientations"][-1]
            second_node = ""
            second_last_node = ""
            if len(path_info["nodes"]) > 2:
                second_node = path_info["nodes"][1] + path_info["orientations"][1]
                second_last_node = (
                    path_info["nodes"][-2] + path_info["orientations"][-2]
                )
            if not found_up_stream and (
                new_first_node == second_node or new_second_node == first_node
            ):
                similar_path["up_stream"] = i
                found_up_stream = True
            if not found_down_stream and (
                new_last_node == second_last_node or new_second_last_node == last_node
            ):
                similar_path["down_stream"] = i
                found_down_stream = True
            if found_up_stream and found_down_stream:
                return similar_path

    return similar_path


"""
AM: This method is not used anywhere.
"""
# def order_path_nodes(path_nodes, amr_file, out_dir, threshold=90):
#     """
#     Given that we have a list of nodes that are supposed to represent a given AMR
#     (i.e., AMR path), this method returns their right order to represent the AMR
#     sequence.
#     Curretly, this method is only used for Metacherchant
#     Parameters:
#             path_nodes: the list of nodes representing AMR
#             amr_file: the file containing AMR sequence
#             out_dir: the dir to store some temporary files
#             threshold: used for identity and coverage in bandage+blast
#     Returns:
#             the lists of sorted nodes and their corresponding orientation
#     """
#     path_nodes_info = []
#     no_path_nodes = []
#     for i, node in enumerate(path_nodes):
#         node_info_list = []
#         # write the sequence into a fasta file
#         query_file = create_fasta_file(
#             node.sequence, out_dir, comment=">" + node.name + "\n", file_name="query"
#         )
#         # run blast query for alignement
#         blast_file_name = os.path.join(out_dir, "blast.csv")
#         blast_file = open(blast_file_name, "w")
#         subprocess.run(
#             [
#                 "blastn",
#                 "-query",
#                 query_file,
#                 "-subject",
#                 amr_file,
#                 "-task",
#                 "blastn",
#                 "-outfmt",
#                 "10",
#                 "-max_target_seqs",
#                 "10",
#                 "-evalue",
#                 "0.5",
#                 "-perc_identity",
#                 str(threshold - 1),
#             ],
#             stdout=blast_file,
#             check=True,
#         )
#         blast_file.close()
#         # command = 'blastn -query '+query_file+' -subject '+amr_file+\
#         # 	' -task blastn -outfmt 10 -max_target_seqs 10 -evalue 0.5 -perc_identity '+\
#         # 	str(threshold)+' > '+ blast_file_name
#         # os.system(command)
#         with open(blast_file_name, "r") as file1:
#             myfile = csv.reader(file1)
#             for row in myfile:
#                 identity = int(float(row[2]))
#                 coverage = int(float(row[3]) / len(node.sequence) * 100)
#                 if identity >= threshold and coverage >= threshold:
#                     node_info = {
#                         "name": node.name,
#                         "c_start": int(row[8]),
#                         "c_end": int(row[9]),
#                     }
#                     node_info_list.append(node_info)
#         if not node_info_list:
#             no_path_nodes.append(i)
#             LOG.error(node.name + " was not found in the sequence of " + amr_file)
#         path_nodes_info.append(node_info_list)
#     # order nodes
#     start_list = []
#     orientations = []
#     for path_info in path_nodes_info:
#         if len(path_info) > 0:
#             if path_info[0]["c_start"] < path_info[0]["c_end"]:
#                 start_list.append(path_info[0]["c_start"])
#                 orientations.append("+")
#             else:
#                 start_list.append(path_info[0]["c_end"])
#                 orientations.append("-")
#         else:
#             start_list.append(-1)
#             orientations.append("/")
#     # start_list =[e[0]['c_start'] if len(e)>0 else -1 for e in path_nodes_info ]
#     sorted_indeces = sorted(range(len(start_list)), key=lambda k: start_list[k])
#     sorted_path_nodes = []
#     sorted_orientations = []
#     for index in sorted_indeces:
#         if index not in no_path_nodes:
#             sorted_path_nodes.append(path_nodes[index])
#             sorted_orientations.append(orientations[index])
#
#     return sorted_path_nodes, sorted_orientations

"""
AM: This method is not used anywhere.
"""
# def extract_amr_align_from_file(gfa_file):
#     """
#     Retrieve the list of segments that represent AMR path
#     This method is use for Metacherchant in which such nodes ends with '_start'
#     Parameters:
#             gfa_file: the assembly graph file
#     Return:
#             the list of graph segments representing the AMR path
#     """
#     myGraph = gfapy.Gfa.from_file(gfa_file)
#     path_nodes = []
#     for segment in myGraph.segments:
#         if segment.name.endswith("_start"):
#             path_nodes.append(segment)
#     return path_nodes


def neighborhood_sequence_extraction(
    gfa_file,
    length,
    output_dir,
    threshold=95,
    seq_name_prefix="ng_sequences_",
    # output_name = 'ng_sequences',
    path_node_threshold=10,
    max_kmer_size=55,
    time_out_counter=-1,
    assembler="metaspades",
    amr_seq_align_info="",
):
    """
    The core function to extract the sequences/paths preceding and following the AMR sequence
    in the assembly graph
    Parameters:
            gfa_file:	the GFA file containing the assembly graph
            length:		the length of all sequences around the AMR gene to be extracted
            output_dir: the output directory to store files
            threshold: 	threshold for identity and coverage
            seq_name_prefix: used for naming output file
            path_node_threshold: the threshold used for recursive pre_path and post_path
                    search as long as the length of the path is less that this threshold
            max_kmer_size: the maximum kmer used by the assembler
            assembler: the assembler used to generate the assembly graph
            amr_seq_align_info: a tuple containing the name of file containing amr sequence
                    and the amr alignment info
    Return:
            the name of file containing the list of extracted sequences/paths
    """
    # @Somayeh: what is going on here? amr_seq_align_info isn't passed to this
	# @Fin: You kept only the parallel case here and the last argument is passed when
	# the method is called: lists = p.map(p_extraction, amr_seq_align_info)
    # function
    amr_file, amr_paths_info = amr_seq_align_info
    LOG.debug("amr_file = " + amr_file)
    output_name = seq_name_prefix + os.path.splitext(os.path.basename(amr_file))[0]
    seq_output_dir = os.path.join(output_dir, "sequences")
    os.makedirs(seq_output_dir, exist_ok=True)

    seq_file = os.path.join(
        seq_output_dir,
        output_name
        + "_"
        + str(length)
        + "_"
        + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
        + ".txt",
    )

    if not amr_paths_info:
        os.makedirs(os.path.join(output_dir, "alignment_files"), exist_ok=True)
        # remained_len_thr = length - (length*path_seq_len_percent_threshod/100.0)
        # find all AMR paths in the assembly graph
        """
        AM: I'm certain this branch is impossible to reach.
        """
        found, amr_paths_info = find_amr_related_nodes(
            amr_file,
            gfa_file,
            os.path.join(output_dir, "alignment_files"),
            threshold,
            output_name,
            ""
        )
        if not found:
            LOG.error("no alignment was found for " + amr_file)
            # import pdb; pdb.set_trace()
            return "", ""

    # csv file for path_info
    os.makedirs(os.path.join(output_dir, "paths_info"), exist_ok=True)
    paths_info_file = os.path.join(
        output_dir,
        "paths_info",
        output_name
        + "_"
        + str(length)
        + "_"
        + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
        + ".csv",
    )
    with open(paths_info_file, "a") as fd:
        writer = csv.writer(fd)
        writer.writerow(["sequence", "node", "coverage", "start", "end"])

    # Extract the sequenc of AMR neighborhood
    LOG.debug(f"Calling extract_neighborhood_sequence for {os.path.basename(amr_file)}...")
    seq_counter = 0
    sequence_lists = []
    #path_lists = []
    #path_info_lists = []
    checked_amr_paths_info = []
    seq_info_list = []
    for i, amr_path_info in enumerate(amr_paths_info):
        seq_info = {
            "pre_seq": None,
            "pre_path": None,
            "pre_len": None,
            "post_seq": None,
            "post_path": None,
            "post_len": None,
        }
        similar_path = check_if_similar_ng_extractions_exist(
            amr_path_info, checked_amr_paths_info
        )
        # LOG.info('amr_path_info: '+str(amr_path_info)+' checked_amr_paths_info: '+str(checked_amr_paths_info)+' similar_path'+str(similar_path))
        if (
            similar_path["up_stream"] > -1
            and similar_path["up_stream"] == similar_path["down_stream"]
        ):
            # no need to save any information for this path
            continue
        elif similar_path["up_stream"] == -1 and similar_path["down_stream"] == -1:
            (
                sequence_list,
                path_list,
                path_info_list,
                seq_info,
            ) = extract_neighborhood_sequence(
                gfa_file,
                length,
                amr_path_info,
                path_node_threshold,
                output_name,
                max_kmer_size,
                seq_info,
                output_dir,
                time_out_counter,
                assembler,
            )
        elif similar_path["up_stream"] > -1 and similar_path["down_stream"] > -1:
            # we have extracted both its upstream and downstream but in different paths
            u_id = similar_path["up_stream"]
            d_id = similar_path["down_stream"]
            seq_info = {
                "pre_seq": seq_info_list[u_id]["pre_seq"],
                "pre_path": seq_info_list[u_id]["pre_path"],
                "pre_len": seq_info_list[u_id]["pre_len"],
                "post_seq": seq_info_list[d_id]["post_seq"],
                "post_path": seq_info_list[d_id]["post_path"],
                "post_len": seq_info_list[d_id]["post_len"],
            }
            (
                sequence_list,
                path_list,
                path_info_list,
                seq_info,
            ) = extract_neighborhood_sequence(
                gfa_file,
                length,
                amr_path_info,
                path_node_threshold,
                output_name,
                max_kmer_size,
                seq_info,
                output_dir,
                time_out_counter,
                assembler,
            )
        elif similar_path["up_stream"] > -1:
            # we need to extract its downstream
            u_id = similar_path["up_stream"]
            seq_info["pre_seq"] = seq_info_list[u_id]["pre_seq"]
            seq_info["pre_path"] = seq_info_list[u_id]["pre_path"]
            seq_info["pre_len"] = seq_info_list[u_id]["pre_len"]
            (
                sequence_list,
                path_list,
                path_info_list,
                seq_info,
            ) = extract_neighborhood_sequence(
                gfa_file,
                length,
                amr_path_info,
                path_node_threshold,
                output_name,
                max_kmer_size,
                seq_info,
                output_dir,
                time_out_counter,
                assembler,
            )

        elif similar_path["down_stream"] > -1:
            # we need to extract its upstream
            d_id = similar_path["down_stream"]
            seq_info["post_seq"] = seq_info_list[d_id]["post_seq"]
            seq_info["post_path"] = seq_info_list[d_id]["post_path"]
            seq_info["post_len"] = seq_info_list[d_id]["post_len"]
            (
                sequence_list,
                path_list,
                path_info_list,
                seq_info,
            ) = extract_neighborhood_sequence(
                gfa_file,
                length,
                amr_path_info,
                path_node_threshold,
                output_name,
                max_kmer_size,
                seq_info,
                output_dir,
                time_out_counter,
                assembler,
            )
        if not sequence_list:
            continue
        seq_info_list.append(seq_info)
        sequence_lists.append(sequence_list)
        # path_lists.append(path_list)
        # @Somayeh: this list is never used?
		# @Fin: you're right!
        #path_info_lists.append(path_info_list)
        checked_amr_paths_info.append(amr_path_info)
        # if not sequence_list:
        # 	continue
        write_sequences_to_file(sequence_list, path_list, seq_file)
        LOG.info(
            "Extracted neighborhoods written to " + seq_file
        )
        seq_counter = write_paths_info_to_file(
            path_info_list, paths_info_file, seq_counter
        )
    return seq_file, paths_info_file

def sequence_neighborhood_main(
        params,
        gfa_file,
        amr_seq_align_info,
        debug: bool
):
    """
    The core function to extract the neighborhood of AMRs
    Parameters:
        params: the list pf parameters imported from params.py
		bandage_path: the path to bandage executable file
        gfa_file: the file containing the assembly graph
        amr_seq_align_info: the alignment info (AMR alignment in the graph)
    Return:
        the list of files containing extracted sequence and the details of nodes representing them
    """

    # seq_files = []
    # path_info_files = []
    LOG.info(f"Extracting neighborhood sequences with length = {params.neighbourhood_length}")

    # remove paths from GFA file
    # AM: We don't want to modify the users input data so this now goes to a new file
    """
    AM: Since we don't want to modify the users input data, this is now being written
    to a new file in the output directory. This could be in a temporary directory
    instead.
    """
    path_gfa_new = Path(params.output_dir) / f'{gfa_file.stem}_no_paths.gfa'
    delete_lines_started_with("P", gfa_file, path_gfa_new)
    gfa_file = path_gfa_new

    sequence_dir = os.path.join(
        params.output_dir,
        SEQ_DIR_NAME,
        f'{SEQ_DIR_NAME}_{params.neighbourhood_length}',
    )
    os.makedirs(sequence_dir, exist_ok=True)

    # If running single threaded do not add any overhead using multiprocessing pool
    if params.num_cores == 1:
        lists = list()
        for x in amr_seq_align_info:
            lists.append(neighborhood_sequence_extraction(gfa_file,
                                                          params.neighbourhood_length,
                                                          sequence_dir,
                                                          params.min_target_identity,
                                                          SEQ_NAME_PREFIX,
                                                          1000,  # should this really be an option? path_node_threshold
                                                          params.max_kmer_size,
                                                          params.extraction_timeout,
                                                          params.assembler, x))
    else:
        p_extraction = partial(
            neighborhood_sequence_extraction,
            gfa_file,
            params.neighbourhood_length,
            sequence_dir,
            params.min_target_identity,
            SEQ_NAME_PREFIX,
            1000,  # should this really be an option? path_node_threshold
            params.max_kmer_size,
            params.extraction_timeout,
            params.assembler
        )
        with Pool(params.num_cores) as p:
            lists = list(p.imap(p_extraction, amr_seq_align_info))
    seq_files, path_info_files = zip(*lists)

    # AM: To clean up the file we created earlier
    os.remove(path_gfa_new)

    if debug:
        try_dump_to_disk(
            {
                'seq_files': seq_files,
                'path_info_files': path_info_files
            },
            Path(sequence_dir) / 'debug_sequence_neighborhood_main.json'
        )

    return seq_files, path_info_files


if __name__ == "__main__":
    pass

# text = 'This code is used to find the context of a given AMR gene'
# parser = argparse.ArgumentParser(description=text)
# parser.add_argument('--amr','-A', type=str, default = '',
# 	help = 'the path of the file containing the AMR gene sequence')
# parser.add_argument('--gfa', '-G', type=str, default='',
# 	help = 'the GFA file containing the graph')
# parser.add_argument('--length', '-L', type = int, default=1000,
# 	help = 'the length of AMR gene\'s neighbourhood to be extracted')
# parser.add_argument('--output_dir','-O', type=str, default = OUT_DIR,
# 	help = 'the output directory to store the results')
# args = parser.parse_args()

# if not os.path.exists(args.output_dir):
# 	os.makedirs(args.output_dir)

# if not args.amr:
# 	logging.error('please enter the path for the AMR gene file')
# 	import pdb; pdb.set_trace()
# 	sys.exit()

# if not args.gfa:
# 	logging.error('please enter the path for the assembly file')
# 	import pdb; pdb.set_trace()
# 	sys.exit()

# neighborhood_sequence_extraction(args.gfa, args.length, args.output_dir,
# 							amr_seq_align_info = (arg.amr,''))
