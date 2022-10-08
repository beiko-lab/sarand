"""
File:		utils.py
Author:		Somayeh Kafaie
Date:		March 2021
Purpose:	containing all utility functions used in other files
"""

import sys
import os
import errno
import argparse
import logging
import re
import csv
import datetime
import collections
import shutil
import subprocess
from pathlib import Path


def extract_name_from_file_name(file_name):
    """ """
    return os.path.splitext(os.path.basename(file_name))[0]


def amr_name_from_comment(amr_comment):
    """ """
    amr_name = (
        amr_comment.split("[")[0]
        .split("|")[-1]
        .strip()
        .replace(" ", "_")
        .replace("'", ";")
        .replace("/", "]")
    )
    return amr_name


def amr_name_from_title(amr_title):
    """ """
    return amr_title.strip().replace(" ", "_").replace("'", ";").replace("/", "]")


def restricted_amr_name_from_modified_name(amr_name):
    """ """
    amr_name1 = amr_name.replace(";", "SS")
    amr_name1 = "".join(
        e for e in amr_name1 if e.isalpha() or e.isnumeric() or e == "_" or e == "-"
    )
    return amr_name1


def retreive_original_amr_name(amr_name):
    """ """
    return amr_name.replace(";", "'").replace("]", "/")


def create_fasta_file(seq, output_dir, comment="> sequence:\n", file_name="temp"):
    """
    To create a fasta file for a sequence
    Parameters:
        seq: the sequence to be written into the file
        output_dir: the output directory address
        comment: the comment to be written into fasta file
        file_name: the name of the fasta file
    Return:
        the address of the fasta file
    """
    myfile_name = os.path.join(output_dir, file_name + ".fasta")
    if os.path.isfile(myfile_name):
        os.remove(myfile_name)
    myfile = open(myfile_name, "w")
    myfile.write(comment)
    if not comment.endswith("\n"):
        myfile.write("\n")
    myfile.write(seq)
    if not seq.endswith("\n"):
        myfile.write("\n")
    myfile.close()
    return myfile_name


def initialize_logger(log_file_path):
    """
    Initialise the log handler using a specific log_file_path
    """
    logging.basicConfig(level=logging.DEBUG)
    # logging.basicConfig()
    log_formatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]	%(message)s")
    root_logger = logging.getLogger()
    file_handler = logging.FileHandler(log_file_path)
    file_handler.setFormatter(log_formatter)
    root_logger.addHandler(file_handler)
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setFormatter(log_formatter)
    root_logger.addHandler(console_handler)


def retrieve_AMR(file_path):
    """
    To read the AMR gene from the text file.
    Parameters:
        file_path:	the address of the file containing the AMR gene
    Return:
        the sequence of the AMR gene in lower case
    """
    amr_name = ""
    with open(file_path) as fp:
        for i, line in enumerate(fp):
            # skip comment line
            if line.startswith(">"):
                amr_name = amr_name_from_comment(line[:-1])
                continue
            return line, amr_name


def reverse_sign(sign):
    """
    To reverse the sign (+/-)
    Parameetr:
        sign: either + or -
    Return:
        the reverse sign
    """
    if sign == "-":
        return "+"
    elif sign == "+":
        return "-"
    else:
        logging.error("ERROR: ivalid sign!")
        sys.exit(1)


def find_node_orient(node):
    """
    To remove specific characters and return the last character of what remains
    as the orient of the node
    """
    return re.sub("[]}]", "", node)[-1]


def find_node_name(node):
    """
    To remove specific characters and return the rest except the last character
    as the node name
    """
    return re.sub("[]{}[]", "", node)[:-1]


def find_node_name_orient(node):
    """
    To remove specific characters and return the rest as the name+orient of the node
    """
    return re.sub("[]{}[]", "", node)


def exist_in_path(path, mynode):
    """
    To check if a given node exists in the path
    Parameters:
        path:	a list of nodes
        mynde:	the node to check if it is in the path
    Return:
        True if mynode is in the path; False otherwise
    """
    for i, node in enumerate(path):
        if find_node_name_orient(node) == mynode:
            return i
    return -1


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
            logging.error(message)
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
        logging.error(message)
        sys.exit(1)


def run_RGI(
    input_file, output_dir, seq_description, include_loose=False, delete_rgi_files=False
):
    """
    To run RGI and annotate AMRs in the sequence
    # To ensure consistency between Prokka and RGI findings, we annotate found proteins
    # by Prokka (instead of annotationg DNA sequences from scratch)
    Parameters:
        input_file: the file contating proteins annotated by Prokka
        output_dir:  the path for the output directory
        seq_description: a small description of the sequence used for naming
        include_loose: Whether to include loose annotations
    Return:
        the list of extracted annotation information for the sequence
    """
    rgi_dir = os.path.join(output_dir, "rgi_dir")
    if not os.path.exists(rgi_dir):
        try:
            os.makedirs(rgi_dir)
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
            pass

    output_file_name = os.path.join(
        rgi_dir,
        "rgi_output_"
        + seq_description
        + "_"
        + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M"),
    )
    # remove any potential * from the sequence
    delete_a_string_from_file("*", input_file)
    arg_list = [
        "rgi",
        "main",
        "--input_sequence",
        input_file,
        "--output_file",
        output_file_name,
        "--input_type",
        "protein",
        "--clean",
        "--exclude_nudge",
    ]
    if include_loose:
        arg_list.append("--include_loose")
    rgi_command = subprocess.run(arg_list, stdout=subprocess.PIPE, check=True)
    logging.info(rgi_command.stdout.decode("utf-8"))
    seq_info_list = []
    if os.path.isfile(output_file_name + ".txt"):
        with open(output_file_name + ".txt", newline="") as rgi_file:
            rgi_reader = csv.reader(rgi_file, delimiter="\t")
            next(rgi_reader)
            for row in rgi_reader:
                seq_info = {
                    "ORF_ID": row[0],
                    "gene": row[8].strip(),
                    "prediction_type": row[5].strip(),
                    "best_identities": float(row[9]),
                    "family": row[16].strip(),
                }
                seq_info_list.append(seq_info)
    else:
        logging.error("ERROR: RGI didn't run successfully!")
        sys.exit(1)
    # delete temp files
    if delete_rgi_files and os.path.isfile(output_file_name + ".txt"):
        os.remove(output_file_name + ".txt")
    if delete_rgi_files and os.path.isfile(output_file_name + ".json"):
        os.remove(output_file_name + ".json")

    return seq_info_list


def annotate_sequence(
    seq,
    seq_description,
    output_dir,
    no_RGI=False,
    RGI_include_loose=False,
    delete_prokka_dir=False,
):
    """
    To run Prokka for a sequence and extract required information from its
        generated output files
    Parameters:
        seq:	the sequence to be annotated
        seq_description: a small description of the sequence used for naming
        output_dir:  the path for the output directory
        no_RGI:	RGI annotations incorporated for AMR annotation
    Return:
        the list of extracted annotation information for the sequence
    """
    # write the sequence into a temporary file
    seq_file_name = create_fasta_file(
        seq,
        "",
        file_name="temp_"
        + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        + seq_description,
    )
    pid = os.getpid()
    prokka_dir = (
        "prokka_dir_"
        + seq_description
        + "_"
        + str(pid)
        + "_"
        + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
    )
    prefix_name = "mygenome_" + seq_description
    arg_list = [
        "prokka",
        "--metagenome",
        "--outdir",
        prokka_dir,
        "--prefix",
        prefix_name,
        "--fast",
        "--notrna",
        seq_file_name,
    ]
    #cwd = os.getcwd()
    #PROKKA_COMMAND_PREFIX = 'docker run -v '+cwd+':/data staphb/prokka:latest '
    #pre_list = PROKKA_COMMAND_PREFIX.strip().split(" ")
    #arg_list = pre_list + arg_list
    prokka_command = subprocess.run(arg_list, stdout=subprocess.PIPE, check=True)
    logging.info(prokka_command.stdout.decode("utf-8"))
    # move prokka directory to the right address
    shutil.move(prokka_dir, os.path.join(output_dir, prokka_dir))
    prokka_dir = os.path.join(output_dir, prokka_dir)
    RGI_output_list = None
    if not no_RGI:
        RGI_output_list = run_RGI(
            os.path.join(prokka_dir, prefix_name + ".faa"),
            output_dir,
            seq_description,
            RGI_include_loose,
            delete_prokka_dir,
        )

    # Go over Prokka's output files and extract required information
    seq_info = []
    with open(os.path.join(prokka_dir, prefix_name + ".tsv"), "r") as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        # skip the header
        next(reader)
        for row in reader:
            mygene = row[3].strip()
            split_gene = mygene.split("_")
            if len(split_gene) == 2 and split_gene[1].isnumeric():
                mygene = split_gene[0]
            gene_info = {
                "locus_tag": row[0].strip(),
                "gene": mygene,
                "length": row[2].strip(),
                "product": row[6].strip(),
                "start_pos": None,
                "end_pos": None,
                "prokka_gene_name": mygene,
                "RGI_prediction_type": None,
                "coverage": None,
                "family": None,
                "seq_value": seq[:-1],
                "seq_name": None,
                "target_amr": None,
            }
            seq_info.append(gene_info)
    counter = 0
    with open(os.path.join(prokka_dir, prefix_name + ".tbl"), "r") as read_obj:
        for line in read_obj:
            if line[0].isdigit():
                cells = line.split("\t")
                seq_info[counter]["start_pos"] = int(cells[0])
                seq_info[counter]["end_pos"] = int(cells[1])
                counter += 1

    # incorporate RGI findings into Prokka's
    if RGI_output_list:
        for item in RGI_output_list:
            for gene_info in seq_info:
                if item["ORF_ID"].split(" ")[0] == gene_info["locus_tag"]:
                    gene_info["gene"] = item["gene"]
                    gene_info["RGI_prediction_type"] = item["prediction_type"]
                    gene_info["family"] = item["family"]
                    break

    # remove temporary files and folder
    if os.path.isfile(seq_file_name):
        os.remove(seq_file_name)
    if delete_prokka_dir:
        try:
            shutil.rmtree(prokka_dir)
        except OSError as e:
            logging.error("Error: %s - %s." % (e.filename, e.strerror))

    return seq_info


def split_up_down_info(sequence, seq_info):
    """ """
    amr_start = -1
    amr_end = -1
    index = 0
    up_info = []
    down_info = []
    while index < len(sequence):
        if sequence[index].islower() and amr_start == -1:
            amr_start = index
        elif sequence[index].isupper() and amr_start > -1:
            amr_end = index - 1
            break
        index += 1
    # if there is no downstream
    if amr_end == -1 and sequence[-1].islower():
        amr_end = len(sequence) - 1
    elif amr_end == -1 or amr_start == -1:
        logging.error("No AMR sequence (lower case string) was found in " + sequence)
        import pdb

        pdb.set_trace()
    # import pdb;pdb.set_trace()
    # find the gene has the most overlap with the found range
    overlap_thr = 50
    found = False
    amr_info = []
    for gene_info in seq_info:
        start, end = min(gene_info["start_pos"], gene_info["end_pos"]), max(
            gene_info["start_pos"], gene_info["end_pos"]
        )
        if end < amr_start:
            up_info.append(gene_info)
        elif start > amr_end:
            down_info.append(gene_info)
        else:
            # added by 1 because in string indecesstarts from 0
            diff = max((amr_start + 1 - start), 0) + max((end - (amr_end + 1)), 0)
            if ((1 - (float(diff) / (end - start))) * 100) > overlap_thr:
                found = True
                gene_info["target_amr"] = "yes"
                amr_info = gene_info
            elif start < amr_start:
                up_info.append(gene_info)
            else:
                down_info.append(gene_info)

    return found, amr_info, up_info, down_info, seq_info


def compare_two_sequences(
    subject,
    query,
    output_dir,
    threshold=90,
    switch_allowed=True,
    return_file=False,
    subject_coverage=True,
    blast_ext="",
):
    """
    To compare one sequence (shorter sequence) against the other one (longer sequence) using blastn
    """
    # make sure subject is the longer sequence
    if switch_allowed and len(subject) < len(query):
        subject, query = query, subject
    # write the query sequence into a fasta file
    query_file_name = os.path.join(output_dir, "query.fasta")
    with open(query_file_name, "w") as query_file:
        query_file.write("> query \n")
        query_file.write(query)
    # write the query sequence into a fasta file
    subject_file_name = os.path.join(output_dir, "subject.fasta")
    with open(subject_file_name, "w") as subject_file:
        subject_file.write("> subject \n")
        subject_file.write(subject)
    # run blast query for alignement
    blast_file_name = os.path.join(output_dir, "blast" + blast_ext + ".csv")
    blast_file = open(blast_file_name, "w")
    subprocess.run(
        [
            "blastn",
            "-query",
            query_file_name,
            "-subject",
            subject_file_name,
            "-task",
            "blastn-short",
            "-outfmt",
            "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp",
        ],
        stdout=blast_file,
        check=True,
    )
    blast_file.close()

    if return_file:
        return blast_file_name

    with open(blast_file_name, "r") as file1:
        myfile = csv.reader(file1)
        for row in myfile:
            identity = int(float(row[2]))
            coverage = int(float(row[3]) / len(subject) * 100)
            q_coverage = int(float(row[12]))
            if subject_coverage and identity >= threshold and coverage >= threshold:
                return True
            if (
                not subject_coverage
                and identity >= threshold
                and q_coverage >= threshold
            ):
                return True
    return False


def unnamed_genes_are_siginificantly_similar(
    gene_info1, gene_info2, output_dir, threshold=90
):
    """ """
    if gene_info1["gene"] != "" or gene_info2["gene"] != "":
        return False
    start1, end1 = min(gene_info1["start_pos"], gene_info1["end_pos"]), max(
        gene_info1["start_pos"], gene_info1["end_pos"]
    )
    seq1 = gene_info1["seq_value"][start1 - 1 : end1 - 1]
    start2, end2 = min(gene_info2["start_pos"], gene_info2["end_pos"]), max(
        gene_info2["start_pos"], gene_info2["end_pos"]
    )
    seq2 = gene_info2["seq_value"][start2 - 1 : end2 - 1]
    return compare_two_sequences(seq1, seq2, output_dir, threshold)


def seqs_annotation_are_identical(seq_info1, seq_info2, out_dir, threshold=90):
    """ """
    if len(seq_info1) == len(seq_info2):
        identical_rows = 0
        for i, gene_info1 in enumerate(seq_info1):
            gene_info2 = seq_info2[i]
            if (
                gene_info1["gene"] == gene_info2["gene"] and gene_info1["gene"] != ""
            ) or (
                gene_info1["gene"] == gene_info2["gene"]
                and unnamed_genes_are_siginificantly_similar(
                    gene_info1, gene_info2, out_dir, threshold
                )
            ):
                identical_rows += 1
        if identical_rows == len(seq_info1):
            return True
    return False


def similar_seq_annotation_already_exist(
    seq_info_list, all_seq_info_lists, out_dir, threshold=90
):
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
        if seqs_annotation_are_identical(seq_info_list, seq_list, out_dir, threshold):
            found = True
            break

    return found


def extract_info_from_overlap_file(overlap_file_name):
    """ """
    heads = []
    member_lists = []
    unique_amr_list = []
    with open(overlap_file_name, "r") as read_obj:
        for line in read_obj:
            if ":" in line:
                items = line[:-1].split(":")
                if len(items[1]) > 0:
                    heads.append(items[0])
                    members = items[1].split(", ")
                    member_lists.append(members)
                else:
                    unique_amr_list.append(items[0])
    return heads, member_lists, unique_amr_list


def extract_unique_align_files(all_align_files, unique_amr_files):
    """ """
    amr_align_files = []
    if all_align_files:
        for amr_file in unique_amr_files:
            found_it = False
            amr_name = extract_name_from_file_name(amr_file)
            for align_file in all_align_files:
                if os.path.basename(align_file).startswith(amr_name + "_align"):
                    found_it = True
                    amr_align_files.append(align_file)
                    break
            if not found_it:
                logging.error("no alignment was found for " + amr_file)
                import pdb

                pdb.set_trace()
    return amr_align_files


def extract_nodes_in_path(path):
    """
    Parameters:
        path:	a list of nodes with -/+ tail and comma separated : e.g., '(1363) 69-, 2193+ (1786)'
    Return:
        node_list:	list of node numbers -> e.g., [69, 2193]
        orientation_list: list of orientation of nodes -> e.g., [-, +]
        start_pos:	where in the first node the sequence has started -> e.g., 1363
        end_pos:	where in the last node the sequence ended  -> e.g., 1786
    """
    start_pos = 0
    end_pos = 0
    if path.startswith("("):
        index = path.find(")")
        start_pos = int(path[1:index])
    if path.endswith(")"):
        index = path.rfind("(")
        end_pos = int(path[index + 1 : -1])
    # Remove text between ()
    path = (re.sub("\((.*?)\)", "", path)).strip()
    node_list = []
    orientation_list = []
    nodes = path.split(",")
    for node in nodes:
        if "-" in node:
            orientation_list.append("-")
        else:
            orientation_list.append("+")
        node = re.sub("[+-]", "", node.split()[0])
        node_list.append(node)
    return node_list, orientation_list, start_pos, end_pos


def read_path_info_from_align_file(align_file, threshold=95):
    paths_info = []
    found = False
    with open(align_file) as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        # skip the header
        next(reader)
        for row in reader:
            coverage = float(re.sub("[%]", "", row[3]))
            identity = float(re.sub("[%]", "", row[5]))
            if int(coverage) >= threshold and int(identity) >= threshold:
                found = True
                cell_info = row[1].strip()
                nodes, orientation_list, start_pos, end_pos = extract_nodes_in_path(
                    cell_info
                )
                path_info = {
                    "nodes": nodes,
                    "orientations": orientation_list,
                    "start_pos": start_pos,
                    "end_pos": end_pos,
                }
                paths_info.append(path_info)
    if not found:
        logging.info("ERROR: no path info was found in " + align_file)
    return found, paths_info


def read_path_info_from_align_file_with_multiple_amrs(align_file, threshold=99):
    """ """
    paths_info_list = collections.defaultdict(list)
    with open(align_file) as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        # skip the header
        next(reader)
        for row in reader:
            amr_name = restricted_amr_name_from_modified_name(
                row[0].split("|")[-1].strip().replace(" ", "_").replace("'", ";")
            )
            coverage = float(re.sub("[%]", "", row[3]))
            identity = float(re.sub("[%]", "", row[5]))
            if int(coverage) >= threshold and int(identity) >= threshold:
                cell_info = row[1].strip()
                nodes, orientation_list, start_pos, end_pos = extract_nodes_in_path(
                    cell_info
                )
                path_info = {
                    "nodes": nodes,
                    "orientations": orientation_list,
                    "start_pos": start_pos,
                    "end_pos": end_pos,
                }
                paths_info_list[amr_name].append(path_info)
    return paths_info_list


def extract_path_info_for_amrs(all_align_files, unique_amr_files, amr_count, threshold):
    """ """
    unique_amr_path_list = []
    if len(all_align_files) == amr_count:
        amr_align_files = extract_unique_align_files(all_align_files, unique_amr_files)
        for align_file in amr_align_files:
            found, paths_info = read_path_info_from_align_file(align_file, threshold)
            if found:
                unique_amr_path_list.append(paths_info)
            else:
                logging.error(align_file + " file was not found or was empty!")
                import pdb

                pdb.set_trace()
    else:
        paths_info_group_list = []
        for align_file in all_align_files:
            paths_info_group = read_path_info_from_align_file_with_multiple_amrs(
                align_file, threshold
            )
            paths_info_group_list.append(paths_info_group)
        unique_amr_path_list = []
        for amr_file in unique_amr_files:
            amr_found = False
            restricted_amr_name = extract_name_from_file_name(amr_file)
            for paths_info_group in paths_info_group_list:
                if restricted_amr_name in paths_info_group:
                    amr_found = True
                    path_info = paths_info_group[restricted_amr_name]
                    unique_amr_path_list.append(path_info)
            if not amr_found:
                logging.error(
                    "ERROR: no path info was found for " + restricted_amr_name
                )
                import pdb

                pdb.set_trace()
    return unique_amr_path_list


def delete_lines_started_with(ch, filename):
    """
    To delete all the lines in a text file that starts with a given character
    Parameters:
        ch: the character
        filename: the text file
    """
    # command = "sed -i '/^P/d' " + file_name
    # os.system(command)
    file1 = open(filename, "r")
    # file2 = open('temp.txt', 'w')
    file2 = open("temp_" + os.path.basename(filename), "w")
    for line in file1.readlines():
        if not (line.startswith(ch)):
            file2.write(line)
    file1.close()
    file2.close()
    os.rename("temp_" + os.path.basename(filename), filename)


def delete_a_string_from_file(ch, filename):
    """
    To delete a given character or string from a file
    Parameters:
        ch: the character or string to be deleted
        filename: the text file
    """
    with open(filename, "r") as infile, open(
        "temp_" + os.path.basename(filename), "w"
    ) as outfile:
        data = infile.read()
        data = data.replace(ch, "")
        outfile.write(data)
    os.rename("temp_" + os.path.basename(filename), filename)


def check_file(path: str) -> Path:
    """
    Check an input file exists and is readable
    """
    path = Path(path)
    if path.exists() and path.is_file():
        return path.resolve()
    else:
        raise argparse.ArgumentTypeError(
            f"{path} can't be found, please double check the path"
        )


def check_dependencies(programs):
    """
    Check all dependencies exist and work
    """
    missing = False
    for program in programs:
        try:
            program_name = program.split()[0]
            output = subprocess.run(
                program,
                shell=True,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                encoding="utf-8",
            )
            #logging.debug(f"Tool {program_name} is installed: {output.stdout.strip()}")
            print(f"Tool {program_name} is installed: {output.stdout.strip()}")
        except:
            #logging.error(f"Tool {program_name} is not installed")
            print(f"Tool {program_name} is not installed")
            missing = True

    if missing:
        #logging.error("One or more dependencies are missing please install")
        print("One or more dependencies are missing, please install")
        sys.exit(1)
    else:
        #logging.debug("All dependencies found")
        print("All dependencies found")


def validate_range(value_type, minimum, maximum):
    """Determine whether arg value is within minimum and maximum range"""

    def range_checker(arg):
        """
        argparse type function to determine value is within specified range
        """
        if value_type is float:
            try:
                val = float(arg)
            except ValueError:
                raise argparse.ArgumentTypeError("Must be a float")
        elif value_type is int:
            try:
                val = int(arg)
            except ValueError:
                raise argparse.ArgumentTypeError("Must be an int")

        if val < minimum or val > maximum:
            raise argparse.ArgumentTypeError(f"must be in range [{minimum}-{maximum}]")
        return val

    return range_checker
