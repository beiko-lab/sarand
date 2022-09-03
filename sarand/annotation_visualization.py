"""
File:		annotation_visualization.py
Author:		Somayeh Kafaie
Date:		September 2020
Purpose:	To visualize sequences annotations and make their comparison easier

To run:
	python annotation_visualization.py --csvfile <annotation file path>
    --title <the image title> --output <the output image name>

"""

from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
from PIL import Image
import argparse
from csv import DictReader
import logging
import shutil
import os


def show_images(image_list, main_title, output, cols=1, title_list=None):
    """
    Display a list of images in a single figure with matplotlib
    Parameters:
            image_list: List of np.arrays compatible with plt.imshow
            main_title: The main title of the generated figure
            output: The address of the generated image to be saved in
            cols (Default = 1): Number of columns in figure (number of rows is
                    set to np.ceil(n_images/float(cols)))
            title_list: List of titles corresponding to each image. Must have
            the same length as titles
    """
    assert (title_list is None) or (len(image_list) == len(title_list))
    n_images = len(image_list)
    if title_list is None:
        title_list = ["Image (%d)" % i for i in range(1, n_images + 1)]
    fig = plt.figure()
    for n, (image, title) in enumerate(zip(image_list, title_list)):
        _ = fig.add_subplot(math.ceil(n_images / float(cols)), cols, n + 1)
        plt.imshow(image)
        ax = plt.gca()
        ax.axes.xaxis.set_ticks([])
        ax.axes.yaxis.set_ticks([])
        plt.xlabel(title, size=12)
    fig.suptitle(main_title, size=18)
    fig.set_size_inches(np.array(fig.get_size_inches()) * n_images)
    plt.savefig(output)


def extract_annotation_from_csv(input_csv_file):
    """
    To extract the annotations of sequences from the csv file
    Parameters:
            input_csv_file: the input csv file containing annotations
    Return:
            The list of annotaions of all sequences, their length and their name (title)
    """
    seq_info_list = []
    seq_length_list = []
    seq_info = []
    title_list = []
    with open(input_csv_file, "r") as myfile:
        myreader = DictReader(myfile)
        old_seq = ""
        has_row = False
        for row in myreader:
            has_row = True
            if row["coverage"] != "":
                coverage = " (" + str(round(float(row["coverage"]), 2)) + ")"
            else:
                coverage = ""
            gene_info = {
                "name": row["gene"],
                "start_pos": int(row["start_pos"]),
                "end_pos": int(row["end_pos"]),
                "coverage": coverage,
            }
            cur_seq = row["seq_name"]
            if cur_seq != old_seq:
                if seq_info:
                    seq_info_list.append(seq_info)
                seq_info = []
                old_seq = cur_seq
                title_list.append(cur_seq)
                seq_length_list.append(int(row["seq_length"]))
            seq_info.append(gene_info)
        seq_info_list.append(seq_info)
    if not has_row:
        logging.error("there is no sequence to visualize!")
        return [], [], []
    if len(seq_info_list) > 20:
        logging.error("there are more than 20 sequences to visualize!")
        return [], [], []
    return seq_info_list, seq_length_list, title_list


def visualize_annotation(input_csv_file, output, title=""):
    """
    The core function to read the csv file containing annotations into image_list
    object and visualize them in an image
    Parameetrs:
            input_csv_file:	the csv file containing the annotation for all extracted sequences
            output:			the address of the generated image to be saved in
            title:			the title used in the generated image
    """
    seq_info_list, seq_length_list, title_list = extract_annotation_from_csv(
        input_csv_file
    )
    if not seq_info_list:
        return
    # create a temporary folder to store tmporary images
    temp_dir = "tmp_image/"
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    image_list = []
    for counter, seq_info in enumerate(seq_info_list):
        if not seq_info:
            break
        features = []
        for gene_info in seq_info:
            features.append(
                GraphicFeature(
                    start=gene_info["start_pos"],
                    end=gene_info["end_pos"],
                    strand=+1,
                    color="#ffd700",
                    label=gene_info["name"] + gene_info["coverage"],
                )
            )
        record = GraphicRecord(
            sequence_length=seq_length_list[counter], features=features
        )
        ax, _ = record.plot(figure_width=10)
        ax.figure.savefig(
            os.path.join(temp_dir, "temp" + str(counter) + ".jpg"), bbox_inches="tight"
        )
        img = Image.open(os.path.join(temp_dir, "temp" + str(counter) + ".jpg"))
        image_list.append(img)
    show_images(
        image_list=image_list, main_title=title, output=output, title_list=title_list
    )
    # delete temporary directory
    if os.path.exists(temp_dir):
        try:
            shutil.rmtree(temp_dir)
        except OSError as e:
            logging.error("Error: %s - %s." % (e.filename, e.strerror))


def main(args):
    """ """
    if args.csvfile:
        visualize_annotation(args.csvfile, args.output, args.title)
    else:
        logging.error(
            "please enter the path for the csv file containing the sequences annotation"
        )
        sys.exit()


if __name__ == "__main__":

    text = "This code is used to visualize the sequences annotation extracted from a csv file"
    parser = argparse.ArgumentParser(description=text)
    parser.add_argument(
        "--csvfile",
        "-C",
        type=str,
        default="/media/Data/PostDoc/Dalhousie/Work/Test1/Experiments/Set30_13Sep2020/\
    salmonella_10-1/annotation_1000/annotation_detail.csv",
        help="the path of the csv file containing the sequences annotation",
    )
    parser.add_argument(
        "--title", "-T", type=str, default="", help="the title for the generated image"
    )
    parser.add_argument(
        "--output",
        "-O",
        type=str,
        default="gene_comparison.png",
        help="the address of the output file",
    )
    args = parser.parse_args()
    main(args)
