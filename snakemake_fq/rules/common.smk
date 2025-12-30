import re
import pysam
import csv
import pickle
import glob
import numpy as np
import pandas as pd

from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

gff_columns = ["chr", 1, "type", "start", "end", 2, "strand", 3, "info"]

def get_fq_path_m(wildcards):
    sample_df = pd.read_csv("/NAS/wg_pxn/Project/TEcap/sample_sheet/SampleInfo.txt", sep="\t")
    fq_folder = sample_df.loc[sample_df["SampleName"]==wildcards.sample_name, "fastq_dir"].item().split(",")
    return fq_folder[0]

def get_fq_in_string(wildcards):
    sample_df = pd.read_csv("/NAS/wg_pxn/Project/TEcap/sample_sheet/SampleInfo.txt", sep="\t")
    fq_folder = sample_df.loc[sample_df["SampleName"]==wildcards.sample_name, "fastq_dir"].item().split(",")
    s_value = " ".join([i+"/*.fastq.gz" for i in fq_folder])
    return s_value

def count_exon_number(gff_df):
    """count exon number in every transcript"""
    num_dict = {}
    current_trans = None
    exon_count = 0
    for idx, row in gff_df.iterrows():
        if not current_trans:
            current_trans = row["trans_id"]
            continue
        if row["trans_id"] == current_trans:
            exon_count += 1
        else:
            num_dict[current_trans] = exon_count
            current_trans = row["trans_id"]
            exon_count = 0
        if idx == gff_df.shape[0]-1:
            num_dict[current_trans] = exon_count
    return num_dict
