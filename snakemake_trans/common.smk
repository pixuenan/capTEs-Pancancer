import re
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import pysam
import csv
import pickle
import glob
import numpy as np
import pandas as pd

from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sklearn.cluster import DBSCAN


direction_dict = {"forward": "reverse", "reverse": "forward"}

gff_columns = ["chr", 1, "type", "start", "end", 2, "strand", 3, "info"]

kozak1 = list("GCCGCCATGG")
kozak2 = list("GCCACCATGG")

def get_trans_info(gff_df):
    """Extract exon number and strand information of every transcript."""
    trans_info_dict = {}
    trans_idx = []
    line_num = gff_df.shape[0]
    for idx, row in gff_df.iterrows():
        if not idx == 0 and row["type"] == "transcript" or idx == line_num-1:
            trans_df = gff_df.loc[trans_idx]
            exon_num = trans_df.shape[0]-1
            strand = trans_df["strand"].unique().tolist()[0]
            trans_id = gff_df.loc[idx-1, "trans_id"]
            trans_info_dict.setdefault(trans_id, {}).update({"exon_num": exon_num, "strand": strand})
            trans_idx = [idx]
        else:
            trans_idx.append(idx)
    return trans_info_dict

def check_exon_order(exon_id, exon_num, strand):
    """Check exon order in the transcript."""
    if strand == ".":
        if exon_id in [1, exon_num]:
            exon_order = "first"
        else:
            exon_order = "middle"
    elif strand == "+":
        if exon_id == 1:
            exon_order = "first"
        elif exon_id == exon_num:
            exon_order = "last"
        else:
            exon_order = "middle"
    elif strand == "-":
        if exon_id == exon_num:
            exon_order = "first"
        elif exon_id == 1:
            exon_order = "last"
        else:
            exon_order = "middle" 
    return exon_order
    
def write_trans2fa_500(select_gff, ref_records_dict, output_fa_prefix):
    """write gff file to fa, every fasta file contains 500 records"""
    result_dict = {}
    current_trans = None
    seq_ls = []
    for idx, row in select_gff.iterrows():
        if not current_trans:
            current_trans = row["trans_id"]
            continue
        if row["trans_id"] == current_trans:
            seq_ls.append(str(ref_records_dict[row["chr"]].seq[row["start"]-1:row["end"]]))
            continue
        else:
            result_dict[current_trans] = "".join(seq_ls)
            current_trans = row["trans_id"]
            seq_ls = []
            continue
        if idx == select_gff.shape[0]-1:
            result_dict[current_trans] = "".join(seq_ls)
    file_num = int(len(result_dict)/500)+1
    trans_ls = list(result_dict.keys())
    for i in range(file_num):
        output_fa = output_fa_prefix+".{}.fa".format(i)
        with open(output_fa, "w+") as handle:
            for key in trans_ls[i*500:min((i+1)*500,len(trans_ls))]:
                value = result_dict[key]
                SeqIO.write(SeqRecord(Seq(value), id=key, description=""), handle, "fasta")

def check_reverse_ele(ele_ls, dir_ls, level, info_dict):
    reve_pair_ls = []
    group_pair_ls = []
    for idx, ele in enumerate(ele_ls):
        if idx == len(ele)-1:
            continue
        te_type = info_dict[ele][level]
        te_dir = dir_ls[idx]
        for n_idx, n_ele in enumerate(ele_ls[idx+1:]):
            n_te_type = info_dict[n_ele][level]
            n_te_dir = dir_ls[idx+1+n_idx]
            if n_te_dir != te_dir and n_te_type == te_type:
                group_pair_ls.append((te_type, info_dict[ele]['repFamily'], info_dict[n_ele]['repFamily'],
                                      info_dict[ele]['repName'], info_dict[n_ele]['repName']))
                reve_pair_ls.append((ele, n_ele))
    return reve_pair_ls, group_pair_ls

def build_read_exon(row):
    """construct exon coordinate from flair all_corrected.bed """
    exon_len_ls = list(map(int, row["exon_len"].split(",")[:-1]))
    start_diff_ls = list(map(int, row["start_diff"].split(",")[:-1]))
    exon_ls = []
    for idx, value in enumerate(exon_len_ls):
        if value <= 0:
            continue
        exon_ls.append([row["start"]+start_diff_ls[idx]+1, row["start"]+start_diff_ls[idx]+value])
    start, end = exon_ls[0][0], exon_ls[-1][-1]
    splice_ls = []
    for idx, value in enumerate(exon_ls):
        if idx == len(exon_ls)-1:
            continue
        n_value = exon_ls[idx+1]
        splice_ls.append([value[-1], n_value[0]])
    read_exon_ls = [[start]]
    for value in splice_ls:
        read_exon_ls[-1].extend([value[0], row["read_id"]])
        read_exon_ls.append([value[1]])
    read_exon_ls[-1].extend([end, row["read_id"]])
    read_df = pd.DataFrame(read_exon_ls)
    read_df["chr"] = row["chr"]
    read_df["strand"] = row["strand"]
    return read_df

def construct_exon(row):
    """construct exon coordinate from flair all_corrected.bed """
    exon_len_ls = list(map(int, row["exon_len"].split(",")[:-1]))
    start_diff_ls = list(map(int, row["start_diff"].split(",")[:-1]))
    exon_ls = []
    for idx, value in enumerate(exon_len_ls):
        if value < 0:
            continue
        exon_ls.append([row["chr"], row["start"]+start_diff_ls[idx], 
                        row["start"]+start_diff_ls[idx]+value, row["read_id"]])  
    return exon_ls

def score_kozak_seq(in_str):
    kozak1_sum = sum([i==kozak1[idx] for idx, i in enumerate(list(in_str))])
    kozak2_sum = sum([i==kozak2[idx] for idx, i in enumerate(list(in_str))])
    score = max([kozak1_sum, kozak2_sum])
    return score

def is_valid_protein_seq(sequence):
    """Checks if a protein sequence is valid.
    Args:
        sequence (str): The protein sequence.
    Returns:
        bool: True if the sequence is valid, False otherwise.
    """
    valid_amino_acids = set("ARNDCEQGHILKMFPSTWYV")
    return all(aa in valid_amino_acids for aa in sequence)

def get_z_score_read(read_df, trans_p_dict, read_id, trans_id):
    indicator = read_df.loc[read_id, trans_id]
    denominator = sum([read_df.loc[read_id, t_id]*trans_p_dict[t_id] for t_id in trans_p_dict.keys()])
    return trans_p_dict[trans_id]*indicator/denominator

def get_trans_p(trans_read_z_dict, read_number, trans_id):
    return sum(trans_read_z_dict[trans_id].values())/read_number

def clean_overlap_trans(df, read_dict, trans_dict):
    """select between transcripts supported by the same set of reads"""
    remove_ls = []
    trans_ls = df["trans_id"].unique().tolist()
    trans_ele_dict = {i: get_trans_ele(read_dict, trans_dict[i]['path'][:-1]) for i in trans_ls}
    skip_idx_ls = [len(trans_ls)-1]
    for idx, trans_id in enumerate(trans_ls):
        if idx in skip_idx_ls:
            continue
        trans_ele = trans_ele_dict[trans_id]
        cur_trans_df = df[df["trans_id"]==trans_id]
        for n_idx, n_trans_id in enumerate(trans_ls[idx+1:]):
            n_trans_df = df[df["trans_id"]==n_trans_id]
            diff_valid = compare_trans_dis(cur_trans_df, n_trans_df)
            n_trans_ele = trans_ele_dict[n_trans_id]
            if len(trans_ele.difference(n_trans_ele)) or not diff_valid:
                continue
            diff_part = select_trans(trans_dict[trans_id]["path"], trans_dict[n_trans_id]["path"], trans_id, n_trans_id, trans_dict)
            if diff_part == "double":
                remove_ls.extend([trans_id, n_trans_id])
                skip_idx_ls.append(n_idx+idx+1)
            elif diff_part == 1:
                remove_ls.append(trans_id)
                break
            elif diff_part == 2:
                remove_ls.append(n_trans_id)
                skip_idx_ls.append(n_idx+idx+1)
    return df[~df["trans_id"].isin(remove_ls)]

def get_trans_ele(read_dict, trans_part_ls):
    ele_ls = []
    for part in trans_part_ls:
        if part in read_dict.keys():
            ele_ls.extend(read_dict[part])
        else:
            ele_ls.append(part)
    return set(ele_ls)

def get_te_subfamily(type_df, te_dict, region_ls, sample_id):
    """summary te type information and intergenic information for the transcripts"""
    te_sub_ls = []
    te_class_ls = []
    te_family_ls = []
    te_id_ls = []
    intergenic_flag_ls = []
    region_id = None
    for idx, row in type_df.iterrows():
        if "-".join(idx.split("_")[:3])!=region_id:
            region_id = "-".join(idx.split("_")[:3])
            with open("07.remapping/{}/obj/{}.trans.path.dict.obj".format(sample_id, region_id), "rb") as handle:
                path_dict = pickle.load(handle)
            flag = False
            for region in region_ls:
                if region in idx:
                    flag = True
                    break
        intergenic_flag_ls.append(flag)
        ele_ls = path_dict[idx]["path"]
        sub_ls, class_ls, family_ls, te_ls = [], [], [], []
        for ele in ele_ls:
            ele_prefix = "_".join(ele.split("_")[:-1])
            if ele_prefix in ["Alu", "long_TE"]:
                te_subtype = te_dict[ele]["repName"]
                te_class = te_dict[ele]["repClass"]
                te_family = te_dict[ele]["repFamily"]
                sub_ls.append(te_subtype)
                class_ls.append(te_class)
                family_ls.append(te_family)
                te_ls.append(ele)
        te_sub_ls.append("|".join(sub_ls))
        te_class_ls.append("|".join(set(class_ls)))
        te_family_ls.append("|".join(set(family_ls)))
        te_id_ls.append("|".join(te_ls))
    type_df["te_subfamily"] = te_sub_ls
    type_df["te_class"] = te_class_ls
    type_df["te_family"] = te_family_ls
    type_df["intergenic"] = intergenic_flag_ls
    type_df["te_id"] = te_id_ls
    return type_df

def write_trans2fa(select_gff, ref_records_dict, output_fa):
    """write gff file to fa"""
    result_dict = {}
    current_trans = None
    seq_ls = []
    for idx, row in select_gff.iterrows():
        if not current_trans:
            current_trans = row["trans_id"]
            continue
        if row["trans_id"] == current_trans:
            seq_ls.append(str(ref_records_dict[row["chr"]].seq[row["start"]-1:row["end"]]))
            continue
        else:
            result_dict[current_trans] = "".join(seq_ls)
            current_trans = row["trans_id"]
            seq_ls = []
            continue
        if idx == select_gff.shape[0]-1:
            result_dict[current_trans] = "".join(seq_ls)
    with open(output_fa, "w+") as handle:
        for key, value in result_dict.items():
            SeqIO.write(SeqRecord(Seq(value), id=key, description=""), handle, "fasta") 


def construct_bed(region_ls):
    region_df = pd.DataFrame([i.split("_") for i in region_ls])
    region_df[1] = region_df[1].astype(int)
    region_df.sort_values(by=[0,1], inplace=True)
    return region_df

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

def count_junction_length(in_s):
    pattern = re.compile(r"\d+[A-Z]")
    result = pattern.findall(in_s)

    length, times = 0, 0
    for value in result:
        if "N" in value:
            length += int(value[:-1])
            times += 1
    return length, times

def reform_info_column(end_trans_df):
    re_end_trans_df = end_trans_df.reset_index()
    re_end_trans_df["info"] = re_end_trans_df.apply(lambda x: x["type"]=="exon" and \
                                                        'gene_id "{}"; transcript_id "{}"; exon_number "{}";'.format(x.gene_id, x.trans_id, int(x.name)-1) \
                                                        or 'gene_id "{}"; transcript_id "{}";'.format(x.gene_id, x.trans_id), axis=1)
    return re_end_trans_df[gff_columns]

def order_trans_df(df):
    trans_order = df[df["type"]=="transcript"].sort_values(by=["start"])["trans_id"]
    new_trans_df_ls = []
    for trans_id in trans_order:
        new_trans_df_ls.append(df[df["trans_id"]==trans_id])
    ordered_gff = pd.concat(new_trans_df_ls)
    return ordered_gff

def order_gff_all(gff_df):
    trans_order = gff_df[gff_df["type"]=="transcript"].sort_values(by=["chr", "start"])["trans_id"]
    new_trans_df_ls = []
    for trans_id in trans_order:
        new_trans_df_ls.append(gff_df[gff_df["trans_id"]==trans_id])
    ordered_gff = pd.concat(new_trans_df_ls)
    return ordered_gff

def get_support_read(chr_name, region_start, region_end, bam_content, select_te_df,
                     select_te_dict, trans_ls, path_dict, trans_dict):
    """build matrix of reads transcripts relationship"""
    global for_ls, rev_ls
    for_ls, rev_ls = [], []
    for read in bam_content.fetch(chr_name, region_start, region_end):
        if read.reference_start < region_start or read.reference_end > region_end:
            continue
        if read.is_secondary:
            continue
        left_clip, right_clip = get_read_softclip(read)
        subject1 = "".join(["T" for i in range(max([20, left_clip-20]))])
        subject2 = "".join(["A" for i in range(max([20, right_clip-20]))])
        if left_clip > 20 and left_clip < 400 and read.is_reverse and \
            seq2_aligner.align("AGCATACGA", read.seq[:min([left_clip, 30])]).score >= 7 and \
            seq2_aligner.align(subject1, read.seq[left_clip-15:left_clip+15]).score <25 and \
            seq2_aligner.align(subject1, read.seq[:left_clip]).score>=max([15, left_clip-25]):
            rev_ls.append(read.qname)
        if right_clip > 20 and right_clip < 400 and not read.is_reverse and \
            seq2_aligner.align("TCGTATGCT", read.seq[-min([right_clip, 30]):]).score >= 7 and \
            seq2_aligner.align(subject2, read.seq[-right_clip-15:-right_clip+15]).score <25 and \
            seq2_aligner.align(subject2, read.seq[-right_clip:]).score>=max([15, right_clip-25]):
            for_ls.append(read.qname)

    read_dict = {}
    te_dict = {}
    ## get reads besides cut site
    for idx, row in select_te_df.iterrows():
        te_ob = TE(row["te_id"])
        cor_ls = select_te_dict[row["te_id"]]
        te_ob.get_cut_coordinate(cor_ls, row["start"], row["end"])
        read_content = bam_content.fetch(chr_name, row["start"], row["end"])
        read_dict = te_ob.get_cut_read_ls(read_content, read_dict)
        te_dict[row["te_id"]] = te_ob
    ## record polyA reads at TE
    for idx, row in select_te_df.iterrows():
        te_obj = te_dict[row["te_id"]]
        te_obj.get_end_read_ls(read_dict)
        te_obj.clean_polyA_in_end_reads(for_ls, rev_ls)
        te_dict[te_obj.id] = te_obj

    trans_read_dict = {}
    read_dict = {"unspliced": [], "spliced": []}
    for trans_id in trans_ls:
        trans_ele = path_dict[trans_id]["path"]
        trans_ele_valid = []
        direc = trans_ele[-1]
        ele_num = len(trans_ele)-1
        for idx, i in enumerate(trans_ele[:-1]):
            if not "_".join(i.split("_")[:-1]) in ["Alu", "long_TE"]:
                trans_ele_valid.append(i)
            ## record unspliced end reads supporting transcript
            else:
                te_obj = te_dict[i]
                if idx==0:
                    read_ls = direc=="forward" and te_obj.reve_end_ls or te_obj.forw_end_ls
                elif idx==ele_num-1:
                    read_ls = direc=="forward" and te_obj.forw_poly_read or te_obj.reve_poly_read
                ## for TE sit between two elements, the end reads were not count
                else:
                    continue
                read_dict["unspliced"].extend(read_ls)
                for j in read_ls:
                    trans_read_dict.setdefault(trans_id, {}).update({j:1})
        ## record spliced reads supporting transcript
        for i in trans_ele_valid:
            read_dict["spliced"].extend(trans_dict[i])
            for j in trans_dict[i]:
                trans_read_dict.setdefault(trans_id, {}).update({j:1})

    read_df = pd.DataFrame(trans_read_dict).fillna(0)
    return read_df, read_dict


def get_candidate_path(chr_name, region_start, region_end, bam_content, select_te_df,
                       select_te_dict, clean_select_gff, keep_trans_dict, sample_name):
    ## extract read with poly A at the soft clipping part with poly length exceeds 20 bp plus seq2 primer
    global for_ls, rev_ls
    for_ls, rev_ls = [], []
    for read in bam_content.fetch(chr_name, region_start, region_end):
        if read.reference_start < region_start or read.reference_end > region_end:
            continue
        if read.is_secondary:
            continue
        left_clip, right_clip = get_read_softclip(read)
        subject1 = "".join(["T" for i in range(max([20, left_clip-20]))])
        subject2 = "".join(["A" for i in range(max([20, right_clip-20]))])
        if left_clip > 20 and left_clip < 400 and read.is_reverse and \
            seq2_aligner.align("AGCATACGA", read.seq[:min([left_clip, 30])]).score >= 7 and \
            seq2_aligner.align(subject1, read.seq[left_clip-15:left_clip+15]).score <25 and \
            seq2_aligner.align(subject1, read.seq[:left_clip]).score>=max([15, left_clip-25]):
            rev_ls.append(read.qname)
        if right_clip > 20 and right_clip < 400 and not read.is_reverse and \
            seq2_aligner.align("TCGTATGCT", read.seq[-min([right_clip, 30]):]).score >= 7 and \
            seq2_aligner.align(subject2, read.seq[-right_clip-15:-right_clip+15]).score <25 and \
            seq2_aligner.align(subject2, read.seq[-right_clip:]).score>=max([15, right_clip-25]):
            for_ls.append(read.qname)

    read_dict = {}
    te_dict = {}
    remove_idx_ls = []
    ## get reads besides cut site
    for idx, row in select_te_df.iterrows():
        te_ob = TE(row["te_id"])
        cor_ls = select_te_dict[row["te_id"]]
        te_ob.get_cut_coordinate(cor_ls, row["start"], row["end"])
        ## remove TEs that are not target of sgRNA
        if not te_ob.target:
            remove_idx_ls.append(idx)
            continue

        read_content = bam_content.fetch(chr_name, row["start"], row["end"])
        read_dict = te_ob.get_cut_read_ls(read_content, read_dict)
        te_dict[row["te_id"]] = te_ob
    select_te_df.drop(index=remove_idx_ls, inplace=True)
    select_te_df.reset_index(inplace=True)
    ## get transcript link between TE
    for idx, row in select_te_df.iterrows():
        te_obj = te_dict[row["te_id"]]
        te_obj.get_end_read_ls(read_dict)
        te_obj.clean_polyA_in_end_reads(for_ls, rev_ls)
        te_obj.get_spliced_read_ls(read_dict)
        ## check if te directly link to transcript
        te_obj.get_link_transcripts(keep_trans_dict, row)
        te_dict[te_obj.id] = te_obj

    ## link te by END reads which span the cut site
    ## clean END reads
    for idx, row in select_te_df.iterrows():
        if idx == select_te_df.shape[0]-1:
            continue
        local_te_obj = te_dict[row["te_id"]]
        if local_te_obj.clean_forw_end_ls:
            local_te_obj.get_link_te(select_te_df[idx+1:], te_dict, read_dict)
        local_te_obj.get_end_read_coordinate(read_dict)
        te_dict[local_te_obj.id] = local_te_obj

    for idx, row in select_te_df.iterrows():
        if idx == select_te_df.shape[0]-1:
            local_te_obj = te_dict[row["te_id"]]
            local_te_obj.get_end_read_coordinate(read_dict)
            te_dict[local_te_obj.id] = local_te_obj

    ## add reverse direction te link
    for idx, row in select_te_df.iterrows():
        if idx == 0:
            continue
        local_te_obj = te_dict[row["te_id"]]
        for t_idx, t_row in select_te_df[:idx].iterrows():
            pre_te_obj = te_dict[t_row["te_id"]]
            if not local_te_obj.id in pre_te_obj.forw_link_te:
                continue
            local_te_obj.reve_link_te.append(t_row["te_id"])
        te_dict[local_te_obj.id] = local_te_obj

    ## get te link by transcripts
    for idx, row in select_te_df.iterrows():
        local_te_obj = te_dict[row["te_id"]]
        if not local_te_obj.trans_dict:
            continue
        local_te_obj.get_link_te_by_trans(select_te_df, te_dict)
        local_te_obj.get_link_end_trans()
        te_dict[local_te_obj.id] = local_te_obj

    ## select TE with link transcript end with polyA
    trans_polya_ls = []
    candidate_trans_polya_path = []
    for te_id, te_obj in te_dict.items():
        if not te_obj.trans_dict:
            continue
        for trans_id, trans_info in te_obj.trans_dict.items():
            forw_polya_reads, reve_polya_reads = [], []
            ## count reads mapping to transcript with polyA
            if trans_info[0] == "forward":
                forw_polya_reads = set(trans_info[1]).intersection(for_ls)
            else:
                reve_polya_reads = set(trans_info[1]).intersection(rev_ls)
            forw_polya_count, reve_polya_count = len(forw_polya_reads), len(reve_polya_reads)
            if forw_polya_count >= 3 and reve_polya_count >= 3:
                continue
            elif trans_id in te_obj.forw_link_trans:
                if (forw_polya_count >= 2 and forw_polya_count/len(trans_info[1]) >= 0.4) or forw_polya_count >= 5:
                    candidate_trans_polya_path.append([te_id, trans_id, "forward"])
                    trans_polya_ls.append(trans_id)
            elif trans_id in te_obj.reve_link_trans:
                if (reve_polya_count >= 2 and reve_polya_count/len(trans_info[1]) >= 0.4) or reve_polya_count >= 5: 
                    candidate_trans_polya_path.append([te_id, trans_id, "reverse"])
                    trans_polya_ls.append(trans_id)
            elif trans_id in te_obj.forw_link_trans:
                if (forw_polya_count >= 2 and forw_polya_count/len(trans_info[1]) >= 0.4) or forw_polya_count >= 5: 
                    candidate_trans_polya_path.append([te_id, trans_id, "forward"])
                    trans_polya_ls.append(trans_id)

    trans_polya_path = []
    for te_id, trans_id, direction in candidate_trans_polya_path:
        current_path = [te_id, trans_id]
        for i in find_path_end(current_path, direction, te_dict):
            trans_polya_path.append(i+[direction])

    clean_trans_polya_path = []
    ## remove path with elements with polyA besides transcript at the end
    for t_path in trans_polya_path:
        for ele in t_path[:-2]:
            if ele in trans_polya_ls:
                continue
        ele_prefix = [i.split("_")[0] for i in t_path[:-1]]
        if ele_prefix.count("Alu") >= 7:
            continue
        clean_trans_polya_path.append(t_path)

    te_polya_path_ls = []
    for te_id, te_obj in te_dict.items():
        if te_obj.forw_poly_count >= 3:
            current_path = [te_id]
            for i in find_path_end(current_path, "forward", te_dict):
                te_polya_path_ls.append(i+["forward"])
        if te_obj.reve_poly_count >= 3:
            current_path = [te_id]
            for i in find_path_end(current_path, "reverse", te_dict):
                te_polya_path_ls.append(i+["reverse"])

    clean_te_polya_path = []
    for t_path in te_polya_path_ls:
        ## remove path only contain TEs
        flag = True
        ele_prefix = [i.split("_")[0] for i in t_path[:-1]]
        if set(ele_prefix) == {"Alu"}:
            continue
        ## exclude the path with more than 7 Alus
        if ele_prefix.count("Alu") >= 7:
            continue
        for i in t_path[:-1]:
            if not "_".join(i.split("_")[:-1]) in ["Alu", "long_TE"]:
                flag = False
        if flag:
            continue
        ## remove path with intermediate transcripts with polyA
        for ele in t_path[:-2]:
            if ele in trans_polya_ls:
                continue
        clean_te_polya_path.append(t_path)

    trans_path_dict = {}
    ## the element placed at the end of the list was the one with polyA 
    for p_idx, cur_path in enumerate(clean_te_polya_path + clean_trans_polya_path):
        direction = cur_path[-1]
        ## path end with transcript
        if not "_".join(cur_path[-2].split("_")[:-1]) in ["Alu", "long_TE"]:
            end_trans = cur_path[-2]
            end_trans_df = clean_select_gff[clean_select_gff["trans_id"]==end_trans]
            to_added_eles = list(reversed(cur_path[:-2]))
            end_trans_df = trans_extend(to_added_eles, end_trans_df, te_dict, clean_select_gff, direction, trans_polya=True)
        ## path end with TE
        else:
            ## construct transcript firstly by walk to the nearest transcript at the end
            for idx, ele in enumerate(list(reversed(cur_path[:-1]))):
                if not "_".join(ele.split("_")[:-1]) in ["Alu", "long_TE"]:
                    end_trans = ele
                    break
            end_trans_df = clean_select_gff[clean_select_gff["trans_id"]==end_trans]
            to_added_eles = list(cur_path[-idx-1:-1])
            end_trans_df = trans_extend(to_added_eles, end_trans_df, te_dict, clean_select_gff, direction_dict[direction])
            to_added_eles = list(reversed(cur_path[:-idx-2]))
            end_trans_df = trans_extend(to_added_eles, end_trans_df, te_dict, clean_select_gff, direction, trans_polya=True)
        new_id = "{}_{}_{}_{}_{}".format(chr_name, region_start, region_end, p_idx, sample_name)
        new_gene_id = "{}_{}_{}".format(chr_name, region_start, region_end)
        end_trans_df["trans_id"] = new_id
        end_trans_df["gene_id"] = new_gene_id
        if "level_0" in end_trans_df.columns:
            end_trans_df.drop(columns=["index", "level_0"], inplace=True)
        else:
            end_trans_df.drop(columns=["index"], inplace=True)
        end_trans_df["strand"] = direction == "forward" and "+" or "-"
        new_end_trans_df = reform_info_column(end_trans_df)
        trans_path_dict.setdefault(new_id, {}).update({"path": cur_path})
        trans_path_dict.setdefault(new_id, {}).update({"gff": new_end_trans_df})
    return te_dict, trans_path_dict

def get_valid_transcript(trans_fa, bam_content, gff_df):
    """collect valid transcript information from remapping result"""
    trans_records = SeqIO.parse(trans_fa, "fasta")
    trans_records_dict = SeqIO.to_dict(trans_records)
    trans_len_dict = {i:len(trans_records_dict[i].seq) for i in trans_records_dict.keys()}

    keep_trans_dict = {}
    for trans_id in trans_len_dict.keys():
        trans_len = trans_len_dict[trans_id]
        valid_read = get_valid_read(trans_id, trans_len, bam_content)
        if len(valid_read) >= 3:
            exon_num = gff_df[(gff_df["trans_id"]==trans_id) & (gff_df["type"]=="exon")].shape[0]
            start = gff_df.loc[(gff_df["trans_id"]==trans_id) & (gff_df["type"]=="transcript"), "start"].tolist()[0]
            end = gff_df.loc[(gff_df["trans_id"]==trans_id) & (gff_df["type"]=="transcript"), "end"].tolist()[0]
            keep_trans_dict[trans_id] = {"read_ls": valid_read, "exon_num": exon_num, "start": start, "end": end}
    return keep_trans_dict

def get_valid_read(trans_id, trans_len, bam_content):
    """Get reads that remapped to the transcript and pass the filter"""
    keep_ls = []
    for read in bam_content.fetch(trans_id):
        ## filter out read with distance to the end of the transcript greater than 1 kb
        if read.reference_start > 100 and trans_len-read.reference_end > 100:
            continue
        ## filter out read 150bp longer than transcript or longer than 10% of the lenght of the transcript
        len_dif = read.infer_read_length() - trans_len
        if len_dif > 150 or len_dif > trans_len*0.1:
            continue
        ## filter out read with left or right soft clipping longer than 100 bp
        if read.qstart > 100 or read.query_length-read.qend > 100:
            continue
        ## filter out read that cannot span 80% of the transcript
        if read.qlen/trans_len < 0.8:
            continue
        flag = False
        for pair in read.cigartuples:
            ## filter out reads with deletion longer than 50bp
            if pair[0] == 2 and pair[1] > 50:
                flag = True
                break
            ## filter out reads with insertion longer than 50bp
            if pair[0] == 1 and pair[1] > 50:
                flag = True
                break
        if flag:
            continue
        keep_ls.append(read.qname)
    return keep_ls

def clean_transcripts(trans_ls, path_gff_ls, dist_thre=30):
    remove_trans_ls = []
    [i.reset_index(inplace=True) for i in path_gff_ls]
    for idx, cur_trans_df in enumerate(path_gff_ls):
        if 0 in cur_trans_df["start"].tolist() or 0 in cur_trans_df["end"].tolist():
            remove_trans_ls.append(idx)
        if idx == len(trans_ls)-1:
            continue
        ## remove single exon transcript
        if cur_trans_df.shape[0] == 2:
            remove_trans_ls.append(idx)
            continue
        for n_idx, n_trans_df in enumerate(path_gff_ls[idx+1:]):
            ## remove transcript with identical exon number
            if cur_trans_df.shape[0] == n_trans_df.shape[0] and n_trans_df.shape[0] > 2:
                trans_distance = abs(cur_trans_df["start"]-n_trans_df["start"]).sum()+abs(cur_trans_df["end"]-n_trans_df["end"]).sum()
                if trans_distance < dist_thre:
                    remove_trans_ls.append(n_idx+idx+1)
    return remove_trans_ls

def clean_transcripts_te(trans_ls, path_gff_ls, path_dict, dist_thre=30):
    remove_trans_ls = []
    [i.reset_index(inplace=True) for i in path_gff_ls]
    ## record te element in each trans
    te_ls = []
    for trans_id in trans_ls:
        te_ele = []
        path_ele = path_dict[trans_id]["path"]
        for i in path_ele:
            if "_".join(i.split("_")[:-1]) in ["Alu", "long_TE"]:
                te_ele.append(i)
        te_ls.append(set(te_ele))
    for idx, cur_trans_df in enumerate(path_gff_ls):
        if 0 in cur_trans_df["start"].tolist() or 0 in cur_trans_df["end"].tolist():
            remove_trans_ls.append(idx)
        if idx == len(trans_ls)-1:
            continue
        ## remove single exon transcript
        if cur_trans_df.shape[0] == 2:
            remove_trans_ls.append(idx)
            continue
        for n_idx, n_trans_df in enumerate(path_gff_ls[idx+1:]):
            if te_ls[idx] != te_ls[idx+n_idx+1]:
                continue
            ## remove transcript with identical exon number
            if cur_trans_df.shape[0] == n_trans_df.shape[0] and n_trans_df.shape[0] > 2:
                trans_distance = abs(cur_trans_df["start"]-n_trans_df["start"]).sum()+abs(cur_trans_df["end"]-n_trans_df["end"]).sum()
                if trans_distance < dist_thre:
                    remove_trans_ls.append(n_idx+idx+1)
    return remove_trans_ls

def get_valid_transcript_info(bam_file, fa_file, gff_df):
    if os.path.getsize(bam_file) <= 500:
        keep_trans_dict = {}
        select_gff = pd.DataFrame()
    else:
        bam_content = pysam.AlignmentFile(bam_file, "rb")
        keep_trans_dict = get_valid_transcript(fa_file, bam_content, gff_df)
        select_gff = gff_df.loc[gff_df["trans_id"].isin(keep_trans_dict.keys())]
    return keep_trans_dict, select_gff

