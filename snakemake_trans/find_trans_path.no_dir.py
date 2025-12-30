import re
import os
import sys
import argparse
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

def merge_exon_intervals(gff_df):
    """
    Merges overlapping intervals in a pandas DataFrame with 'start' and 'end' columns.
    """
    df = gff_df[gff_df["type"]=="exon"]
    df = df.sort_values(by='start')  # Sort by start for efficient merging
    merged_intervals = []
    current_start = df['start'].iloc[0]
    current_end = df['end'].iloc[0]
    for i in range(1, len(df)):
        start = df['start'].iloc[i]
        end = df['end'].iloc[i]

        # Overlap condition: current interval overlaps with the next interval
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged_intervals.append([current_start, current_end])
            current_start = start
            current_end = end
    merged_intervals.append([current_start, current_end])  # Append the last interval
    return merged_intervals

def get_z_score_read(read_df, trans_p_dict, read_id, trans_id):
    indicator = read_df.loc[read_id, trans_id]
    denominator = sum([read_df.loc[read_id, t_id]*trans_p_dict[t_id] for t_id in trans_p_dict.keys()])
    return trans_p_dict[trans_id]*indicator/denominator

def cal_trans_read_dict(trans_ls, path_dict, te_dict, trans_dict):
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

def get_read_mapq(bam_c):
    map_q = {}
    for read in bam_content.fetch():
        map_q[read.qname] = read.mapq
    return map_q

def order_trans_df(df):
    trans_order = df[df["type"]=="transcript"].sort_values(by=["start"])["trans_id"]
    new_trans_df_ls = []
    for trans_id in trans_order:
        new_trans_df_ls.append(df[df["trans_id"]==trans_id])
    ordered_gff = pd.concat(new_trans_df_ls)
    return ordered_gff

def select_trans(trans_1_part, trans_2_part, trans_1_id, trans_2_id, path_dict):
    """return the transcript number to remove"""
    trans_1_uniq = list(set(trans_1_part).difference(set(trans_2_part)))
    trans_2_uniq = list(set(trans_2_part).difference(set(trans_1_part)))
    if len(trans_1_uniq) == 1 and len(trans_2_uniq) == 1:
        if set(trans_1_uniq[0]+trans_2_uniq[0]) == {'forward', 'reverse'}:
            return "double"
        else:
            ele_dict = {trans_1_uniq[0]: 2, trans_2_uniq[0]: 1}
            ## select the ENST trans
            for ele in [trans_1_uniq[0], trans_2_uniq[0]]:
                if ele.startswith("ENST"):
                    return ele_dict[ele]
            ## both stringtie, select the one with high exon number and short length
            if trans_1_uniq[0].startswith("STRG") and trans_2_uniq[0].startswith("STRG"):
                uniq_1_gff = path_dict[trans_1_id]["gff"]
                uniq_2_gff = path_dict[trans_2_id]["gff"]
                trans_1_exon_num = uniq_1_gff.shape[0]-1
                trans_2_exon_num = uniq_2_gff.shape[0]-1
                if trans_1_exon_num == trans_2_exon_num:
                    trans_1_len = uniq_1_gff.loc[0,"end"] - uniq_1_gff.loc[0,"start"]
                    trans_2_len = uniq_2_gff.loc[0,"end"] - uniq_2_gff.loc[0,"start"]
                    return trans_1_len >= trans_2_len and 1 or 2
                return trans_1_exon_num <= trans_2_exon_num and 1 or 2
            ## select the stringtie trans
            for ele in [trans_1_uniq[0], trans_2_uniq[0]]:
                if ele.startswith("STRG"):
                    return ele_dict[ele]
    else:
        return None

def compare_trans_dis(a_trans_df, b_trans_df):
    """Compare two transcripts, return True
    if they differed at the fisrt or last exon,
    otherwise False"""
    cur_trans_df = a_trans_df.reset_index()
    n_trans_df = b_trans_df.reset_index()
    exon_num = cur_trans_df.shape[0]-1
    start_diff = cur_trans_df["start"] - n_trans_df["start"]
    end_diff = cur_trans_df["end"] - n_trans_df["end"]
    mid_dis = start_diff[2:exon_num-1].sum()+end_diff[2:exon_num-1].sum()
    end_dis = start_diff[1]+start_diff[exon_num]+end_diff[1]+end_diff[exon_num]
    if mid_dis==0:
        if end_dis < 30:
            return True
        else:
            return False
    else:
        return False

def build_read_splice(row):
    """construct exon coordinate from flair all_corrected.bed """
    exon_len_ls = list(map(int, row["exon_len"].split(",")[:-1]))
    start_diff_ls = list(map(int, row["start_diff"].split(",")[:-1]))
    exon_ls = []
    for idx, value in enumerate(exon_len_ls):
        exon_ls.append([row["start"]+start_diff_ls[idx]+1, row["start"]+start_diff_ls[idx]+value])
    start, end = exon_ls[0][0], exon_ls[-1][-1]
    splice_ls = []
    for idx, value in enumerate(exon_ls):
        if idx == len(exon_ls)-1:
            continue
        n_value = exon_ls[idx+1]
        splice_ls.append([value[-1], n_value[0]])
    return {"start": start, "end": end, "splice_ls": splice_ls}

def reform_info_column(end_trans_df):
    re_end_trans_df = end_trans_df.reset_index()
    re_end_trans_df["info"] = re_end_trans_df.apply(lambda x: x["type"]=="exon" and \
                                                    'gene_id "{}"; transcript_id "{}"; exon_number "{}";'.format(x.gene_id, x.trans_id, int(x.name)-1) \
                                                    or 'gene_id "{}"; transcript_id "{}";'.format(x.gene_id, x.trans_id), axis=1)
    return re_end_trans_df[gff_columns]

def count_junction_length(in_s):
    pattern = re.compile(r"\d+[A-Z]")
    result = pattern.findall(in_s)

    length, times = 0, 0
    for value in result:
        if "N" in value:
            length += int(value[:-1])
            times += 1
    return length, times

def trans_extend(to_added_eles, end_trans_df, te_dict, clean_select_gff, direction, trans_polya=False):
    ele_num = len(to_added_eles)
    for idx, ele in enumerate(to_added_eles):
        end_trans_index = end_trans_df.index
        ## skip the TE that in the middle of two transcripts
        if not "_".join(ele.split("_")[:-1]) in ["Alu", "long_TE"] \
           and idx!=ele_num-1 and not "_".join(to_added_eles[idx+1][:2].split("_")[:-1]) in ["Alu", "long_TE"]:
                continue
        elif "_".join(ele.split("_")[:-1]) in ["Alu", "long_TE"]:
            te_id = ele
            ## foward direction, extend the first exon
            if direction == "forward":
                te_start = te_dict[te_id].reve_end_start
                end_trans_df.loc[end_trans_index[0], "start"] = te_start
                end_trans_df.loc[end_trans_index[1], "start"] = te_start
            ## reverse direction, extend the last exon
            elif direction == "reverse":
                te_end = te_dict[te_id].forw_end_end
                end_trans_df.loc[end_trans_index[-1], "end"] = te_end
                end_trans_df.loc[end_trans_index[0], "end"] = te_end
        else:
            cur_trans_gff = clean_select_gff[clean_select_gff["trans_id"]==ele]
            cur_trans_gff_index = cur_trans_gff.index
            ## foward direction, extend the first exon
            if direction == "forward":
                cur_trans_last_exon_start = cur_trans_gff["start"].tolist()[-1]
                cur_trans_start = cur_trans_gff["start"].tolist()[0]
                end_trans_df.loc[end_trans_index[1], "start"] = cur_trans_last_exon_start
                end_trans_df.loc[end_trans_index[0], "start"] = cur_trans_start
                end_trans_df = pd.concat([pd.DataFrame(end_trans_df.loc[end_trans_index[0]]).T,
                                          cur_trans_gff.loc[cur_trans_gff_index[1:-1]],
                                          end_trans_df.loc[end_trans_index[1:]]])
            ## reverse direction, extend the last exon
            elif direction == "reverse":
                cur_trans_first_exon_end = cur_trans_gff["end"].tolist()[1]
                cur_trans_end = cur_trans_gff["end"].tolist()[0]
                end_trans_df.loc[end_trans_index[-1], "end"] = cur_trans_first_exon_end
                end_trans_df.loc[end_trans_index[0], "end"] = cur_trans_end
                end_trans_df = pd.concat([end_trans_df, cur_trans_gff.loc[cur_trans_gff_index[2:]]])
    return end_trans_df.reset_index()

def get_read_softclip(read):
    left_soft, right_soft = 0, 0
    if read.cigartuples[0][0] == 4:
        left_soft = read.cigartuples[0][1]
    if read.cigartuples[-1][0] == 4:
        right_soft = read.cigartuples[-1][1]
    return left_soft, right_soft

def find_path_end(current_path, direction, te_dict):
    ## end the path with more than 7 Alus
    ele_prefix = [i.split("_")[0] for i in current_path]
    if ele_prefix.count("Alu") >= 7:
        return current_path
    te_id = current_path[0]
    current_te_obj = te_dict[te_id]
    if direction == "reverse":
        if current_te_obj.clean_forw_end_ls:
            yield current_path
        if current_te_obj.forw_link_end_trans:
            for trans_id in current_te_obj.forw_link_end_trans:
                te_path = current_path[:]
                te_path.insert(0, trans_id)
                yield te_path
        if current_te_obj.forw_link_te:
            for te_id in current_te_obj.forw_link_te:
                te_path = current_path[:]
                te_path.insert(0, te_id)
                yield from find_path_end(te_path, direction, te_dict)
        if current_te_obj.forw_trans_link_te:
            for te_id, trans_id in current_te_obj.forw_trans_link_te:
                te_path = current_path[:]
                te_path.insert(0, trans_id)
                te_path.insert(0, te_id)
                yield from find_path_end(te_path, direction, te_dict)
    elif direction == "forward":
        if current_te_obj.clean_reve_end_ls:
            yield current_path
        if current_te_obj.reve_link_end_trans:
            for trans_id in current_te_obj.reve_link_end_trans:
                te_path = current_path[:]
                te_path.insert(0, trans_id)
                yield te_path
        if current_te_obj.reve_link_te:
            for te_id in current_te_obj.reve_link_te:
                te_path = current_path[:]
                te_path.insert(0, te_id)
                yield from find_path_end(te_path, direction, te_dict)
        if current_te_obj.reve_trans_link_te:
            for te_id, trans_id in current_te_obj.reve_trans_link_te:
                te_path = current_path[:]
                te_path.insert(0, trans_id)
                te_path.insert(0, te_id)
                yield from find_path_end(te_path, direction, te_dict)

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

def clean_gff(path_dict):
    """remove redundant transcript by DBSCAN clustering"""
    trans_num_collect_dict = {}
    for trans_id, info_dict in path_dict.items():
        exon_num = info_dict["gff"].shape[0]-1
        ## remove transcript with start or end 0
        gff_df = info_dict["gff"]
        trans_df = gff_df[gff_df['type']=="transcript"]
        t_idx = trans_df.index.tolist()[0]
        if trans_df.loc[t_idx, "start"] == 0 or trans_df.loc[t_idx, "end"] == 0:
            continue
        trans_num_collect_dict.setdefault(exon_num, []).append(trans_id)

    keep_trans_ls = []
    for exon_num in trans_num_collect_dict.keys():
        trans_cor_metrix = []
        trans_ls = trans_num_collect_dict[exon_num]
        for trans_id in trans_ls:
            cor_df = path_dict[trans_id]["gff"][["start", "end"]]
            cor_df.drop(0, inplace=True)
            cor_ls = []
            for idx, row in cor_df.iterrows():
                cor_ls.extend(row.tolist())
            trans_cor_metrix.append(cor_ls)
        db = DBSCAN(eps=0.2, min_samples=2).fit(trans_cor_metrix)
        labels = db.labels_
        for label in set(labels):
            if label == -1:
                for idx, value in enumerate(trans_ls):
                    if labels[idx] == -1:
                        keep_trans_ls.append(trans_ls[idx])
                continue
            identical_trans_id_ls = []
            for idx, value in enumerate(trans_ls):
                if labels[idx] == label:
                    identical_trans_id_ls.append(trans_ls[idx])
            keep_trans_ls.append(identical_trans_id_ls[0])

    keep_trans_df = [path_dict[i]["gff"] for i in keep_trans_ls]
    clean_path_gff = pd.concat(keep_trans_df).reset_index()
    clean_path_gff.drop(columns=["index"], inplace=True)
    return clean_path_gff

class READ:
    def __init__(self, read_id):
        self.id = read_id

    def get_read_info(self, read):
        self.start = read.reference_start
        self.end = read.reference_end
        self.spliced = count_junction_length(read.cigarstring)[0] > 50 and True or False

class TE:
    def __init__(self, te_id):
        self.id = te_id
        self.trans_dict = {}
        self.target = True
        self.te_link_trans = []
        self.forw_link_te, self.reve_link_te = [], []
        self.forw_trans_link_te, self.reve_trans_link_te = [], []
        self.forw_link_trans, self.reve_link_trans = [], []
        self.forw_link_end_trans, self.reve_link_end_trans = [], []
        self.forw_poly_count, self.reve_poly_count = 0, 0

    def get_cut_coordinate(self, cor_ls, te_start, te_end):
        self.start, self.end = te_start, te_end
        flag = False
        result_cor_ls = []
        for cut_cor in cor_ls:
            cut_start, cut_end = cut_cor[0], cut_cor[1]
            if cut_start+cut_end==0:
                flag = True
            elif self.id.startswith("Alu"):
                result_cor_ls.append((te_start, te_end))
            else:
                result_cor_ls.append((min([cut_start-50, te_start]), max([cut_end+50, te_end])))
        if flag:
            self.target = False
        else:
            self.cor_ls = result_cor_ls

    def get_cut_read_ls(self, read_content, read_dict, supp_reads):
        self.forw_read_ls = []
        self.reve_read_ls = []
        for read in read_content:
            if read.mapq < 10:
                continue
            if read.is_secondary or read.qname in supp_reads:
                continue
            ## search for reads starts or ends within sgRNA match region
            for te_cor in self.cor_ls:
                te_start, te_end = te_cor[0], te_cor[1]
                if read.reference_start >= te_start and read.reference_start <= te_end and not read.is_reverse:
                    self.forw_read_ls.append(read.qname)
                    read_obj = READ(read.qname)
                    read_obj.get_read_info(read)
                    read_dict.setdefault(read.qname, read_obj)
                if read.reference_end <= te_end and read.reference_end >= te_start and read.is_reverse:
                    self.reve_read_ls.append(read.qname)
                    read_obj = READ(read.qname)
                    read_obj.get_read_info(read)
                    read_dict.setdefault(read.qname, read_obj)
        return read_dict

    def get_end_read_ls(self, read_dict):
        ## check if te had END reads
        self.forw_end_ls = [i for i in self.forw_read_ls if not read_dict[i].spliced]
        self.reve_end_ls = [i for i in self.reve_read_ls if not read_dict[i].spliced]
        self.clean_forw_end_ls = self.forw_end_ls[:]
        self.clean_reve_end_ls = self.reve_end_ls[:]

    def get_end_read_coordinate(self, read_dict):
        forw_end_ls = [read_dict[i].end for i in self.clean_forw_end_ls]
        reve_start_ls = [read_dict[i].start for i in self.clean_reve_end_ls]
        self.forw_end_end = forw_end_ls and np.bincount(forw_end_ls).argmax() or self.end
        self.reve_end_start = reve_start_ls and np.bincount(reve_start_ls).argmax() or self.start

    def get_spliced_read_ls(self, read_dict):
        self.forw_sp_read_ls = [i for i in self.forw_read_ls if read_dict[i].spliced]
        self.reve_sp_read_ls = [i for i in self.reve_read_ls if read_dict[i].spliced]

    def get_link_transcripts(self, trans_dict, row):
        for trans_id, trans_info in trans_dict.items():
            if trans_info["exon_num"] < 2:
                continue
            forw_trans_ls = set(self.forw_sp_read_ls).intersection(set(trans_info["read_ls"]))
            reve_trans_ls = set(self.reve_sp_read_ls).intersection(set(trans_info["read_ls"]))
            if len(forw_trans_ls) >= 3 and trans_info["start"] > row["start"]:
                self.trans_dict[trans_id] = ("forward", forw_trans_ls)
                self.forw_link_trans.append(trans_id)
            if len(reve_trans_ls) >= 3 and trans_info["end"] < row["end"]:
                self.trans_dict[trans_id] = ("reverse", reve_trans_ls)
                self.reve_link_trans.append(trans_id)

    def get_cut_range(self, read_dict):
        local_forw_start_ls = [read_dict[i].start for i in self.clean_forw_end_ls]
        local_reve_end_ls = [read_dict[i].end for i in self.clean_reve_end_ls]
        start_bincount = np.bincount(local_forw_start_ls)
        end_bincount = np.bincount(local_reve_end_ls)
        if start_bincount.any() and end_bincount.any():
            cut_site = start_bincount.max() >= end_bincount.max() and start_bincount.argmax() or end_bincount.argmax()
            self.cut_range = (cut_site-20, cut_site+20)
        elif start_bincount.any():
            cut_site = start_bincount.argmax()
            self.cut_range = (cut_site-20, cut_site+20)
        elif end_bincount.any():
            cut_site = end_bincount.argmax()
            self.cut_range = (cut_site-20, cut_site+20)
        else:
            self.cut_range = (0,1)

    def get_link_te(self, te_df, te_dict, read_dict):
        ## link the te in forward direction by END reads which span the cut site
        ## clean END reads
        for t_idx, t_row in te_df.iterrows():
            target_te_obj = te_dict[t_row["te_id"]]
            if not target_te_obj.reve_end_ls:
                continue
            local_cut_ls = set(range(self.cut_range[0], self.cut_range[1]+1))
            target_cut_ls = set(range(target_te_obj.cut_range[0], target_te_obj.cut_range[1]+1))
            target_rev_link_reads = [i for i in target_te_obj.clean_reve_end_ls if read_dict[i].start in local_cut_ls]
            local_for_link_reads = [i for i in self.clean_forw_end_ls if read_dict[i].end in target_cut_ls]
            self.clean_forw_end_ls = list(set(self.clean_forw_end_ls).difference(set(local_for_link_reads)))
            target_te_obj.clean_reve_end_ls = list(set(target_te_obj.clean_reve_end_ls).difference(set(target_rev_link_reads)))
            te_dict[t_row["te_id"]] = target_te_obj
            if len(local_for_link_reads) >= 3 or len(target_rev_link_reads) >= 3:
                self.forw_link_te.append(t_row["te_id"])

    def get_link_te_by_trans(self, te_df, te_dict):
        ## link te by transcript
        self.te_link_trans = []
        idx = te_df[te_df["te_id"]==self.id].index.item()
        for trans_id, trans_info in self.trans_dict.items():
            if trans_info[0] == "reverse":
                for t_idx, t_row in te_df[:idx].iterrows():
                    target_te_obj = te_dict[t_row["te_id"]]
                    if not target_te_obj.trans_dict:
                        continue
                    if trans_id in target_te_obj.trans_dict.keys() and target_te_obj.trans_dict[trans_id][0]=="forward":
                        self.reve_trans_link_te.append((target_te_obj.id, trans_id))
                        self.te_link_trans.append(trans_id)
            elif trans_info[0] == "forward" and idx < te_df.shape[0]-1:
                for t_idx, t_row in te_df[idx+1:].iterrows():
                    target_te_obj = te_dict[t_row["te_id"]]
                    if not target_te_obj.trans_dict:
                        continue
                    if trans_id in target_te_obj.trans_dict.keys() and target_te_obj.trans_dict[trans_id][0]=="reverse":
                        self.forw_trans_link_te.append((target_te_obj.id, trans_id))
                        self.te_link_trans.append(trans_id)

    def get_link_end_trans(self):
        """exclude transcript that link TE in the linked transcript with end"""
        self.forw_link_end_trans = set(self.forw_link_trans).difference(self.te_link_trans)
        self.reve_link_end_trans = set(self.reve_link_trans).difference(self.te_link_trans)

def get_candidate_path(chr_name, region_start, region_end, bam_content, select_te_df,
                       select_te_dict, clean_select_gff, keep_trans_dict, sample_name, hc_trans_id_max, supp_reads):

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
        read_dict = te_ob.get_cut_read_ls(read_content, read_dict, supp_reads)
        te_dict[row["te_id"]] = te_ob
    select_te_df.drop(index=remove_idx_ls, inplace=True)
    select_te_df.reset_index(inplace=True)
    ## get transcript link between TE
    for idx, row in select_te_df.iterrows():
        te_obj = te_dict[row["te_id"]]
        te_obj.get_end_read_ls(read_dict)
        te_obj.get_cut_range(read_dict)
        te_obj.get_spliced_read_ls(read_dict)
        ## check if te directly link to transcript
        te_obj.get_link_transcripts(keep_trans_dict, row)
        te_dict[te_obj.id] = te_obj

    ## link te by END reads which span the cut site
    ## clean END reads
    for idx, row in select_te_df.iterrows():
        local_te_obj = te_dict[row["te_id"]]
        if idx != select_te_df.shape[0]-1 and local_te_obj.clean_forw_end_ls:
            local_te_obj.get_link_te(select_te_df[idx+1:], te_dict, read_dict)
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
    ## select TE with end reads
    te_path_ls = []
    for te_id, te_obj in te_dict.items():
        if not te_obj.clean_forw_end_ls and not te_obj.clean_reve_end_ls:
            continue
        if len(te_obj.clean_forw_end_ls) >= 3 or (te_obj.reve_link_te or te_obj.reve_link_trans and not te_obj.forw_end_ls):
            current_path = [te_id]
            for i in find_path_end(current_path, "forward", te_dict):
                te_path_ls.append(i+["forward"])
        if len(te_obj.clean_reve_end_ls) >= 3 or (te_obj.forw_link_te or te_obj.forw_link_trans and not te_obj.reve_end_ls):
            current_path = [te_id]
            for i in find_path_end(current_path, "reverse", te_dict):
                te_path_ls.append(i+["reverse"])

    trans_path_dict = {}
    cal_path_ls = []
    ## the element placed at the end of the list was the one with polyA
    for p_idx, cur_path in enumerate(te_path_ls):
        ele_s = "|".join(cur_path[:-1])
        if ele_s in cal_path_ls:
            continue
        cal_path_ls.append(ele_s)
        direction = cur_path[-1]
        ## path end with transcript
        if not "_".join(cur_path[-2].split("_")[:-1]) in ["Alu", "long_TE"]:
            end_trans = cur_path[-2]
            end_trans_df = clean_select_gff[clean_select_gff["trans_id"]==end_trans]
            to_added_eles = list(reversed(cur_path[:-2]))
            end_trans_df = trans_extend(to_added_eles, end_trans_df, te_dict, clean_select_gff, direction, trans_polya=True)
        ## path end with TE
        else:
            end_trans = None
            ## construct transcript firstly by walk to the nearest transcript at the end
            for idx, ele in enumerate(list(reversed(cur_path[:-1]))):
                if not "_".join(ele.split("_")[:-1]) in ["Alu", "long_TE"]:
                    end_trans = ele
                    break
            if end_trans:
                end_trans_df = clean_select_gff[clean_select_gff["trans_id"]==end_trans]
                to_added_eles = list(cur_path[-idx-1:-1])
                end_trans_df = trans_extend(to_added_eles, end_trans_df, te_dict, clean_select_gff, direction_dict[direction])
                to_added_eles = list(reversed(cur_path[:-idx-2]))
                end_trans_df = trans_extend(to_added_eles, end_trans_df, te_dict, clean_select_gff, direction)
            else: 
                ## all elements are TE
                fir_ele, last_ele = cur_path[0], cur_path[-2]
                fir_s, fir_e = te_dict[fir_ele].reve_end_start, te_dict[fir_ele].forw_end_end
                last_s, last_e = te_dict[last_ele].reve_end_start, te_dict[last_ele].forw_end_end
                trans_start = min([fir_s, last_s])
                trans_end = max([fir_e, last_e])
                ## skip TE without start or end reads
                if trans_start == 0 or trans_end == 0:
                    continue
                trans_gff_dict = [[0, chr_name, "TE", "transcript", trans_start, trans_end, ".", ".", "."],
                                  [1, chr_name, "TE", "exon", trans_start, trans_end, ".", ".", "."]]
                end_trans_df = pd.DataFrame(trans_gff_dict, columns=["index", "chr", 1, "type", "start", "end", 2, 3, "info"])
        ## remove transcript with start or end 0
        t_idx = end_trans_df.index.tolist()[0]
        if end_trans_df.loc[t_idx, "start"] == 0 or end_trans_df.loc[t_idx, "end"] == 0:
            continue
        new_id = "{}_{}_{}_{}_{}".format(chr_name, region_start, region_end, p_idx+1+hc_trans_id_max, sample_name)
        new_gene_id = "{}_{}_{}".format(chr_name, region_start, region_end)
        end_trans_df["trans_id"] = new_id
        end_trans_df["gene_id"] = new_gene_id
        if "level_0" in end_trans_df.columns:
            end_trans_df.drop(columns=["index", "level_0"], inplace=True)
        else:
            end_trans_df.drop(columns=["index"], inplace=True)
        end_trans_df["strand"] = "."
        new_end_trans_df = reform_info_column(end_trans_df)
        trans_path_dict.setdefault(new_id, {}).update({"path": cur_path})
        trans_path_dict.setdefault(new_id, {}).update({"gff": new_end_trans_df})
    return te_dict, trans_path_dict

def get_trans_df(gff_df, seq_n, start, end):
    trans_df = gff_df[gff_df["type"]=="transcript"]
    select_trans_ls = trans_df.loc[(trans_df["chr"]==seq_n) & (trans_df["start"] > start) \
    & (trans_df["end"] < end), "trans_id"]
    return gff_df[gff_df["trans_id"].isin(select_trans_ls)]

def build_trans_splice(df):
    """construct exon coordinate from transcript gff"""
    start, end = 0, 0
    splice_ls = []
    for idx, row in df.iterrows():
        if idx==0:
            start = row["start"]
        if idx==df.shape[0]-1:
            end = row["end"]
            continue
        splice_ls.append([row["end"], df.loc[idx+1, "start"]])
    return {"start": start, "end": end, "splice_ls": splice_ls}

def find_hc_region(gff_df):
    trans_ls = gff_df['trans_id'].unique().tolist()
    trans_df = gff_df[(gff_df['trans_id'].isin(trans_ls)) & (gff_df["type"]=="transcript")]
    trans_df.set_index("trans_id", inplace=True)
    hc_region = []
    for trans_id in trans_ls:
        t_start, t_end = trans_df.loc[trans_id, "start"], trans_df.loc[trans_id, "end"]
        if not hc_region:
            hc_region.append([t_start, t_end])
        else:
            last_hc = hc_region[-1]
            ## merge region
            if t_start < last_hc[0]+2000:
                new_end = max([last_hc[-1], t_end])
                hc_region.pop()
                hc_region.append([last_hc[0], new_end])
            else:
                hc_region.append([t_start, t_end])
    return hc_region

def check_hc_inter(trans_start, trans_end, hc_region_ls):
    """check transcript intersection with high confidence region"""
    flag = False
    for hc_region in hc_region_ls:
        hc_start, hc_end = hc_region
        hc_inter = set(range(trans_start, trans_end)).intersection(set(range(hc_start, hc_end)))
        ## remove transcript that more than half of which located within high confidence region
        if len(hc_inter) >= 0.5 * (trans_end - trans_start):
            flag = True
    return flag

def select_te_hc(select_te_df, hc_region_ls):
    """select te outside the high confidence region"""
    keep_te_ls = []
    for idx, row in select_te_df.iterrows():
        flag = False
        for hc_region in hc_region_ls:
            hc_start, hc_end = hc_region
            hc_inter = set(range(row["start"], row["end"])).intersection(set(range(hc_start, hc_end)))
            if hc_inter:
               flag = True
        if flag:
            continue
        keep_te_ls.append(row["te_id"])
    return keep_te_ls

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Find link between TE and transcript.')
    parser.add_argument('-rf', '--remap_fd')
    parser.add_argument('-gf', '--gff_fd')
    parser.add_argument('-s', '--sample_name')
    parser.add_argument('-r', '--region_id')
    parser.add_argument('-b', '--te_bed')
    parser.add_argument('-d', '--te_dict')
    args = parser.parse_args()
    stringtie_gff_obj = "{}/{}/{}.R.gff.obj".format(args.gff_fd, args.sample_name, args.sample_name)
    flair_gff_obj = "{}/{}/{}.flair.gff.obj".format(args.gff_fd, args.sample_name, args.sample_name)
    inconsis_gff_obj = "{}/{}/{}.inconsist.gff.obj".format(args.gff_fd, args.sample_name, args.sample_name)
    bam = "{}/{}/bam/{}.spliced.bam".format(args.remap_fd, args.sample_name, args.region_id)
    read_bed = "{}/{}/bed/{}.read.bed".format(args.remap_fd, args.sample_name, args.region_id)
    trans_out_obj = "{}/{}/obj/{}.trans.no_dir.dict.obj".format(args.remap_fd, args.sample_name, args.region_id)
    hc_obj = "{}/{}/obj/{}.trans.dict.obj".format(args.remap_fd, args.sample_name, args.region_id)
    info_ls = args.region_id.split("-")
    chr_name, region_start, region_end = info_ls[0], int(info_ls[1]), int(info_ls[2])

    with open(stringtie_gff_obj, "rb") as handle:
        stringtie_gff = pickle.load(handle)
    with open(flair_gff_obj, "rb") as handle:
        flair_gff = pickle.load(handle)
    with open(inconsis_gff_obj, "rb") as handle:
        inconsis_gff = pickle.load(handle)
    ## read high confidence result and extract high confidence region
    with open(hc_obj, "rb") as handle:
        hc_dict = pickle.load(handle)
    hc_trans_ls = hc_dict and hc_dict['gff']['trans_id'].unique().tolist() or []
    hc_trans_id_max = hc_trans_ls and max([int(i.split("_")[3]) for i in hc_trans_ls]) or 0
    hc_region_ls = hc_dict and find_hc_region(hc_dict['gff']) or []
    hc_exon_region_ls = hc_dict and merge_exon_intervals(hc_dict['gff']) or []

    flair_select_gff = get_trans_df(flair_gff, chr_name, region_start, region_end)
    stringtie_select_gff = get_trans_df(stringtie_gff, chr_name, region_start, region_end)
    inconsis_select_gff = get_trans_df(inconsis_gff, chr_name, region_start, region_end)
    select_gff = pd.concat([stringtie_select_gff, flair_select_gff, inconsis_select_gff])
    ## build splice information of read outside high confidence region
    df = pd.read_csv(read_bed, sep="\t", names=[1,2,3,4,"start","end", "read_id", "mapq", "strand", "s", "e", "info", "exon_num",
                                                            "exon_len", "start_diff"])
    filtered_df = df[df["exon_num"]>1]
    read_dict = {}
    read_q_dict = {}
    for idx, row in filtered_df.iterrows():
        read_dict.setdefault(row['exon_num'], {}).update({row["read_id"]: build_read_splice(row)})
        read_q_dict[row['read_id']] = row['mapq']

    ## clean transcripts by read splice information
    if select_gff.empty:
        trans_ls = []
    else:
        trans_ls = select_gff["trans_id"].unique().tolist()
    trans_info_df = select_gff[(select_gff["trans_id"].isin(trans_ls)) & (select_gff["type"]=="transcript")]
    trans_info_df.set_index("trans_id", inplace=True)
    trans_read_dict = {}
    for trans_id in trans_ls:
        trans_df = select_gff.loc[(select_gff["trans_id"]==trans_id) & (select_gff["type"]=="exon"), ["start", "end"]].reset_index()
        exon_num = trans_df.shape[0]
        if not exon_num in read_dict.keys():
            continue
        cur_trans_info_df = trans_info_df.loc[trans_id]
        if isinstance(cur_trans_info_df, pd.Series):
            trans_start, trans_end = cur_trans_info_df["start"], cur_trans_info_df["end"]
        else:
            trans_start, trans_end = cur_trans_info_df["start"].tolist()[0], cur_trans_info_df["end"].tolist()[0]
        if check_hc_inter(trans_start, trans_end, hc_region_ls):
            continue
        trans_splice = build_trans_splice(trans_df)
        support_ls = []
        for key, value in read_dict[exon_num].items():
            if read_q_dict[key] < 10:
                continue
            ## check if every splice site on read is consistent with transcript
            if value["splice_ls"]==trans_splice["splice_ls"]:
                if value["start"] >= trans_splice["start"]-50 and value["end"] <= trans_splice["end"]+50 and \
                   value["start"] <= trans_splice["start"]+50 and value["end"] >= trans_splice["end"]-50:
                    support_ls.append(key)
        trans_read_dict[trans_id] = support_ls

    keep_trans_dict = {}
    for trans_id, value in trans_read_dict.items():
        if len(value) < 2:
            continue
        trans_df = select_gff.loc[select_gff["trans_id"]==trans_id]
        exon_num = trans_df[trans_df["type"]=="exon"].shape[0]
        start = trans_df.loc[trans_df["type"]=="transcript", "start"].tolist()[0]
        end = trans_df.loc[trans_df["type"]=="transcript", "end"].tolist()[0]
        keep_trans_dict[trans_id] = {"read_ls": value, "exon_num": exon_num, "start": start, "end": end} 
    if select_gff.empty:
        with open(trans_out_obj, "wb") as handle:
            pickle.dump(None, handle)
    else:
        clean_select_gff = select_gff[select_gff["trans_id"].isin(keep_trans_dict.keys())]
     
        # link transcript and TE end reads
        bam_content = pysam.AlignmentFile(bam, "rb")
        with open(args.te_dict, "rb") as handle:
            te_cor_dict = pickle.load(handle)
        te_df = pd.read_csv(args.te_bed, sep="\t", names=["chr", "start", "end", "strand", "repName","repClass","repFamily", "te_id"])
        
        select_te_df = te_df[(te_df["chr"]==chr_name) & (te_df["start"]>=region_start) & (te_df["end"]<=region_end)]
        keep_te_ls = select_te_hc(select_te_df, hc_exon_region_ls)
        select_te_dict = {i:te_cor_dict[i] for i in keep_te_ls}
        # record supplementary reads
        supp_reads = []
        read_content = bam_content.fetch()
        for read in read_content:
            if read.is_supplementary:
                supp_reads.append(read.qname)
        te_dict, trans_path_dict = get_candidate_path(chr_name, region_start, region_end, bam_content, select_te_df[select_te_df["te_id"].isin(keep_te_ls)], 
                                                      select_te_dict, clean_select_gff, keep_trans_dict, args.sample_name, hc_trans_id_max, supp_reads) 
        if not trans_path_dict:
            with open(trans_out_obj, "wb") as handle:
                pickle.dump(None, handle)
        else:
            # remove redudant transcripts
            df = clean_gff(trans_path_dict)
            df["trans_id"] = df["info"].apply(lambda x: x.split(";")[1].split()[1][1:-1])
            clean_df = clean_overlap_trans(df, trans_read_dict, trans_path_dict)
     
            # order transcript
            ordered_df = order_trans_df(clean_df)
      
            summary_dict = {"gff": ordered_df, "trans_path": trans_path_dict}
            with open(trans_out_obj, "wb") as handle:
                pickle.dump(summary_dict, handle)
