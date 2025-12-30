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

seq2_aligner = Align.PairwiseAligner()
seq2_aligner.target_internal_gap_score = -2
seq2_aligner.target_internal_extend_gap_score = -5
seq2_aligner.mismatch_score = -1
seq2_aligner.mode = "local"

gff_columns = ["chr", 1, "type", "start", "end", 2, "strand", 3, "info"]

def get_z_score_read(read_df, trans_p_dict, read_id, trans_id):
    indicator = read_df.loc[read_id, trans_id]
    denominator = sum([read_df.loc[read_id, t_id]*trans_p_dict[t_id] for t_id in trans_p_dict.keys()])
    return trans_p_dict[trans_id]*indicator/denominator

def get_trans_p(trans_read_z_dict, read_number, trans_id):
    return sum(trans_read_z_dict[trans_id].values())/read_number

def get_trans_sreads(read_df):
    """Apply EM algorithm on read transcript matrix to get reads supporting transcript"""
    trans_read_z_dict = {}
    read_num, trans_num = read_df.shape
    trans_p_dict = {i: 1/trans_num for i in read_df.columns.values}
    ## speed calculation by clustering reads
    db = DBSCAN(eps=0.2, min_samples=2).fit(read_df)
    read_idx_dict = {value: key for key, value in read_df.reset_index()["index"].to_dict().items()}

    for i in range(200):
        dup_dict = {}
        cur_trans_p_dict = {key: value for key, value in trans_p_dict.items()}
        for read_id in read_df.index.tolist():
            cluster_num = db.labels_[read_idx_dict[read_id]]
            if not cluster_num in dup_dict.keys():
                for trans_id in read_df.columns.values:
                    trans_read_z_dict.setdefault(trans_id, {}).update({read_id: get_z_score_read(read_df, trans_p_dict, read_id, trans_id)})
                if cluster_num != -1:
                    dup_dict[cluster_num] = read_id
            else:
                cal_read = dup_dict[cluster_num]
                for trans_id in read_df.columns.values:
                    trans_read_z_dict.setdefault(trans_id, {}).update({read_id: trans_read_z_dict[trans_id][cal_read]})
        for trans_id in read_df.columns.values:
            trans_p_dict[trans_id] = get_trans_p(trans_read_z_dict, read_num, trans_id)
        diff = sum([abs(cur_trans_p_dict[i]-trans_p_dict[i]) for i in trans_p_dict.keys()])
        if diff < 1e-6:
            break
    trans_read_z_df = pd.DataFrame(trans_read_z_dict)
    trans_ls = trans_read_z_df.columns.tolist()
    return trans_read_z_df

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
           and idx!=ele_num-1 and to_added_eles[idx+1][:2]!="transcript":
                continue
        elif "_".join(ele.split("_")[:-1]) in ["Alu", "long_TE"]:
            te_id = ele
            ## foward direction, extend the first exon
            if direction == "forward":
                te_start = trans_polya and te_dict[te_id].reve_end_start or te_dict[te_id].reve_poly_start
                end_trans_df.loc[end_trans_index[0], "start"] = te_start
                end_trans_df.loc[end_trans_index[1], "start"] = te_start
            ## reverse direction, extend the last exon
            elif direction == "reverse":
                te_end = trans_polya and te_dict[te_id].forw_end_end or te_dict[te_id].forw_poly_end
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

    def get_cut_read_ls(self, read_content, read_dict):
        self.forw_read_ls = []
        self.reve_read_ls = []
        for read in read_content:
            if read.mapq < 10:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
            ## search for reads start within sgRNA match region
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
        self.forw_end_end = forw_end_ls and np.bincount(forw_end_ls).argmax() or 0
        self.reve_end_start = reve_start_ls and np.bincount(reve_start_ls).argmax() or 0
        forw_poly_end_ls = [read_dict[i].end for i in self.forw_poly_read]
        reve_poly_start_ls = [read_dict[i].start for i in self.reve_poly_read]
        self.forw_poly_end = forw_poly_end_ls and np.bincount(forw_poly_end_ls).argmax() or 0
        self.reve_poly_start = reve_poly_start_ls and np.bincount(reve_poly_start_ls).argmax() or 0

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

    def clean_polyA_in_end_reads(self, forw_ls, reve_ls):
        """exclude polyA reads in end reads"""
        forw_polya_reads = set(self.clean_forw_end_ls).intersection(set(forw_ls))
        reve_polya_reads = set(self.clean_reve_end_ls).intersection(set(reve_ls))
        forw_polya_count, reve_polya_count = len(forw_polya_reads), len(reve_polya_reads)
        self.clean_forw_end_ls = list(set(self.clean_forw_end_ls).difference(set(forw_ls)))
        self.clean_reve_end_ls = list(set(self.clean_reve_end_ls).difference(set(reve_ls)))
        self.forw_poly_count = forw_polya_count
        self.reve_poly_count = reve_polya_count
        self.forw_poly_read = forw_polya_reads
        self.reve_poly_read = reve_polya_reads

def get_candidate_path(chr_name, region_start, region_end, bam_content, select_te_df,
                       select_te_dict, clean_select_gff, keep_trans_dict, sample_name):
    ## extract read with poly A at the soft clipping part with poly length exceeds 20 bp plus seq2 primer
    global for_ls, rev_ls
    for_ls, rev_ls = [], []
    for read in bam_content.fetch(chr_name, region_start, region_end):
        if read.mapq < 10:
            continue
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
        te_obj.get_cut_range(read_dict)
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
            elif trans_id in te_obj.reve_link_trans:
                if (reve_polya_count >= 2 and reve_polya_count/len(trans_info[1]) >= 0.33) or reve_polya_count >= 5:
                    candidate_trans_polya_path.append([te_id, trans_id, "reverse"])
                    trans_polya_ls.append(trans_id)
            elif trans_id in te_obj.forw_link_trans:
                if (forw_polya_count >= 2 and forw_polya_count/len(trans_info[1]) >= 0.33) or forw_polya_count >= 5:
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
    trans_out_obj = "{}/{}/obj/{}.trans.dict.obj".format(args.remap_fd, args.sample_name, args.region_id)
    info_ls = args.region_id.split("-")
    chr_name, region_start, region_end = info_ls[0], int(info_ls[1]), int(info_ls[2])

    with open(stringtie_gff_obj, "rb") as handle:
        stringtie_gff = pickle.load(handle)
    with open(flair_gff_obj, "rb") as handle:
        flair_gff = pickle.load(handle)
    with open(inconsis_gff_obj, "rb") as handle:
        inconsis_gff = pickle.load(handle)
    flair_select_gff = get_trans_df(flair_gff, chr_name, region_start, region_end)
    stringtie_select_gff = get_trans_df(stringtie_gff, chr_name, region_start, region_end)
    inconsis_select_gff = get_trans_df(inconsis_gff, chr_name, region_start, region_end)
    select_gff = pd.concat([stringtie_select_gff, flair_select_gff, inconsis_select_gff])
    ## build read splice information 
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
    trans_read_dict = {}
    for trans_id in trans_ls:
        trans_df = select_gff.loc[(select_gff["trans_id"]==trans_id) & (select_gff["type"]=="exon"), ["start", "end"]].reset_index()
        exon_num = trans_df.shape[0]
        if not exon_num in read_dict.keys():
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
        if len(value) < 3:
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
     
        ## link transcript and TE end reads
        bam_content = pysam.AlignmentFile(bam, "rb")
        with open(args.te_dict, "rb") as handle:
            te_cor_dict = pickle.load(handle)
        te_df = pd.read_csv(args.te_bed, sep="\t", names=["chr", "start", "end", "strand", "repName","repClass","repFamily", "te_id"])
 
        select_te_df = te_df[(te_df["chr"]==chr_name) & (te_df["start"]>=region_start) & (te_df["end"]<=region_end)]
        select_te_dict = {i:te_cor_dict[i] for i in select_te_df["te_id"]}
        te_dict, trans_path_dict = get_candidate_path(chr_name, region_start, region_end, bam_content, select_te_df, select_te_dict, clean_select_gff, keep_trans_dict, args.sample_name) 
 
        if not trans_path_dict:
            with open(trans_out_obj, "wb") as handle:
                pickle.dump(None, handle)
        else:
            ## remove redudant transcripts
            df = clean_gff(trans_path_dict)
            df["trans_id"] = df["info"].apply(lambda x: x.split(";")[1].split()[1][1:-1])
            clean_df = clean_overlap_trans(df, trans_read_dict, trans_path_dict)
     
            ## order transcript
            ordered_df = order_trans_df(clean_df)
      
            ## em estimation
            trans_ls = clean_df['trans_id'].unique().tolist()
            read_df, read_dict = cal_trans_read_dict(trans_ls, trans_path_dict, te_dict, trans_read_dict)
            
            multi_num = 1
            ## downsample reads
            if read_df.shape[0] >= 40000:
               read_ls = read_df.index.tolist()
               rng = np.random.default_rng()
               new_read_ls = rng.choice(read_ls, 20000, replace=False)
               multi_num = read_df.shape[0]/20000
               read_df = read_df.loc[new_read_ls]
               
            if read_df.shape[1] >= 300:
                top_300_trans = read_df.sum().sort_values(ascending=False).index[:300]
                read_sum = read_df[top_300_trans].sum(axis=1)
                remove_read = read_sum[read_sum==0].index
                keep_read = read_df.index.difference(remove_read)
                read_df = read_df.loc[keep_read, top_300_trans]
            trans_read_df = get_trans_sreads(read_df)
            summary_dict = {"trans_reads_df": trans_read_df, "read_info": read_dict, "gff": ordered_df, "trans_path": trans_path_dict, "trans_read": trans_read_dict, "multi_num": multi_num}
            with open(trans_out_obj, "wb") as handle:
                pickle.dump(summary_dict, handle)
      
      
      
      
      
      
