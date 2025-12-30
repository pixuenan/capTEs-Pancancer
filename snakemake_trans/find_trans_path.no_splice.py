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

def clean_transcripts(path_dict, dist_thre=30):
    remove_trans_ls = []
    path_gff_ls = [i["gff"] for i in path_dict.values()]
    for idx, trans_df in enumerate(path_gff_ls):
        if 0 in trans_df["start"].tolist() or 0 in trans_df["end"].tolist():
            remove_trans_ls.append(idx)
        if idx == len(path_gff_ls)-1:
            continue
        cur_trans_df = trans_df[trans_df["type"]!="transcript"]
        for n_idx, next_trans_df in enumerate(path_gff_ls[idx+1:]):
            ## remove transcript with identical exon number
            n_trans_df = next_trans_df[next_trans_df["type"]!="transcript"]
            if cur_trans_df.shape[0] == n_trans_df.shape[0]:
                trans_distance = abs(cur_trans_df["start"]-n_trans_df["start"]).sum()+abs(cur_trans_df["end"]-n_trans_df["end"]).sum()
                if trans_distance < dist_thre:
                    remove_trans_ls.append(n_idx+idx+1)

    keep_trans_ls = set(list(range(len(path_gff_ls)))).difference(set(remove_trans_ls))
    keep_trans_df = [path_gff_ls[i] for i in keep_trans_ls]
    clean_path_gff = pd.concat(keep_trans_df).reset_index()
    clean_path_gff.drop(columns=["index"], inplace=True)
    return clean_path_gff

def order_trans_df(df):
    trans_order = df[df["type"]=="transcript"].sort_values(by=["start"])["trans_id"]
    new_trans_df_ls = []
    for trans_id in trans_order:
        new_trans_df_ls.append(df[df["trans_id"]==trans_id])
    ordered_gff = pd.concat(new_trans_df_ls)
    return ordered_gff

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

def trans_extend(to_added_eles, end_trans_df, te_dict, direction, trans_polya=False):
    ele_num = len(to_added_eles)
    for idx, ele in enumerate(to_added_eles):
        end_trans_index = end_trans_df.index
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
    return end_trans_df.reset_index()

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

def get_trans_ele(read_dict, trans_part_ls):
    ele_ls = []
    for part in trans_part_ls:
        if part in read_dict.keys():
            ele_ls.extend(read_dict[part])
        else:
            ele_ls.append(part)
    return set(ele_ls)

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

    def get_cut_read_ls(self, read_content, read_dict, supp_reads):
        self.forw_read_ls = []
        self.reve_read_ls = []
        for read in read_content:
            if read.mapq < 10:
                continue
            if read.is_secondary or read.qname in supp_reads:
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
                       select_te_dict, sample_name, supp_reads):

    read_dict = {}
    te_dict = {}
    remove_idx_ls = []
    ## get reads besides cut site
    for idx, row in select_te_df.iterrows():
        te_obj = TE(row["te_id"])
        cor_ls = select_te_dict[row["te_id"]]
        te_obj.get_cut_coordinate(cor_ls, row["start"], row["end"])
        ## remove TEs that are not target of sgRNA
        if not te_obj.target:
            remove_idx_ls.append(idx)
            continue

        read_content = bam_content.fetch(chr_name, row["start"], row["end"])
        read_dict = te_obj.get_cut_read_ls(read_content, read_dict, supp_reads)
        te_obj.get_end_read_ls(read_dict)
        te_obj.get_cut_range(read_dict)
        te_dict[row["te_id"]] = te_obj
    select_te_df.drop(index=remove_idx_ls, inplace=True)
    select_te_df.reset_index(inplace=True)

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

    ## select TE with end reads
    te_path_ls = []
    for te_id, te_obj in te_dict.items():
        if not te_obj.forw_end_ls and not te_obj.reve_end_ls:
            continue
        if len(te_obj.forw_end_ls) >= 3:
            current_path = [te_id]
            for i in find_path_end(current_path, "forward", te_dict):
                te_path_ls.append(i+["forward"])
        if len(te_obj.reve_end_ls) >= 3:
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
        ## path end with TE
        end_trans = None
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
        new_id = "{}_{}_{}_{}_{}".format(chr_name, region_start, region_end, p_idx, sample_name)
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

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Find link between TE outside splice regions')
    parser.add_argument('-rf', '--remap_fd')
    parser.add_argument('-cr', '--chrom_region')
    parser.add_argument('-gf', '--gff_fd')
    parser.add_argument('-chr', '--chr')
    parser.add_argument('-bam', '--bam_fd')
    parser.add_argument('-s', '--sample_name')
    parser.add_argument('-b', '--te_bed')
    parser.add_argument('-d', '--te_dict')
    args = parser.parse_args()
    bam = "{}/{}.hg38.bam".format(args.bam_fd, args.sample_name)
    trans_out_obj = "{}/{}/obj/{}.trans.no_splice.dict.obj".format(args.remap_fd, args.sample_name, args.chr)
    region_te_obj = "{}/{}/{}.no_splice.te.dict.obj".format(args.gff_fd, args.sample_name, args.sample_name)

    with open(region_te_obj, "rb") as handle:
        region_te_dict = pickle.load(handle)
    with open(args.te_dict, "rb") as handle:
        te_cor_dict = pickle.load(handle)
    with open(args.chrom_region, "rb") as handle:
        chr_region_dict = pickle.load(handle)
    bam_content = pysam.AlignmentFile(bam, "rb")
    te_df = pd.read_csv(args.te_bed, sep="\t", names=["chr", "start", "end", "strand", "repName","repClass","repFamily", "te_id"])
    chr_n, chr_idx = args.chr.split("_")
    chr_region = chr_region_dict[args.chr]
    df_ls = []
    path_dict = {}
    for key, value in region_te_dict.items():
        info_ls = key.split("-")
        if info_ls[0] != chr_n or not value:
            continue
        if int(info_ls[1])>chr_region[1] or int(info_ls[2])<chr_region[0]:
            continue
        chr_name, region_start, region_end = info_ls[0], int(info_ls[1]), int(info_ls[2])
        select_te_df = te_df[te_df["te_id"].isin(value)]
        select_te_dict = {i:te_cor_dict[i] for i in value}
        # record supplementary reads
        supp_reads = []
        read_content = bam_content.fetch(chr_name, region_start, region_end)
        for read in read_content:
            if read.is_supplementary:
                supp_reads.append(read.qname) 

        # link transcript and TE end reads
        te_dict, trans_path_dict = get_candidate_path(chr_name, region_start, region_end, bam_content, select_te_df, 
                                                      select_te_dict, args.sample_name, supp_reads) 
        if trans_path_dict:
            # remove redudant transcripts
            clean_df = clean_transcripts(trans_path_dict)
            df_ls.append(clean_df)
            path_dict.update(trans_path_dict)
    if not df_ls:
        with open(trans_out_obj, "wb") as handle:
            pickle.dump(None, handle)
    else:
        df = pd.concat(df_ls).reset_index()
        df["trans_id"] = df["info"].apply(lambda x: x.split(";")[1].split()[1][1:-1])
        
        # order transcript
        ordered_df = order_trans_df(df)
        
        summary_dict = {"gff": ordered_df, "trans_path": path_dict}
        with open(trans_out_obj, "wb") as handle:
            pickle.dump(summary_dict, handle)
      
      
      
      
      
      
