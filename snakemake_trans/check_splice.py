"""
Script to check reads supporting splice site.
Xuenan Pi
2024.11.26
"""

from multiprocessing.pool import ThreadPool
import pandas as pd
import subprocess
import argparse
import datetime
import pickle
import glob
import re
import os

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
    return {"start": start, "end": end, "splice_ls": splice_ls, "strand": row["strand"]}

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Check reads support splice site within samples.')
    parser.add_argument('-b', '--bed_fold', help="Path to the folder containing sample bed files")
    parser.add_argument('-s', '--splice_fold', help="""Path to the folder to save splice results. 
                                                       Every run of the program will generate result 
                                                       file named {splice_id}.dict.obj in this folder.
                                                       Do not output to or delete any file in this folder.""")
    parser.add_argument('-f', '--splice_file', help="""Path to the file containing splice site to search. 
                                                       Every row contains one splice site id in the format
                                                       of {chr}_{start}_{end}.""")
    parser.add_argument('-o', '--out_prefix', help="""Prefix of the output file. Could be in the format of
                                                       {folder_name}/prefix.""")
    parser.add_argument('-t', '--threads', help="Number of threads")
    args = parser.parse_args()
    current_datetime = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    site_file_ls = glob.glob("{}/*.dict.obj".format(args.splice_fold))
    cal_site_ls = [os.path.basename(i).split(".")[0] for i in site_file_ls]
    with open(args.splice_file, "r+") as handle:
        lines = handle.readlines()
    splice_ls = [i.strip() for i in lines]
    ## remove splice site existed from calculation
    to_cal_ls = []
    bed_cont = []
    for i in splice_ls:
        if i in cal_site_ls:
            continue
        to_cal_ls.append(i)
        bed_cont.append(re.sub("_", "\t", i)+"\n")

    ## calculate splice site supporting read number
    if to_cal_ls:
        bed_out = "{}/{}.bed".format(args.splice_fold, current_datetime)
        with open(bed_out, "w+") as handle:
            handle.writelines(bed_cont)
        file_ls = glob.glob("{}/*.all.sorted.clean.bed".format(args.bed_fold))
        summary_dict = {}
        pool = ThreadPool(int(args.threads))
        processes = []
        sample_ls = []
        for filename in file_ls:
            sample_name = os.path.basename(filename).split(".")[0]
            sample_ls.append(sample_name)
            read_bed_file = "{}/{}.{}.read.bed".format(args.splice_fold, current_datetime, sample_name)
            bed_cmd = "bedtools intersect -a {} -b {} -wb > {}".format(bed_out, filename, read_bed_file)
            processes.append(pool.apply_async(subprocess.Popen(bed_cmd, shell=True)))
        pool.close()
        pool.join()
        for sample_name in sample_ls:
            read_bed_file = "{}/{}.{}.read.bed".format(args.splice_fold, current_datetime, sample_name)
            #subprocess.check_call(bed_cmd, shell=True)
            df = pd.read_csv(read_bed_file, sep="\t", usecols=[0,1,2,4,6,7,8,12,13,14], names=["r_chr", "r_s", "r_e", "start", "read_id", "mapq", "strand", "exon_num", "exon_len", "start_diff"])
            df["region_id"] = df.apply(lambda x: "{}_{}_{}".format(x["r_chr"], x["r_s"], x["r_e"]), axis=1)
            for splice_id in to_cal_ls:
                splice_df = df[df["region_id"] == splice_id]
                if splice_df.empty:
                    summary_dict.setdefault(splice_id, {}).setdefault(sample_name, {})
                    continue
                hit_dict = {}
                for idx, row in splice_df.iterrows():
                    region_start, region_end = row['r_s'], row['r_e']
                    flag = False
                    if row["mapq"] < 10:
                        continue
                    sp_dict = build_read_splice(row)
                    for i in sp_dict["splice_ls"]:
                        if i == [region_start, region_end]:
                            flag = True
                    if flag:
                        hit_dict[row["read_id"]] = sp_dict
                summary_dict.setdefault(splice_id, {}).setdefault(sample_name, hit_dict)
            rm_cmd = "rm {}".format(read_bed_file)
            subprocess.check_call(rm_cmd, shell=True)
        subprocess.check_call("rm {}".format(bed_out), shell=True)
        ## write read information for each splice site into a dictionary object 
        for splice_id, info_dict in summary_dict.items():
            with open("{}/{}.dict.obj".format(args.splice_fold, splice_id), "wb") as handle:
                pickle.dump(info_dict, handle)

    ## summary information for every splice site
    count_dict = {}
    result_dict = {}
    for splice_id in splice_ls:
        with open("{}/{}.dict.obj".format(args.splice_fold, splice_id), "rb") as handle:
            info_dict = pickle.load(handle)
            result_dict[splice_id] = info_dict
            for sample_name, value in info_dict.items():
                count_dict.setdefault(splice_id, {}).setdefault(sample_name, len(value))
    count_df = pd.DataFrame(count_dict).reset_index()
    count_df.sort_values(by="index", inplace=True)
    count_df.rename(columns={"index": "sample_id"}, inplace=True)
    out_count = "{}.count.tsv".format(args.out_prefix)
    count_df.to_csv(out_count, sep="\t", index=False)
    with open("{}.result.dict.obj".format(args.out_prefix), "wb") as handle:
        pickle.dump(result_dict, handle)
