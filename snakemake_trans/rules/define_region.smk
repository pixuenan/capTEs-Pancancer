rule expand_trans:
    input:
        config["path"]["gff"]+"{sample_name}/{sample_name}.R.gff.obj"
    output:
        config["path"]["gff"]+"{sample_name}/{sample_name}.R.trans.ex20bp.bed"
    run:
        with open(input[0], "rb") as handle:
            gff_df = pickle.load(handle)
        trans_df = gff_df[gff_df["type"]=="transcript"]
        trans_df.sort_values(by=["chr", "start"], inplace=True)
        trans_df["e_start"] = trans_df["start"].apply(lambda x: x-20 > 0 and x-20 or 1)
        trans_df["e_end"] = trans_df["end"]+20
        trans_df[["chr", "e_start", "e_end", "trans_id"]].to_csv(output[0], sep="\t", index=False, header=False)

rule merge_trans:
    input:
        config["path"]["gff"]+"{sample_name}/{sample_name}.R.trans.ex20bp.bed"
    output:
        config["path"]["gff"]+"{sample_name}/{sample_name}.R.trans.ex20bp.merged.bed"
    shell:
        "bedtools merge -c 4 -o collapse -i {input} > {output}"

rule select_trans:
    input:
        gff_obj=config["path"]["gff"]+"{sample_name}/{sample_name}.R.gff.obj",
        bed=config["path"]["gff"]+"{sample_name}/{sample_name}.R.trans.ex20bp.merged.bed"
    output:
        config["path"]["gff"]+"{sample_name}/{sample_name}.R.trans.ex20bp.merged.selected.bed"
    run:
        with open(input.gff_obj, "rb") as handle:
            gff_df = pickle.load(handle)
        trans_exon_dict = count_exon_number(gff_df)
        bed_df = pd.read_csv(input.bed, sep="\t", names=["chr", "start", "end", "trans_id"])
        bed_df["exon_num"] = bed_df["trans_id"].apply(lambda x: ",".join([str(trans_exon_dict[i]) for i in x.split(",")]))
        bed_df["valid"] = bed_df["exon_num"].apply(lambda x: set(x.split(","))!={"1"} and True or False)
        bed_df.loc[bed_df["valid"]==True, ["chr", "start", "end", "trans_id"]].to_csv(output[0], sep="\t", index=False, header=False)

rule select_trans_flair:
    input:
        gff_obj=config["path"]["gff"]+"{sample_name}/{sample_name}.flair.gff.obj",
        bed=config["path"]["gff"]+"{sample_name}/{sample_name}.flair.trans.ex20bp.merged.bed"
    output:
        config["path"]["gff"]+"{sample_name}/{sample_name}.flair.trans.ex20bp.merged.selected.bed"
    run:
        with open(input.gff_obj, "rb") as handle:
            gff_df = pickle.load(handle)
        trans_exon_dict = count_exon_number(gff_df)
        bed_df = pd.read_csv(input.bed, sep="\t", names=["chr", "start", "end", "trans_id"])
        bed_df["exon_num"] = bed_df["trans_id"].apply(lambda x: ",".join([str(trans_exon_dict[i]) for i in x.split(",")]))
        bed_df["valid"] = bed_df["exon_num"].apply(lambda x: set(x.split(","))!={"1"} and True or False)
        bed_df.loc[bed_df["valid"]==True, ["chr", "start", "end", "trans_id"]].to_csv(output[0], sep="\t", index=False, header=False)
