def get_nosplice_region_ls_wd(wildcards):
    with open(config["ref"]["hg38_chrom_region"], "rb") as handle:
        chrom_region_dict = pickle.load(handle)
    subfix_ls = chrom_region_dict.keys()
    return [config["path"]["remapping"]+"{}/obj/{}.trans.no_splice.dict.obj".format(wildcards.sample_name, i) for i in subfix_ls]

rule no_splice_region:
    input:
        config["path"]["gff"]+"{sample_name}/{sample_name}.R.flair.trans.ex20bp.selected.merged.bed"
    params:
        config["ref"]["hg38_chrom"]
    output:
        config["path"]["gff"]+"{sample_name}/{sample_name}.R.flair.no_splice.ex20bp.bed"
    shell:
        """
        bedtools subtract -a {params} -b {input} > {output}
        """

rule no_splice_te:
    input:
        config["path"]["gff"]+"{sample_name}/{sample_name}.R.flair.no_splice.ex20bp.bed"
    params:
        config["ref"]["te_bed"],
    output:
        config["path"]["gff"]+"{sample_name}/{sample_name}.R.flair.no_splice.ex20bp.TE.bed"
    shell:
        """
        bedtools intersect -a {input} -b {params} -wa -wb > {output}
        """

rule no_splice_te_dict:
    input:
        config["path"]["gff"]+"{sample_name}/{sample_name}.R.flair.no_splice.ex20bp.TE.bed"
    params:
        config["path"]["gff"]+"{sample_name}/{sample_name}.all.read.TE.mq10.exon.bed"
    output:
        config["path"]["gff"]+"{sample_name}/{sample_name}.no_splice.te.dict.obj"
    run:
        te_read_df = pd.read_csv(params[0], sep="\t", usecols=[3,11], names=["read_id", "te_id"])
        te_read_dict = {}
        for idx, row in te_read_df.iterrows():
            te_read_dict.setdefault(row["te_id"], 1)
        region_te_df = pd.read_csv(input[0], sep="\t", usecols=[0,1,2,10], names=["chr", "start", "end", "te_id"])
        region_te_df["region_id"] = region_te_df.apply(lambda x: "{}-{}-{}".format(x["chr"], x["start"], x["end"]), axis=1)
        grouped = region_te_df.groupby("region_id")
        region_te_dict = {}
        for group_name, group in grouped:
            te_ls = list(set(group['te_id']))
            for te_id in te_ls:
                if te_id in te_read_dict:
                    region_te_dict.setdefault(group_name, []).append(te_id)
        with open(output[0], "wb") as handle:
            pickle.dump(region_te_dict, handle)

rule no_splice_gff:
    input:
        config["path"]["gff"]+"{sample_name}/{sample_name}.no_splice.te.dict.obj"
    params:
        chr="{chr}",
        s="{sample_name}",
        b=config["ref"]["te_bed"],
        d=config["ref"]["te_cor_dict"],
        cr=config["ref"]["hg38_chrom_region"],
        rf=config["path"]["remapping"],
        gf=config["path"]["gff"],
        bam="/NAS/wg_pxn/share/TEcap/04.bam/",
        script=config["tool"]["find_path_no_splice"]
    output:
        config["path"]["remapping"]+"{sample_name}/obj/{chr}.trans.no_splice.dict.obj",
    shell:
        """
        python {params.script} -rf {params.rf} -cr {params.cr} -gf {params.gf} -s {params.s} -chr {params.chr} -bam {params.bam} -b {params.b} -d {params.d}
        """

rule link_no_splice_gff:
    input:
        file_ls=get_nosplice_region_ls_wd
    output:
        config["path"]["gff"]+"{sample_name}/{sample_name}.link.sr.clean.no_splice.gff"
    run:
        trans_df_ls = []
        for filename in input:
            with open(filename, "rb") as handle:
                trans_dict = pickle.load(handle)
            if not trans_dict:
                continue
            df = trans_dict['gff']
            rm_ls = []
            trans_ls = df["trans_id"].unique().tolist()
            for trans_id in trans_ls:
                trans_df = df[df["trans_id"]==trans_id] 
                trans_df_ls.append(trans_df)

        ## rename information column
        df_ls = []
        for trans_df in trans_df_ls:
            trans_df["gene_id"] = trans_df["trans_id"].apply(lambda x: "_".join(x.split("_")[:-1]))
            ## drop transcripts with one exon
            new_info_ls = []
            trans_idx = trans_df.index[0]
            for idx, row in trans_df.iterrows():
                if row["type"] == "transcript":
                    new_info = 'gene_id "{}"; transcript_id "{}";'.format(row["gene_id"], row["trans_id"])
                else:
                    new_info = 'gene_id "{}"; transcript_id "{}"; exon_number "{}";'.format(row["gene_id"], row["trans_id"], idx-trans_idx)
                new_info_ls.append(new_info)
            trans_df["new_info"] = new_info_ls
            df_ls.append(trans_df)

        rename_df = pd.concat(df_ls)
        rename_df.drop(columns=["info"], inplace=True)
        rename_df.rename(columns={"new_info": "info"}, inplace=True)
        remove_ls = rename_df.loc[rename_df["start"]>=rename_df["end"], "trans_id"].unique().tolist()
        rename_df[2] = '.'
        rename_df.loc[~rename_df["trans_id"].isin(remove_ls), gff_columns].to_csv(output[0], sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE)

