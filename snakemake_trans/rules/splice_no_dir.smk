def get_nodir_region_ls_wd(wildcards):
    inter_df = pd.read_csv("{}/{}/{}.R.flair.trans.ex20bp.selected.merged.bed".format(config["path"]["gff"], wildcards.sample_name, wildcards.sample_name), sep="\t", names=["r_chr", "r_start", "r_end", "trans_id_ls"])
    inter_df["region"] = inter_df.apply(lambda x: "{}-{}-{}".format(x["r_chr"], x["r_start"], x["r_end"]), axis=1)
    return ["{}{}/obj/{}.trans.no_dir.dict.obj".format(config["path"]["remapping"], wildcards.sample_name, i) for i in inter_df["region"].tolist()]

rule region_process_no_dir:
    input:
        split=config["path"]["remapping"]+"{sample_name}/split_region.txt",
        f_bed=config["path"]["gff"]+"{sample_name}/{sample_name}.all.sorted.clean.bed",
    params:
        region="{region}",
        s="{sample_name}",
        b=config["ref"]["te_bed"],
        d=config["ref"]["te_cor_dict"],
        rf=config["path"]["remapping"],
        gf=config["path"]["gff"],
        script=config["tool"]["find_path_no_dir"]
    output:
        trans_dict=config["path"]["remapping"]+"{sample_name}/obj/{region}.trans.no_dir.dict.obj",
    shell:
        """
        python {params.script} -rf {params.rf} -gf {params.gf} -s {params.s} -r {params.region} -b {params.b} -d {params.d}
        """

rule link_gff:
    input:
        file_ls=get_nodir_region_ls_wd
    output:
        config["path"]["gff"]+"{sample_name}/{sample_name}.link.sr.clean.no_dir.gff"
    run:
        trans_df_ls = []
        for filename in input:
            with open(filename, "rb") as handle:
                trans_dict = pickle.load(handle)
            if not trans_dict:
                continue
            df = trans_dict['gff']
            rm_ls = []
            if len(df["chr"].unique().tolist()) > 1:
                for idx, row in df.iterrows():
                    if row["chr"] != row["trans_id"].split("_")[0]:
                        rm_ls.append(idx)
                df.drop(rm_ls, axis=0, inplace=True)
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

