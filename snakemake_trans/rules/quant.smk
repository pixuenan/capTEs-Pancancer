rule summary_EM:
    input:
        config["path"]["gff"]+"{sample_name}/{sample_name}.link.sr.clean.gff",
    output:
        gff=config["path"]["gff"]+"{sample_name}/{sample_name}.link.sr.clean.EM.gff",
        csv=config["path"]["gff"]+"{sample_name}/{sample_name}.link.sr.clean.EM.csv",
    run:
        gff_df = pd.read_csv(input[0], sep="\t", names=gff_columns)
        gff_df["trans_id"] = gff_df["info"].apply(lambda x: x.split(";")[1].split()[1][1:-1])
        gff_df["region_id"] = gff_df["trans_id"].apply(lambda x: "-".join(x.split("_")[:3]))
        sample_id = os.path.basename(input[0]).split(".")[0]
        region_ls = gff_df["region_id"].unique().tolist()
        df_ls = []
        count_dict = {}
        for region_id in region_ls:
            with open("07.remapping/{}/obj/{}.trans.dict.obj".format(sample_id, region_id), "rb") as handle:
                summary_dict = pickle.load(handle)

            trans_read_z_df = summary_dict['trans_reads_df']
            read_dict = summary_dict['read_info']

            trans_ls = trans_read_z_df.columns.tolist()
            unsp_set = set(read_dict["unspliced"])
            sp_set = set(read_dict["spliced"])
            trans_count_dict = {}
            for idx, row in trans_read_z_df.iterrows():
                donor_trans = trans_ls[row.argmax()]
                if row[donor_trans] < 0.2:
                    continue
                if idx in unsp_set:
                    trans_count_dict.setdefault(donor_trans, {}).setdefault("unspliced", []).append(idx)
                elif idx in sp_set:
                    trans_count_dict.setdefault(donor_trans, {}).setdefault("spliced", []).append(idx)
            for key, value in trans_count_dict.items():
                if not "spliced" in value or not value["spliced"]:
                    continue
                count_dict.setdefault(key, {}).update({"unspliced": "unspliced" in value and len(value["unspliced"]) or 0, 
                                                       "spliced": len(value["spliced"])})
        gff_df.loc[gff_df["trans_id"].isin(count_dict.keys()), gff_columns].to_csv(output.gff, sep="\t", header=False, quoting=csv.QUOTE_NONE, index=False)
        pd.DataFrame(count_dict).T.to_csv(output.csv)


