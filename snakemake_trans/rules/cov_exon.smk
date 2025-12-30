rule read_end:
    input:
        cor=config["path"]["gff"]+"{sample_name}/{sample_name}_all_corrected.bed",
        inconsi=config["path"]["gff"]+"{sample_name}/{sample_name}_all_inconsistent.bed"
    params:
        config["path"]["gff"]+"{sample_name}/{sample_name}_all.bed"
    output:
        config["path"]["gff"]+"{sample_name}/{sample_name}_all.read.mq10.bed"
    run:
        shell("""
        cat {input.cor} {input.inconsi} | bedtools sort -i - | uniq > {params}
        awk -F'\t' '{{if($5>=10){{print}}}}' {params} > {output}
        rm {params}
        """)

rule exon_bed:
    input:
        config["path"]["gff"]+"{sample_name}/{sample_name}_all.read.mq10.bed"
    output:
        temp(config["path"]["gff"]+"{sample_name}/{sample_name}_all.read.mq10.exon.bed")
    run:
        df = pd.read_csv(input[0], sep="\t", names=["chr", "start", "end", "read_id", "mapq", "strand", "s", "e", "info", "exon_num", "exon_len", "start_diff"])
        new_row_ls = []
        for idx, row in df.iterrows():
            new_row_ls.extend(construct_exon(row))
        exon_df = pd.DataFrame(new_row_ls)
	exon_df.to_csv(output[0], sep="\t", index=False, header=False)

rule te_read:
    input:
        config["path"]["gff"]+"{sample_name}/{sample_name}_all.read.mq10.exon.bed"
    output:
        config["path"]["gff"]+"{sample_name}/{sample_name}.all.read.TE.mq10.exon.bed"
    params:
        config['ref']['te_bed']
    run:
        shell("""
        bedtools intersect -a {input} -b {params} -wb > {output}
        """)

rule summary_read:
    input:
        te_read=config["path"]["gff"]+"{sample_name}/{sample_name}.all.read.TE.mq10.exon.bed",
    output:
        config["path"]["gff"]+"{sample_name}/{sample_name}.TE.mq10.exon.count.summary.csv"
    params: 
        bed=config['ref']['te_bed'],
        s_name="{sample_name}"
    run:
        df = pd.read_csv(input.te_read, sep="\t", usecols=[3,11], names=["read_id", "te_id"])
        te_dict = {}
        for idx, row in df.iterrows():
            te_dict.setdefault(row["te_id"], []).append(row["read_id"])
        te_read_count_dict = {}
        for key, value in te_dict.items():
            te_read_count_dict[key] = len(set(value))
        
        te_df = pd.read_csv(params.bed, sep="\t", names=["chr", "start", "end", "strand", "repName","repClass","repFamily", "te_id"])
        te_count_df = pd.DataFrame({"count": te_read_count_dict}).reset_index().rename(columns={"index": "te_id"})
        final_df = te_count_df.merge(te_df, how="left")
        final_df[["te_id", "count"]].to_csv(output[0], index=False)
