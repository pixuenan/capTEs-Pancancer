rule merge_bed:
    input:
        cor=config["path"]["gff"]+"{sample_name}/{sample_name}_all_corrected.bed",
        inconsis=config["path"]["gff"]+"{sample_name}/{sample_name}_all_inconsistent.bed"
    output:
        config["path"]["gff"]+"{sample_name}/{sample_name}.all.sorted.clean.bed"
    run:
        shell("""
        cat {input.cor} {input.inconsis} | awk -F'\t' '{{if($10!=1){{print}}}}' | bedtools sort -i - > {output}
        """)

rule split_region:
    input:
        r_bed=config["path"]["gff"]+"{sample_name}/{sample_name}.R.flair.trans.ex20bp.selected.merged.bed",
        bed=config["path"]["gff"]+"{sample_name}/{sample_name}.all.sorted.clean.bed",
    output:
        config["path"]["remapping"]+"{sample_name}/split_region.txt"
    params:
        config["path"]["remapping"]+"{sample_name}"
    run:
        merge_df = pd.read_csv(input.r_bed, sep="\t", names=["r_chr", "r_start", "r_end"])
        merge_df.drop_duplicates(inplace=True)
        if not os.path.exists("{}/bed".format(params[0])):
            os.mkdir("{}/bed".format(params[0]))
        for idx, row in merge_df.iterrows():
            with open("{}/bed/{}-{}-{}.bed".format(params[0], row["r_chr"], row["r_start"], row["r_end"]), "w+") as handle:
                handle.write("{}\t{}\t{}".format(row["r_chr"], row["r_start"], row["r_end"]))
        with open("{}/split_region.txt".format(params[0]), "w+") as handle:
            handle.write("Done.")

rule region_process:
    input:
        bam=config["path"]["bam"]+"{sample_name}.hg38.bam",
        split=config["path"]["remapping"]+"{sample_name}/split_region.txt",
        f_bed=config["path"]["gff"]+"{sample_name}/{sample_name}.all.sorted.clean.bed",
    params:
        bed=config["path"]["remapping"]+"{sample_name}/bed/{region}.bed",
        folder=config["path"]["remapping"]+"{sample_name}/bam",
        region="{region}",
        s="{sample_name}",
        b=config["ref"]["te_bed"],
        d=config["ref"]["te_cor_dict"],
        rf=config["path"]["remapping"],
        gf=config["path"]["gff"],
        script=config["tool"]["find_path"]
    output:
        trans_dict=config["path"]["remapping"]+"{sample_name}/obj/{region}.trans.dict.obj",
        bam=config["path"]["remapping"]+"{sample_name}/bam/{region}.spliced.bam",
        bed=config["path"]["remapping"]+"{sample_name}/bed/{region}.read.bed",
    shell:
        """
        bedtools intersect -a {params.bed} -b {input.f_bed} -wb > {output.bed}
        region=$(echo {params.region} | awk -F'-' '{{print $1":"$2"-"$3}}')
        samtools view -bh {input.bam} $region > {output.bam}
        samtools index {output.bam}
        python {params.script} -rf {params.rf} -gf {params.gf} -s {params.s} -r {params.region} -b {params.b} -d {params.d}
        """
