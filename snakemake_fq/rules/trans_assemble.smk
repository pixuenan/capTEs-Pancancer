rule stringtie:
    input:
        ancient(config["path"]["bam"]+"{sample_name}.hg38.bam")
    output:
        config["path"]["gff"]+"{sample_name}/{sample_name}.hg38.R.gff"
    params:
        config["ref"]["hg38_gff"]
    threads: 20
    run:
        shell("""
        stringtie -R -p {threads} -G {params} -o {output} {input}
        """)

rule convert2bed:
    input:
        config["path"]["bam"]+"{sample_name}.hg38.bam"
    output:
        temp(config["path"]["bam"]+"{sample_name}.hg38.bed")
    run:
        shell("""
        bam2Bed12 -i {input} > {output} 
        """)

rule build_trans:
    input:
        bed=rules.convert2bed.output,
        fq=config["path"]["bam"]+"fq/{sample_name}.fq"
    output:
        correct=config["path"]["gff"]+"{sample_name}/{sample_name}.isoforms.gtf",
        inconsist=config["path"]["gff"]+"{sample_name}/{sample_name}.inconsist.isoforms.gtf"
    params:
        ref_gff=config["ref"]["hg38_gff"],
        ref_fa=config["ref"]["hg38_fa"],
        prefix=config["path"]["gff"]+"{sample_name}/{sample_name}"
    threads: 20
    run:
        shell("""
        flair correct -q {input.bed} -f {params.ref_gff} -g {params.ref_fa} --output {params.prefix} --threads {threads}
        flair collapse -q {params.prefix}_all_corrected.bed -r {input.fq} -g {params.ref_fa} --filter comprehensive --support 2 --gtf {params.ref_gff} --output {params.prefix} --threads {threads} 
        flair collapse -q {params.prefix}_all_inconsistent.bed -r {input.fq} -g {params.ref_fa} --filter comprehensive --gtf {params.ref_gff} --output {params.prefix}.inconsist --threads {threads} 
        """)

rule convert_obj:
    input:
        flair=config["path"]["gff"]+"{sample_name}/{sample_name}.isoforms.gtf",
        stringtie=config["path"]["gff"]+"{sample_name}/{sample_name}.hg38.R.gff",
        inconsist=config["path"]["gff"]+"{sample_name}/{sample_name}.inconsist.isoforms.gtf"
    output:
        flair=config["path"]["gff"]+"{sample_name}/{sample_name}.flair.gff.obj",
        stringtie=config["path"]["gff"]+"{sample_name}/{sample_name}.R.gff.obj",
        inconsist=config["path"]["gff"]+"{sample_name}/{sample_name}.inconsist.gff.obj"
    run:
        gff_columns = ["chr", 1, "type", "start", "end", 2, "strand", 3, "info"]
        flair_gff = pd.read_csv(input.flair, sep="\t", names=gff_columns)
        flair_gff["trans_id"] = flair_gff["info"].apply(lambda x: x.split(";")[1].split()[1][1:-1])
        stringtie_gff = pd.read_csv(input.stringtie, sep="\t", names=gff_columns, skiprows=2)
        stringtie_gff["trans_id"] = stringtie_gff["info"].apply(lambda x: x.split(";")[1].split()[1][1:-1])
        inconsist_gff = pd.read_csv(input.inconsist, sep="\t", names=gff_columns)
        inconsist_gff["trans_id"] = inconsist_gff["info"].apply(lambda x: x.split(";")[1].split()[1][1:-1])
        with open(output.flair, "wb") as handle:
            pickle.dump(flair_gff, handle)
        with open(output.stringtie, "wb") as handle:
            pickle.dump(stringtie_gff, handle)
        with open(output.inconsist, "wb") as handle:
            pickle.dump(inconsist_gff, handle)

