rule mapping:
    input:
        config["path"]["qc_fastq"]+"{sample_name}/basecall_pass.clean.fq"
    output:
        sam=temp(config["path"]["bam"]+"{sample_name}.hg38.sam"),
        bam=config["path"]["bam"]+"{sample_name}.hg38.bam"
    params:
        config["ref"]["hg38_fa"]
    threads: 20
    run:
        shell("""
        minimap2 -ax splice --MD {params} {input} -t {threads} > {output.sam}
        samtools sort -@ {threads} {output.sam} -o {output.bam}
        samtools index {output.bam}
        """)

rule bam2fq:
    input:
        config["path"]["bam"]+"{sample_name}.hg38.bam"
    output:
        temp(config["path"]["bam"]+"fq/{sample_name}.fq")
    threads: 20
    shell:
        """
        samtools fastq -@ {threads} {input} > {output}
        """
