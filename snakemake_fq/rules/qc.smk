rule concate:
    input:
        get_fq_path_m
    params:
        get_fq_in_string
    output:
        temp(config["path"]["qc_fastq"]+"{sample_name}/basecall_pass.fq.gz")
    shell:
        "cat {params} > {output}" 

rule filter:
    input:
        rules.concate.output
    output:
        temp(config["path"]["qc_fastq"]+"{sample_name}/basecall_pass.clean.fq.gz")
    params:
        config["params"]["nanofilt"]
    threads: 20
    run:
        shell("""
        porechop --extra_end_trim 0 -i {input} -t {threads} | NanoFilt {params} | gzip > {output}
        """)
