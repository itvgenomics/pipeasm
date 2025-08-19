rule check_hifi_files:
    input:
        config['hifi_reads']
    output:
        "results/Trimming_QC/HiFi/{sample}.fastq.gz"
    threads: 1
    shell:
        """
        mkdir -p results/Trimming_QC/HiFi && \
        cp {input} {output}
        """
