rule check_hifi_files:
    input:
        config['hifi_reads']
    output:
        temp("results/Trimming_QC/HiFi/{sample}.fastq.gz")
    shell:
        """
        mkdir -p results/Trimming_QC/HiFi && \
        cp {input} {output}
        """
