rule check_hifi_files:
    input:
        config['hifi_reads']
    output:
        f"results/Trimming_QC/HiFi/{config['sample']}.fastq.gz",
    threads: 1
    shell:
        f"""
        mkdir -p results/Trimming_QC/HiFi && \
        cp {config['hifi_reads']} results/Trimming_QC/HiFi/{config["sample"]}.fastq.gz
        """
