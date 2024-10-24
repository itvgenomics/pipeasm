rule check_ont_files:
    input:
        config['ont_reads']
    output:
        f"results/Trimming_QC/ONT/{config['sample']}_ONT.fastq.gz",
    threads: 1
    shell:
        f"""
        mkdir -p results/Trimming_QC/ONT && \
        cp {config['ont_reads']} results/Trimming_QC/ONT/{config["sample"]}_ONT.fastq.gz
        """