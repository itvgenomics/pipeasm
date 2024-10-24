rule check_hic_files:
    input:
        config['hic_r1'],
        config['hic_r2']
    output:
        "results/Trimming_QC/HiC/{sample}_R1.fastq.gz",
        "results/Trimming_QC/HiC/{sample}_R2.fastq.gz"
    threads: 1
    shell:
        f"""
        mkdir -p results/Trimming_QC/HiC && \
        cp {config['hic_r1']} results/Trimming_QC/HiC/{config["sample"]}_R1.fastq.gz && \
        cp {config['hic_r2']} results/Trimming_QC/HiC/{config["sample"]}_R2.fastq.gz
        """