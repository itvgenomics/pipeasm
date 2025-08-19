rule check_hic_files:
    input:
        r1=config['hic_r1'],
        r2=config['hic_r2']
    output:
        "results/Trimming_QC/HiC/{sample}_R1.fastq.gz",
        "results/Trimming_QC/HiC/{sample}_R2.fastq.gz"
    shell:
        """
        mkdir -p results/Trimming_QC/HiC && \
        cp {input.r1} results/Trimming_QC/HiC/{wildcards.sample}_R1.fastq.gz && \
        cp {input.r2} results/Trimming_QC/HiC/{wildcards.sample}_R2.fastq.gz
        """
