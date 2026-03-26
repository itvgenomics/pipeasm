rule check_ont_files:
    input:
        config['ont_reads']
    output:
        temp("results/Trimming_QC/ONT/{sample}_ONT.fastq.gz"),
    shell:
        """
        mkdir -p results/Trimming_QC/ONT && \
        cp {input} {output}
        """
