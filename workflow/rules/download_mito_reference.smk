rule download_mito_reference:
    output:
        fasta="resources/{sample}.reference.fasta",
        gb="resources/{sample}.reference.gb"
    threads: 1
    log:
        "logs/{sample}.download_mito_reference.log"
    benchmark:
        "benchmarks/{sample}.download_mito_reference.txt"
    singularity:
        f"{config["sif_dir"]}/mitohifi.sif"
    shell:
        """
        findMitoReference.py --species '{config[species]}' \
        --email ncbiapirunner@gmail.com \
        --outfolder resources/ \
        --min_length 16000 >> {log} 2>&1 && \
        mv resources/*.fasta {output.fasta} && \
        mv resources/*.gb {output.gb}
        """
