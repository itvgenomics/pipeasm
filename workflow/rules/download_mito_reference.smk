rule download_mito_reference:
    output:
        fasta="resources/{sample}.reference.fasta",
        gb="resources/{sample}.reference.gb"
    log:
        "logs/{sample}.download_mito_reference.log"
    benchmark:
        "benchmarks/{sample}.download_mito_reference.txt"
    singularity:
        f"{config["sif_dir"]}/mitohifi.sif"
    shell:
        """
        rm -f resources/*.fasta resources/*.gb && \
        findMitoReference.py --species '{config[species]}' \
        --email ncbiapirunner@gmail.com \
        --outfolder resources/ \
        --min_length 12000 >> {log} 2>&1 && \
        mv resources/*.fasta {output.fasta} && \
        mv resources/*.gb {output.gb}
        """
