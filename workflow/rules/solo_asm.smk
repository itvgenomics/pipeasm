rule solo_assembly:
    input:
        hifi="results/Trimming_QC/HiFi/{sample}.trimmed.fastq.gz"
    output:
        "results/Assembly/Contigging/Solo_Asm/{sample}.p_ctg.gfa",
        "results/Assembly/Contigging/Solo_Asm/{sample}.a_ctg.gfa"
    threads:
        config["threads"]
    params:
        purgelevel= config['hifiasm']['purgelevel'],
        similarity= config['hifiasm']['similarity']
    log:
        "logs/{sample}.solo_assembly.log"
    benchmark:
        "benchmarks/{sample}.solo_assembly.txt"
    singularity:
        f"{config["sif_dir"]}/hifiasm.sif"
    shell:
        """
        hifiasm -t {threads} {params} \
        -o results/Assembly/Contigging/Solo_Asm/{wildcards.sample} --primary \
        {input.hifi} \
        >> {log} 2>&1
        """
