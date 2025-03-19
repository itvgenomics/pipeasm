rule phased_assembly:
    input:
        hifi="results/Trimming_QC/HiFi/{sample}.trimmed.fastq.gz",
        hic_r1="results/Trimming_QC/HiC/{sample}_R1.trimmed_paired.fastq.gz",
        hic_r2="results/Trimming_QC/HiC/{sample}_R2.trimmed_paired.fastq.gz",
    output:
        expand("results/Assembly/Contigging/Phased_Asm/{{sample}}.hic.hap{hap}.p_ctg.gfa", hap=["1", "2"])
    threads:
        config["threads"]
    params:
        purgelevel= config['hifiasm']['purgelevel'],
        similarity= config['hifiasm']['similarity']
    log:
        "logs/{sample}.phased_assembly.log"
    benchmark:
        "benchmarks/{sample}.phased_assembly.txt"
    singularity:
        f"{config["sif_dir"]}/hifiasm.sif"
    shell:
        """
        hifiasm -t {threads} {params} \
        -o results/Assembly/Contigging/Phased_Asm/{wildcards.sample} \
        --h1 {input.hic_r1} \
        --h2 {input.hic_r2} \
        {input.hifi} \
        >> {log} 2>&1
        """
