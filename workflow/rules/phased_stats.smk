rule phased_gfastats_fasta_hap1:
    input:
        "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.gfa"
    output:
        "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa"
    threads:
        config["software_threads"]["gfastats"]
    params:
        params= config['gfastats']['params']
    log:
        "logs/{sample}.phased_gfastats_fasta_hap1.log"
    benchmark:
        "benchmarks/{sample}.phased_gfastats_fasta_hap1.txt"
    singularity:
        f"{config["sif_dir"]}/gfastats.sif"
    shell:
        """
        bash -c 'gfastats -f {input} \
        -j {threads} {params} > {output} \
        2> {log}'
        """

rule phased_gfastats_stats_hap1:
    input:
        "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa"
    output:
        "results/Assembly/Genome_Stats/GFAstats/{sample}.hic.hap1.p_ctg.fa.stats"
    threads:
        config["software_threads"]["gfastats"]
    params:
        params= config['gfastats']['params']
    log:
        "logs/{sample}.phased_gfastats_stats_hap1.log"
    benchmark:
        "benchmarks/{sample}.phased_gfastats_stats_hap1.txt"
    singularity:
        f"{config["sif_dir"]}/gfastats.sif"
    shell:
        """
        bash -c 'gfastats -f {input} \
        -j {threads} > {output} 2> {log}'
        """

rule phased_gfastats_fasta_hap2:
    input:
        "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.gfa"
    output:
        "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa"
    threads:
        config["software_threads"]["gfastats"]
    params:
        params= config['gfastats']['params']
    log:
        "logs/{sample}.phased_gfastats_fasta_hap2.log"
    benchmark:
        "benchmarks/{sample}.phased_gfastats_fasta_hap2.txt"
    singularity:
        f"{config["sif_dir"]}/gfastats.sif"
    shell:
        """
        bash -c 'gfastats -f {input} \
        -j {threads} {params} > {output} \
        2> {log}'
        """

rule phased_gfastats_stats_hap2:
    input:
        "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa"
    output:
        "results/Assembly/Genome_Stats/GFAstats/{sample}.hic.hap2.p_ctg.fa.stats"
    threads:
        config["software_threads"]["gfastats"]
    params:
        params= config['gfastats']['params']
    log:
        "logs/{sample}.phased_gfastats_stats_hap2.log"
    benchmark:
        "benchmarks/{sample}.phased_gfastats_stats_hap2.txt"
    singularity:
        f"{config["sif_dir"]}/gfastats.sif"
    shell:
        """
        bash -c 'gfastats -f {input} \
        -j {threads} > {output} 2> {log}'
        """
