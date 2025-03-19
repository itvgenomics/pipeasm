rule scaffolding_gfastats_stats_hap1:
    input:
        "results/Scaffolding/YAHS_Scaffolding/Hap1/{sample}.yahs_scaffolds_final.fa"
    output:
        "results/Scaffolding/Scaffolding_stats/GFAstats/{sample}.yahs_scaffolds_hap1.fa.stats"
    threads:
        config["software_threads"]["gfastats"]
    params:
        params= config['gfastats']['params']
    log:
        "logs/{sample}.scaffolding_gfastats_stats_hap1.log"
    benchmark:
        "benchmarks/{sample}.scaffolding_gfastats_stats_hap1.txt"
    singularity:
        f"{config["sif_dir"]}/gfastats.sif"
    shell:
        """
        bash -c 'gfastats -f {input} \
        -j {threads} > {output} 2> {log}'
        """

rule scaffolding_gfastats_stats_hap2:
    input:
        "results/Scaffolding/YAHS_Scaffolding/Hap2/{sample}.yahs_scaffolds_final.fa"
    output:
        "results/Scaffolding/Scaffolding_stats/GFAstats/{sample}.yahs_scaffolds_hap2.fa.stats"
    threads:
        config["software_threads"]["gfastats"]
    params:
        params= config['gfastats']['params']
    log:
        "logs/{sample}.scaffolding_gfastats_stats_hap2.log"
    benchmark:
        "benchmarks/{sample}.scaffolding_gfastats_stats_hap2.txt"
    singularity:
        f"{config["sif_dir"]}/gfastats.sif"
    shell:
        """
        bash -c 'gfastats -f {input} \
        -j {threads} > {output} 2> {log}'
        """
