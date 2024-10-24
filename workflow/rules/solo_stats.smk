rule solo_gfastats_fasta_primary:
    input:
        "results/Assembly/Contigging/Solo_Asm/{sample}.p_ctg.gfa"
    output:
        "results/Assembly/Contigging/Solo_Asm/{sample}.p_ctg.fa"
    threads:
        config["software_threads"]["gfastats"]
    params:
        params= config['gfastats']['params']
    log:
        "logs/{sample}.solo_gfastats_fasta_primary.log"
    benchmark:
        "benchmarks/{sample}.solo_gfastats_fasta_primary.txt"
    singularity:
        "docker://staphb/gfastats:1.3.6"
    shell:
        """
        bash -c 'gfastats -f {input} \
        -j {threads} {params} > {output} \
        2> {log}'
        """

rule solo_gfastats_stats_primary:
    input:
        "results/Assembly/Contigging/Solo_Asm/{sample}.p_ctg.fa"
    output:
        "results/Assembly/Genome_Stats/GFAstats/{sample}.p_ctg.fa.stats"
    threads:
        config["software_threads"]["gfastats"]
    params:
        params= config['gfastats']['params']
    log:
        "logs/{sample}.solo_gfastats_stats_primary.log"
    benchmark:
        "benchmarks/{sample}.solo_gfastats_stats_primary.txt"
    singularity:
        "docker://staphb/gfastats:1.3.6"
    shell:
        """
        bash -c 'gfastats -f {input} \
        -j {threads} > {output} 2> {log}'
        """

rule solo_gfastats_fasta_alt:
    input:
        "results/Assembly/Contigging/Solo_Asm/{sample}.a_ctg.gfa"
    output:
        "results/Assembly/Contigging/Solo_Asm/{sample}.a_ctg.fa"
    threads:
        config["software_threads"]["gfastats"]
    params:
        params= config['gfastats']['params']
    log:
        "logs/{sample}.solo_gfastats_fasta_alt.log"
    benchmark:
        "benchmarks/{sample}.solo_gfastats_fasta_alt.txt"
    singularity:
        "docker://staphb/gfastats:1.3.6"
    shell:
        """
        bash -c 'gfastats -f {input} \
        -j {threads} {params} > {output} \
        2> {log}'
        """

rule solo_gfastats_stats_alt:
    input:
        "results/Assembly/Contigging/Solo_Asm/{sample}.a_ctg.fa"
    output:
        "results/Assembly/Genome_Stats/GFAstats/{sample}.a_ctg.fa.stats"
    threads:
        config["software_threads"]["gfastats"]
    params:
        params= config['gfastats']['params']
    log:
        "logs/{sample}.solo_gfastats_stats_alt.log"
    benchmark:
        "benchmarks/{sample}.solo_gfastats_stats_alt.txt"
    singularity:
        "docker://staphb/gfastats:1.3.6"
    shell:
        """
        bash -c 'gfastats -f {input} \
        -j {threads} > {output} 2> {log}'
        """