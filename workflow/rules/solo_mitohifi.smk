rule solo_run_mitohifi:
    input:
        hifi_hap1="results/Assembly/Contigging/Solo_Asm/{sample}.p_ctg.fa",
        hifi_hap2="results/Assembly/Contigging/Solo_Asm/{sample}.a_ctg.fa",
        reference_fasta="resources/{sample}.reference.fasta",
        reference_gb="resources/{sample}.reference.gb"
    output:
        "results/Assembly/Mitogenome/Solo_Asm/{sample}.contigs_stats.tsv"
    threads:
        config["threads"]
    log:
        "logs/{sample}.run_mitohifi.log"
    benchmark:
        "benchmarks/{sample}.run_mitohifi.txt"
    params:
        config["geneticcode"]
    singularity:
        f"{config["sif_dir"]}/mitohifi.sif"
    shell:
        """
        cat {input.hifi_hap1} {input.hifi_hap2} >> results/Assembly/Mitogenome/Solo_Asm/concat.fasta && \
        cd results/Assembly/Mitogenome/Solo_Asm && \
        mitohifi.py -t {threads} -c concat.fasta -f ../../../../{input.reference_fasta} -g ../../../../{input.reference_gb} -o {params} >> ../../../../{log} 2>&1 && \
        mv contigs_stats.tsv {wildcards.sample}.contigs_stats.tsv && rm concat.fasta
        """
