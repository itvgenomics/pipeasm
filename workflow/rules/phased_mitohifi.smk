rule phased_run_mitohifi:
    input:
        hifi_hap1="results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa",
        hifi_hap2="results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa",
        reference_fasta="resources/{sample}.reference.fasta",
        reference_gb="resources/{sample}.reference.gb",
        sif="resources/mitohifi.sif"
    output:
        "results/Assembly/Mitogenome/Phased_Asm/{sample}.contigs_stats.tsv"
    threads:
        config["threads"]
    log:
        "logs/{sample}.run_mitohifi.log"
    benchmark:
        "benchmarks/{sample}.run_mitohifi.txt"
    params:
        config["geneticcode"]
    singularity:
        "resources/mitohifi.sif"
    shell:
        """
        cat {input.hifi_hap1} {input.hifi_hap2} >> results/Assembly/Mitogenome/Phased_Asm/concat.fasta && \
        cd results/Assembly/Mitogenome/Phased_Asm && \
        mitohifi.py -t {threads} -c concat.fasta -f ../../../../{input.reference_fasta} -g ../../../../{input.reference_gb} -o {params} >> ../../../../{log} 2>&1 && \
        mv contigs_stats.tsv {wildcards.sample}.contigs_stats.tsv && rm concat.fasta
        """