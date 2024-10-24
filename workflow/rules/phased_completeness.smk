rule phased_compleasm_hap1:
    input:
        fasta="results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa",
        buscodbpath=directory("resources/{sample}_buscodb")
    output:
        "results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap1/{sample}.summary.txt"
    threads:
        config["threads"]
    params:
        buscodb= config['buscodb']
    log:
        "logs/{sample}.phased_compleasm_hap1.log"
    benchmark:
        "benchmarks/{sample}.phased_compleasm_hap1.txt"
    singularity:
        "docker://huangnengcsu/compleasm:v0.2.5"
    shell:
        """
        (bash -c 'compleasm run -a {input.fasta} \
        -o results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap1 \
        -l {params} -t {threads} -L {input.buscodbpath} >> {log} 2>&1 && \
        mv results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap1/summary.txt {output}')
        """

rule phased_compleasm_hap2:
    input:
        fasta="results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa",
        buscodbpath=directory("resources/{sample}_buscodb")
    output:
        "results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap2/{sample}.summary.txt"
    threads:
        config["threads"]
    params:
        buscodb= config['buscodb']
    log:
        "logs/{sample}.phased_compleasm_hap2.log"
    benchmark:
        "benchmarks/{sample}.phased_compleasm_hap2.txt"
    singularity:
        "docker://huangnengcsu/compleasm:v0.2.5"
    shell:
        """
        (bash -c 'compleasm run -a {input.fasta} \
        -o results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap2 \
        -l {params} -t {threads} -L {input.buscodbpath} >> {log} 2>&1 && \
        mv results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap2/summary.txt {output}')
        """