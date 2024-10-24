rule solo_compleasm_hap1:
    input:
        fasta="results/Assembly/Contigging/Solo_Asm/{sample}.p_ctg.fa",
        buscodbpath=directory("resources/{sample}_buscodb")
    output:
        "results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Primary/{sample}.summary.txt"
    threads:
        config["threads"]
    params:
        buscodb= config['buscodb']
    log:
        "logs/{sample}.solo_compleasm_hap1.log"
    benchmark:
        "benchmarks/{sample}.solo_compleasm_hap1.txt"
    singularity:
        "docker://huangnengcsu/compleasm:v0.2.5"
    shell:
        """
        (bash -c 'compleasm run -a {input.fasta} \
        -o results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Primary \
        -l {params} -t {threads} -L {input.buscodbpath} >> {log} 2>&1 && \
        mv results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Primary/summary.txt {output}')
        """

rule solo_compleasm_hap2:
    input:
        fasta="results/Assembly/Contigging/Solo_Asm/{sample}.a_ctg.fa",
        buscodbpath=directory("resources/{sample}_buscodb")
    output:
        "results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Alt/{sample}.summary.txt"
    threads:
        config["threads"]
    params:
        buscodb= config['buscodb']
    log:
        "logs/{sample}.solo_compleasm_hap2.log"
    benchmark:
        "benchmarks/{sample}.solo_compleasm_hap2.txt"
    singularity:
        "docker://huangnengcsu/compleasm:v0.2.5"
    shell:
        """
        (bash -c 'compleasm run -a {input.fasta} \
        -o results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Alt \
        -l {params} -t {threads} -L {input.buscodbpath} >> {log} 2>&1 && \
        mv results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Alt/summary.txt {output}')
        """