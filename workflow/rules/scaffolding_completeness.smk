rule scaffolding_compleasm_hap1:
    input:
        fasta="results/Scaffolding/YAHS_Scaffolding/Hap1/{sample}.yahs_scaffolds_final.fa",
        buscodbpath=directory("resources/{sample}_buscodb")
    output:
        "results/Scaffolding/Scaffolding_stats/Compleasm/Hap1/{sample}.summary.txt"
    params:
        buscodb= config['buscodb']
    log:
        "logs/{sample}.scaffolding_compleasm_hap1.log"
    benchmark:
        "benchmarks/{sample}.scaffolding_compleasm_hap1.txt"
    singularity:
        f"{config["sif_dir"]}/compleasm.sif"
    shell:
        """
        (bash -c 'compleasm run -a {input.fasta} \
        -o results/Scaffolding/Scaffolding_stats/Compleasm/Hap1 \
        -l {params} -t {threads} -L {input.buscodbpath} >> {log} 2>&1 && \
        mv results/Scaffolding/Scaffolding_stats/Compleasm/Hap1/summary.txt {output}')
        """

rule scaffolding_compleasm_hap2:
    input:
        fasta="results/Scaffolding/YAHS_Scaffolding/Hap2/{sample}.yahs_scaffolds_final.fa",
        buscodbpath=directory("resources/{sample}_buscodb")
    output:
        "results/Scaffolding/Scaffolding_stats/Compleasm/Hap2/{sample}.summary.txt"
    params:
        buscodb= config['buscodb']
    log:
        "logs/{sample}.scaffolding_compleasm_hap2.log"
    benchmark:
        "benchmarks/{sample}.scaffolding_compleasm_hap2.txt"
    singularity:
        f"{config["sif_dir"]}/compleasm.sif"
    shell:
        """
        (bash -c 'compleasm run -a {input.fasta} \
        -o results/Scaffolding/Scaffolding_stats/Compleasm/Hap2 \
        -l {params} -t {threads} -L {input.buscodbpath} >> {log} 2>&1 && \
        mv results/Scaffolding/Scaffolding_stats/Compleasm/Hap2/summary.txt {output}')
        """

rule scaffolding_compleasm_prim:
    input:
        fasta="results/Scaffolding/YAHS_Scaffolding/Primary/{sample}.yahs_scaffolds_final.fa",
        buscodbpath=directory("resources/{sample}_buscodb")
    output:
        "results/Scaffolding/Scaffolding_stats/Compleasm/Primary/{sample}.summary.txt"
    params:
        buscodb= config['buscodb']
    log:
        "logs/{sample}.scaffolding_compleasm_prim.log"
    benchmark:
        "benchmarks/{sample}.scaffolding_compleasm_prim.txt"
    singularity:
        f"{config["sif_dir"]}/compleasm.sif"
    shell:
        """
        (bash -c 'compleasm run -a {input.fasta} \
        -o results/Scaffolding/Scaffolding_stats/Compleasm/Primary/ \
        -l {params} -t {threads} -L {input.buscodbpath} >> {log} 2>&1 && \
        mv results/Scaffolding/Scaffolding_stats/Compleasm/Primary/summary.txt {output}')
        """

rule scaffolding_compleasm_alt:
    input:
        fasta="results/Scaffolding/YAHS_Scaffolding/Alternative/{sample}.yahs_scaffolds_final.fa",
        buscodbpath=directory("resources/{sample}_buscodb")
    output:
        "results/Scaffolding/Scaffolding_stats/Compleasm/Alternative/{sample}.summary.txt"
    params:
        buscodb= config['buscodb']
    log:
        "logs/{sample}.scaffolding_compleasm_alt.log"
    benchmark:
        "benchmarks/{sample}.scaffolding_compleasm_alt.txt"
    singularity:
        f"{config["sif_dir"]}/compleasm.sif"
    shell:
        """
        (bash -c 'compleasm run -a {input.fasta} \
        -o results/Scaffolding/Scaffolding_stats/Compleasm/Alternative/ \
        -l {params} -t {threads} -L {input.buscodbpath} >> {log} 2>&1 && \
        mv results/Scaffolding/Scaffolding_stats/Compleasm/Alternative/summary.txt {output}')
        """
