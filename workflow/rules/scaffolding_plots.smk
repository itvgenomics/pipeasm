rule scaffolding_edit_busco_table:
    input:
        "results/Scaffolding/Scaffolding_stats/Compleasm/Hap1/{sample}.summary.txt",
        "results/Scaffolding/Scaffolding_stats/Compleasm/Hap2/{sample}.summary.txt"
    output:
        "results/Scaffolding/Scaffolding_stats/Compleasm/Hap1/{sample}.full_table_busco_format_edit.tsv",
        "results/Scaffolding/Scaffolding_stats/Compleasm/Hap2/{sample}.full_table_busco_format_edit.tsv"
    log:
        "logs/{sample}.scaffolding_edit_busco_table.log"
    threads: 1
    script:
        "../scripts/scaffolding_edit_busco_table.py"

rule scaffolding_run_snailplot_create_hap1:
    input:
        fasta="results/Scaffolding/YAHS_Scaffolding/Hap1/{sample}.yahs_scaffolds_final.fa",
        busco="results/Scaffolding/Scaffolding_stats/Compleasm/Hap1/{sample}.full_table_busco_format_edit.tsv"
    output:
        "results/Scaffolding/Scaffolding_stats/SnailPlot/Hap1/{sample}.bloobtools.create.check"
    log:
        "logs/{sample}.scaffolding_run_snailplot_create_hap1.log"
    benchmark:
        "benchmarks/{sample}.scaffolding_run_snailplot_create_hap1.txt"
    threads: 1
    singularity:
        "docker://genomehubs/blobtoolkit:4.3.5"
    shell:
        """
        blobtools create --replace \
        --fasta {input.fasta} \
        --busco {input.busco} \
        results/Scaffolding/Scaffolding_stats/SnailPlot/Hap1/{wildcards.sample} >> {log} 2>&1 && \
        mv results/Scaffolding/Scaffolding_stats/SnailPlot/Hap1/{wildcards.sample}/* results/Scaffolding/Scaffolding_stats/SnailPlot/Hap1/ && \
        rm -r results/Scaffolding/Scaffolding_stats/SnailPlot/Hap1/{wildcards.sample}/ && \
        touch results/Scaffolding/Scaffolding_stats/SnailPlot/Hap1/{wildcards.sample}.bloobtools.create.check
        """

rule scaffolding_run_snailplot_plot_hap1:
    input:
        fasta="results/Scaffolding/YAHS_Scaffolding/Hap1/{sample}.yahs_scaffolds_final.fa",
        busco="results/Scaffolding/Scaffolding_stats/Compleasm/Hap1/{sample}.full_table_busco_format_edit.tsv",
        create_check="results/Scaffolding/Scaffolding_stats/SnailPlot/Hap1/{sample}.bloobtools.create.check"
    output:
        "results/Scaffolding/Scaffolding_stats/SnailPlot/Hap1/{sample}_Scaffolding_Hap1.snail.png"
    log:
        "logs/{sample}.scaffolding_run_snailplot_plot_hap1.log"
    benchmark:
        "benchmarks/{sample}.scaffolding_run_snailplot_plot_hap1.txt"
    threads: 1
    singularity:
        "docker://genomehubs/blobtoolkit:4.3.5"
    shell:
        """
        blobtools view --plot --view snail results/Scaffolding/Scaffolding_stats/SnailPlot/Hap1/ >> {log} 2>&1 && \
        mv Hap1.snail.png results/Scaffolding/Scaffolding_stats/SnailPlot/Hap1/{wildcards.sample}_Scaffolding_Hap1.snail.png
        """

rule scaffolding_run_snailplot_create_hap2:
    input:
        fasta="results/Scaffolding/YAHS_Scaffolding/Hap2/{sample}.yahs_scaffolds_final.fa",
        busco="results/Scaffolding/Scaffolding_stats/Compleasm/Hap2/{sample}.full_table_busco_format_edit.tsv"
    output:
        "results/Scaffolding/Scaffolding_stats/SnailPlot/Hap2/{sample}.bloobtools.create.check"
    log:
        "logs/{sample}.scaffolding_run_snailplot_create_hap2.log"
    benchmark:
        "benchmarks/{sample}.scaffolding_run_snailplot_create_hap2.txt"
    threads: 1
    singularity:
        "docker://genomehubs/blobtoolkit:4.3.5"
    shell:
        """
        blobtools create --replace \
        --fasta {input.fasta} \
        --busco {input.busco} \
        results/Scaffolding/Scaffolding_stats/SnailPlot/Hap2/{wildcards.sample} >> {log} 2>&1 && \
        mv results/Scaffolding/Scaffolding_stats/SnailPlot/Hap2/{wildcards.sample}/* results/Scaffolding/Scaffolding_stats/SnailPlot/Hap2/ && \
        rm -r results/Scaffolding/Scaffolding_stats/SnailPlot/Hap2/{wildcards.sample}/ && \
        touch results/Scaffolding/Scaffolding_stats/SnailPlot/Hap2/{wildcards.sample}.bloobtools.create.check
        """

rule scaffolding_run_snailplot_plot_hap2:
    input:
        fasta="results/Scaffolding/YAHS_Scaffolding/Hap2/{sample}.yahs_scaffolds_final.fa",
        busco="results/Scaffolding/Scaffolding_stats/Compleasm/Hap2/{sample}.full_table_busco_format_edit.tsv",
        create_check="results/Scaffolding/Scaffolding_stats/SnailPlot/Hap2/{sample}.bloobtools.create.check"
    output:
        "results/Scaffolding/Scaffolding_stats/SnailPlot/Hap2/{sample}_Scaffolding_Hap2.snail.png"
    log:
        "logs/{sample}.scaffolding_run_snailplot_plot_hap2.log"
    benchmark:
        "benchmarks/{sample}.scaffolding_run_snailplot_plot_hap2.txt"
    threads: 1
    singularity:
        "docker://genomehubs/blobtoolkit:4.3.5"
    shell:
        """
        blobtools view --plot --view snail results/Scaffolding/Scaffolding_stats/SnailPlot/Hap2/ >> {log} 2>&1 && \
        mv Hap2.snail.png results/Scaffolding/Scaffolding_stats/SnailPlot/Hap2/{wildcards.sample}_Scaffolding_Hap2.snail.png
        """