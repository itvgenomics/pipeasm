rule solo_edit_busco_table:
    input:
        "results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Primary/{sample}.summary.txt",
        "results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Alt/{sample}.summary.txt"
    output:
        "results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Primary/{sample}.full_table_busco_format_edit.tsv",
        "results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Alt/{sample}.full_table_busco_format_edit.tsv"
    log:
        "logs/{sample}.solo_edit_busco_table.log"
    threads: 1
    script:
        "../scripts/solo_edit_busco_table.py"

rule solo_run_snailplot_create_primary:
    input:
        fasta="results/Assembly/Contigging/Solo_Asm/{sample}.p_ctg.fa",
        busco="results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Primary/{sample}.full_table_busco_format_edit.tsv"
    output:
        "results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/00-Solo-Hap1/{sample}.bloobtools.create.check"
    log:
        "logs/{sample}.solo_run_snailplot_create_primary.log"
    benchmark:
        "benchmarks/{sample}.solo_run_snailplot_create_primary.txt"
    threads: 1
    singularity:
        "docker://genomehubs/blobtoolkit:4.3.5"
    shell:
        """
        blobtools create --replace \
        --fasta {input.fasta} \
        --busco {input.busco} \
        results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/00-Solo-Hap1/{wildcards.sample} >> {log} 2>&1 && \
        mv results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/00-Solo-Hap1/{wildcards.sample}/* results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/00-Solo-Hap1/ && \
        rm -r results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/00-Solo-Hap1/{wildcards.sample}/ && \
        touch results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/00-Solo-Hap1/{wildcards.sample}.bloobtools.create.check
        """

rule solo_run_snailplot_plot_primary:
    input:
        fasta="results/Assembly/Contigging/Solo_Asm/{sample}.p_ctg.fa",
        busco="results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Primary/{sample}.full_table_busco_format_edit.tsv",
        create_check="results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/00-Solo-Hap1/{sample}.bloobtools.create.check"
    output:
        "results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/00-Solo-Hap1/{sample}_Solo_Hap1.snail.png"
    log:
        "logs/{sample}.solo_run_snailplot_plot_primary.log"
    benchmark:
        "benchmarks/{sample}.solo_run_snailplot_plot_primary.txt"
    threads: 1
    singularity:
        "docker://genomehubs/blobtoolkit:4.3.5"
    shell:
        """
        blobtools view --plot --view snail results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/00-Solo-Hap1/ >> {log} 2>&1 && \
        mv 00-Solo-Hap1.snail.png results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/00-Solo-Hap1/{wildcards.sample}_Solo_Hap1.snail.png
        """

rule solo_run_snailplot_create_alt:
    input:
        fasta="results/Assembly/Contigging/Solo_Asm/{sample}.a_ctg.fa",
        busco="results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Alt/{sample}.full_table_busco_format_edit.tsv"
    output:
        "results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/01-Solo-Hap2/{sample}.bloobtools.create.check"
    log:
        "logs/{sample}.solo_run_snailplot_create_alt.log"
    benchmark:
        "benchmarks/{sample}.solo_run_snailplot_create_alt.txt"
    threads: 1
    singularity:
        "docker://genomehubs/blobtoolkit:4.3.5"
    shell:
        """
        blobtools create --replace \
        --fasta {input.fasta} \
        --busco {input.busco} \
        results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/01-Solo-Hap2/{wildcards.sample} >> {log} 2>&1 && \
        mv results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/01-Solo-Hap2/{wildcards.sample}/* results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/01-Solo-Hap2/ && \
        rm -r results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/01-Solo-Hap2/{wildcards.sample}/ && \
        touch results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/01-Solo-Hap2/{wildcards.sample}.bloobtools.create.check
        """

rule solo_run_snailplot_plot_alt:
    input:
        fasta="results/Assembly/Contigging/Solo_Asm/{sample}.a_ctg.fa",
        busco="results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Alt/{sample}.full_table_busco_format_edit.tsv",
        create_check="results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/01-Solo-Hap2/{sample}.bloobtools.create.check"
    output:
        "results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/01-Solo-Hap2/{sample}_Solo_Hap2.snail.png"
    log:
        "logs/{sample}.solo_run_snailplot_plot_alt.log"
    benchmark:
        "benchmarks/{sample}.solo_run_snailplot_plot_alt.txt"
    threads: 1
    singularity:
        "docker://genomehubs/blobtoolkit:4.3.5"
    shell:
        """
        blobtools view --plot --view snail results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/01-Solo-Hap2/ >> {log} 2>&1 && \
        mv 01-Solo-Hap2.snail.png results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/01-Solo-Hap2/{wildcards.sample}_Solo_Hap2.snail.png
        """