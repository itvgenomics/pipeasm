rule phased_edit_busco_table:
    input:
        "results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap1/{sample}.summary.txt",
        "results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap2/{sample}.summary.txt"
    output:
        "results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap1/{sample}.full_table_busco_format_edit.tsv",
        "results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap2/{sample}.full_table_busco_format_edit.tsv"
    log:
        "logs/{sample}.phased_edit_busco_table.log"
    threads: 1
    script:
        "../scripts/phased_edit_busco_table.py"
        

rule phased_run_snailplot_create_hap1:
    input:
        fasta="results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa",
        busco="results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap1/{sample}.full_table_busco_format_edit.tsv"
    output:
        "results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap1/{sample}.bloobtools.create.check"
    log:
        "logs/{sample}.phased_run_snailplot_create_hap1.log"
    benchmark:
        "benchmarks/{sample}.phased_run_snailplot_create_hap1.txt"
    threads: 1
    singularity:
        "docker://genomehubs/blobtoolkit:4.3.5"
    shell:
        """
        blobtools create --replace \
        --fasta {input.fasta} \
        --busco {input.busco} \
        results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap1/{wildcards.sample} >> {log} 2>&1 && \
        mv results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap1/{wildcards.sample}/* results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap1/ && \
        rm -r results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap1/{wildcards.sample}/ && \
        touch results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap1/{wildcards.sample}.bloobtools.create.check
        """

rule phased_run_snailplot_plot_hap1:
    input:
        fasta="results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa",
        busco="results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap1/{sample}.full_table_busco_format_edit.tsv",
        create_check="results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap1/{sample}.bloobtools.create.check"
    output:
        "results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap1/{sample}_Phased_Hap1.snail.png"
    log:
        "logs/{sample}.phased_run_snailplot_plot_hap1.log"
    benchmark:
        "benchmarks/{sample}.phased_run_snailplot_plot_hap1.txt"
    threads: 1
    singularity:
        "docker://genomehubs/blobtoolkit:4.3.5"
    shell:
        """
        blobtools view --plot --view snail results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap1/ >> {log} 2>&1 && \
        mv Phased-Hap1.snail.png results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap1/{wildcards.sample}_Phased_Hap1.snail.png
        """

rule phased_run_snailplot_create_hap2:
    input:
        fasta="results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa",
        busco="results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap2/{sample}.full_table_busco_format_edit.tsv"
    output:
        "results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap2/{sample}.bloobtools.create.check"
    log:
        "logs/{sample}.phased_run_snailplot_create_hap2.log"
    benchmark:
        "benchmarks/{sample}.phased_run_snailplot_create_hap2.txt"
    threads: 1
    singularity:
        "docker://genomehubs/blobtoolkit:4.3.5"
    shell:
        """
        blobtools create --replace \
        --fasta {input.fasta} \
        --busco {input.busco} \
        results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap2/{wildcards.sample} >> {log} 2>&1 && \
        mv results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap2/{wildcards.sample}/* results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap2/ && \
        rm -r results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap2/{wildcards.sample}/ && \
        touch results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap2/{wildcards.sample}.bloobtools.create.check
        """

rule phased_run_snailplot_plot_hap2:
    input:
        fasta="results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa",
        busco="results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap2/{sample}.full_table_busco_format_edit.tsv",
        create_check="results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap2/{sample}.bloobtools.create.check"
    output:
        "results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap2/{sample}_Phased_Hap2.snail.png"
    log:
        "logs/{sample}.phased_run_snailplot_plot_hap2.log"
    benchmark:
        "benchmarks/{sample}.phased_run_snailplot_plot_hap2.txt"
    threads: 1
    singularity:
        "docker://genomehubs/blobtoolkit:4.3.5"
    shell:
        """
        blobtools view --plot --view snail results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap2/ >> {log} 2>&1 && \
        mv Phased-Hap2.snail.png results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap2/{wildcards.sample}_Phased_Hap2.snail.png
        """