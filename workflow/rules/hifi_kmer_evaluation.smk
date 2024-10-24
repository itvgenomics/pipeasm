rule meryl_count:
    input:
        "results/Trimming_QC/HiFi/{sample}.trimmed.fastq.gz"
    output:
        directory("results/Trimming_QC/Meryl_DB/{sample}.trimmed.meryl")
    threads:
        config["threads"]
    params:
        kmer= config['meryl']['kmer']
    log:
        "logs/{sample}.meryl_count.log"
    benchmark:
        "benchmarks/{sample}.meryl_count.txt"
    singularity:
        "docker://itvdsbioinfo/meryl:1.3"
    shell:
        """
        meryl threads={threads} {params} count \
        results/Trimming_QC/HiFi/{wildcards.sample}.trimmed.fastq.gz \
        output results/Trimming_QC/Meryl_DB/{wildcards.sample}.trimmed.meryl \
        >> {log} 2>&1
        """

rule meryl_union:
    input:
        directory("results/Trimming_QC/Meryl_DB/{sample}.trimmed.meryl")
    output:
        directory("results/Trimming_QC/Meryl_DB/{sample}.trimmed.union.meryl")
    threads:
        config["threads"]
    params:
        kmer= config['meryl']['kmer']
    log:
        "logs/{sample}.meryl_union.log"
    benchmark:
        "benchmarks/{sample}.meryl_union.txt"
    singularity:
        "docker://itvdsbioinfo/meryl:1.3"
    shell:
        """
        meryl threads={threads} {params} union-sum \
        output results/Trimming_QC/Meryl_DB/{wildcards.sample}.trimmed.union.meryl \
        results/Trimming_QC/Meryl_DB/{wildcards.sample}.trimmed.meryl \
        >> {log} 2>&1
        """

rule meryl_histo:
    input:
        directory("results/Trimming_QC/Meryl_DB/{sample}.trimmed.union.meryl")
    output:
       "results/Trimming_QC/Meryl_DB/{sample}.trimmed.meryl.histo"
    threads:
        config["threads"]
    params:
        kmer= config['meryl']['kmer']
    log:
        "logs/{sample}.meryl_histo.log"
    benchmark:
        "benchmarks/{sample}.meryl_histo.txt"
    singularity:
        "docker://itvdsbioinfo/meryl:1.3"
    shell:
        """
        bash -c 'meryl threads={threads} {params} histogram \
        results/Trimming_QC/Meryl_DB/{wildcards.sample}.trimmed.union.meryl > \
        results/Trimming_QC/Meryl_DB/{wildcards.sample}.trimmed.meryl.histo' 2> {log}
        """

rule genomescope2_pacbio:
    input:
        "results/Trimming_QC/Meryl_DB/{sample}.trimmed.meryl.histo"
    output:
        "results/Assembly/Genome_Stats/HiFi_GenomeScope2/{sample}_linear_plot.png",
        "results/Assembly/Genome_Stats/HiFi_GenomeScope2/{sample}_log_plot.png",
        "results/Assembly/Genome_Stats/HiFi_GenomeScope2/{sample}_transformed_linear_plot.png",
        "results/Assembly/Genome_Stats/HiFi_GenomeScope2/{sample}_transformed_log_plot.png",
        "results/Assembly/Genome_Stats/HiFi_GenomeScope2/{sample}_model.txt"
    threads:
        config["software_threads"]["genomescope2"]
    params:
        kmer=config["genomescope2"]["kmer"],
        ploidy=config["genomescope2"]["ploidy"]
    log:
        "logs/{sample}.genomescope2_pacbio.log"
    benchmark:
        "benchmarks/{sample}.genomescope2_pacbio.txt"
    singularity:
        "docker://abner12/genomescope:2.0"
    shell:
        """
        genomescope.R {params.ploidy} -i {input} \
        -o results/Assembly/Genome_Stats/HiFi_GenomeScope2/ \
        {params.kmer} -n {wildcards.sample} --testing >> {log} 2>&1
        """

rule smudgeplot:
    input:
        reads="results/Trimming_QC/HiFi/{sample}.trimmed.fastq.gz",
        genomescope="results/Assembly/Genome_Stats/HiFi_GenomeScope2/{sample}_model.txt"
    output:
        "results/Assembly/Genome_Stats/Smudgeplot/{sample}_smudgeplot_log10.pdf",
        "results/Assembly/Genome_Stats/Smudgeplot/{sample}_smudgeplot.pdf"
    threads:
        config["threads"]
    params:
        fastk= config['smudgeplot']['fastk'],
        ploidyplot= config['smudgeplot']['ploidyplot']
    log:
        "logs/{sample}.smudgeplot_count.log"
    benchmark:
        "benchmarks/{sample}.smudgeplot_count.txt"
    singularity:
        "docker://itvdsbioinfo/smudgeplot:0.3.0"
    shell:
        """
        mkdir -p results/Assembly/Genome_Stats/Smudgeplot && \
        FastK {params.fastk} -M64 -T{threads} {input.reads} -Nresults/Assembly/Genome_Stats/Smudgeplot/FastK_Table -Presults/Assembly/Genome_Stats/Smudgeplot >> {log} 2>&1 && \
        PloidyPlot {params.ploidyplot} -T{threads} -oresults/Assembly/Genome_Stats/Smudgeplot/kmerpairs -Presults/Assembly/Genome_Stats/Smudgeplot results/Assembly/Genome_Stats/Smudgeplot/FastK_Table >> {log} 2>&1 && \
        smudgeplot.py plot -t {wildcards.sample} -o results/Assembly/Genome_Stats/Smudgeplot/{wildcards.sample} results/Assembly/Genome_Stats/Smudgeplot/kmerpairs_text.smu >> {log} 2>&1
        """

rule katgcp:
    input:
        reads="results/Trimming_QC/HiFi/{sample}.trimmed.fastq.gz"
    output:
        "results/Assembly/Genome_Stats/KAT/{sample}.mx.png"
    threads:
        config["threads"]
    params:
        kmer= config['katgcp']['kmer'],
    log:
        "logs/{sample}.katgcp.log"
    benchmark:
        "benchmarks/{sample}.katgcp.txt"
    singularity:
        "docker://itvdsbioinfo/kat:2.4.0"
    shell:
        """
        kat gcp -t {threads} {params.kmer} -o results/Assembly/Genome_Stats/KAT/{wildcards.sample} {input.reads} >> {log} 2>&1
        """