rule fastk:
    input:
        "results/Trimming_QC/HiFi/{sample}.trimmed.fastq.gz"
    output:
        "results/Trimming_QC/FastK/{sample}.hist",
        "results/Trimming_QC/FastK/{sample}.ktab",
    params:
        kmer= config['fastk']['kmer'],
        default= config['fastk']['default']
    log:
        "logs/{sample}.fastk.log"
    benchmark:
        "benchmarks/{sample}.fastk.txt"
    singularity:
        f"{config["sif_dir"]}/fkutils.sif"
    shell:
        """
        mkdir -p results/Trimming_QC/FastK && \
        FastK {params.kmer} {params.default} -T{threads} -Ptmp {input} -Nresults/Trimming_QC/FastK/{wildcards.sample} >> {log} 2>&1
        """

rule genescopeFK:
    input:
        "results/Trimming_QC/FastK/{sample}.hist"
    output:
        "results/Assembly/Genome_Stats/GeneScopeFK/{sample}_linear_plot.png",
        "results/Assembly/Genome_Stats/GeneScopeFK/{sample}_log_plot.png",
        "results/Assembly/Genome_Stats/GeneScopeFK/{sample}_transformed_linear_plot.png",
        "results/Assembly/Genome_Stats/GeneScopeFK/{sample}_transformed_log_plot.png",
        "results/Assembly/Genome_Stats/GeneScopeFK/{sample}_model.txt"
    params:
        kmer=config["genescopeFK"]["kmer"],
        ploidy=config["genescopeFK"]["ploidy"],
        outdir="results/Assembly/Genome_Stats/GeneScopeFK"
    log:
        "logs/{sample}.genescopeFK.log"
    benchmark:
        "benchmarks/{sample}.genescopeFK.txt"
    singularity:
        f"{config["sif_dir"]}/fkutils.sif"
    shell:
        """
        mkdir -p {params.outdir} && \
        Histex -G {input} | GeneScopeFK.R -o {params.outdir} {params.kmer} {params.ploidy} 2>{log} && \
        mv {params.outdir}/linear_plot.png {params.outdir}/{wildcards.sample}_linear_plot.png && \
        mv {params.outdir}/log_plot.png {params.outdir}/{wildcards.sample}_log_plot.png && \
        mv {params.outdir}/transformed_linear_plot.png {params.outdir}/{wildcards.sample}_transformed_linear_plot.png && \
        mv {params.outdir}/transformed_log_plot.png {params.outdir}/{wildcards.sample}_transformed_log_plot.png && \
        mv {params.outdir}/model.txt {params.outdir}/{wildcards.sample}_model.txt && \
        mv {params.outdir}/summary.txt {params.outdir}/{wildcards.sample}_summary.txt
        """

rule smudgeplot:
    input:
        reads="results/Trimming_QC/HiFi/{sample}.trimmed.fastq.gz",
        fastk="results/Trimming_QC/FastK/{sample}.ktab"
    output:
        "results/Assembly/Genome_Stats/Smudgeplot/{sample}_smudgeplot_log10.png",
        "results/Assembly/Genome_Stats/Smudgeplot/{sample}_smudgeplot.png"
    params:
        minsize= config['smudgeplot']['minsize'],
        outdir="results/Assembly/Genome_Stats/Smudgeplot"
    log:
        "logs/{sample}.smudgeplot_count.log"
    benchmark:
        "benchmarks/{sample}.smudgeplot_count.txt"
    singularity:
        f"{config["sif_dir"]}/fkutils.sif"
    shell:
        """
        mkdir -p {params.outdir} && \
        smudgeplot hetmers {params.minsize} -t {threads} -tmp tmp -o {params.outdir}/{wildcards.sample} --verbose {input.fastk} >> {log} 2>&1 && \
        smudgeplot all -o {params.outdir}/{wildcards.sample} {params.outdir}/{wildcards.sample}.smu >> {log} 2>&1
        """

rule katgc:
    input:
        "results/Trimming_QC/FastK/{sample}.ktab"
    output:
        "results/Assembly/Genome_Stats/KatGC/{sample}.st.png"
    log:
        "logs/{sample}.katgc.log"
    benchmark:
        "benchmarks/{sample}.katgc.txt"
    singularity:
        f"{config["sif_dir"]}/fkutils.sif"
    shell:
        """
        mkdir -p results/Assembly/Genome_Stats/KatGC/ && \
        cd results/Trimming_QC/FastK/ && \
        KatGC -lfs -T{threads} {wildcards.sample}.ktab ../../../results/Assembly/Genome_Stats/KatGC/{wildcards.sample} >> ../../../{log} 2>&1
        """

rule katcomp:
    input:
        hap1="results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa",
        hap2="results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa"
    output:
        "results/Assembly/Genome_Stats/KatComp/{sample}.st.png"
    params:
        kmer= config['fastk']['kmer'],
        default= config['fastk']['default'],
        outdir="results/Assembly/Genome_Stats/KatComp"
    log:
        "logs/{sample}.katcomp.log"
    benchmark:
        "benchmarks/{sample}.katcomp.txt"
    singularity:
        f"{config["sif_dir"]}/fkutils.sif"
    shell:
        """
        mkdir -p fkutils/KatComp fkutils/KatGC && \
        FastK {params.kmer} -T{threads} -Ptmp {params.default} {input.hap1} -Nfkutils/KatComp/{wildcards.sample}.Hap1 && \
        FastK {params.kmer} -T{threads} -Ptmp {params.default} {input.hap2} -Nfkutils/KatComp/{wildcards.sample}.Hap2 && \
        cd {params.outdir} && \
        KatComp -lfs -T16 {wildcards.sample}.Hap1.ktab {wildcards.sample}.Hap2.ktab {wildcards.sample} >> {log} 2>&1
        """
