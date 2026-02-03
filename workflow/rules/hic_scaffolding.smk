rule yahs_scaffolding_hap1:
    input:
        fasta=lambda wildcards: "results/Decontamination/Contaminants/Phased_Asm_Hap1/{sample}.hic.hap1.p_ctg.clean.fasta" if config["gxdb"] else "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa",
        merged_bam="results/Scaffolding/Initial_Contacts/Hap1/{sample}.merged.bam",
    output:
        "results/Scaffolding/YAHS_Scaffolding/Hap1/{sample}.yahs_scaffolds_final.fa"
    params:
        yahs=config['yahs']
    log:
        "logs/{sample}.yahs_scaffolding_hap1.txt"
    benchmark:
        "benchmarks/{sample}.yahs_scaffolding_hap1.txt"
    singularity:
       f"{config["sif_dir"]}/yahs.sif"
    shell:
        """
        samtools faidx {input.fasta} >> {log} 2>&1 && \
        yahs {params.yahs} -o results/Scaffolding/YAHS_Scaffolding/Hap1/{wildcards.sample}.yahs {input.fasta} {input.merged_bam} >> {log} 2>&1
        """

rule yahs_scaffolding_hap2:
    input:
        fasta=lambda wildcards: "results/Decontamination/Contaminants/Phased_Asm_Hap2/{sample}.hic.hap2.p_ctg.clean.fasta" if config["gxdb"] else "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa",
        merged_bam="results/Scaffolding/Initial_Contacts/Hap2/{sample}.merged.bam",
    output:
        "results/Scaffolding/YAHS_Scaffolding/Hap2/{sample}.yahs_scaffolds_final.fa"
    params:
        yahs=config['yahs']
    log:
        "logs/{sample}.yahs_scaffolding_hap2.txt"
    benchmark:
        "benchmarks/{sample}.yahs_scaffolding_hap2.txt"
    singularity:
       f"{config["sif_dir"]}/yahs.sif"
    shell:
        """
        samtools faidx {input.fasta} >> {log} 2>&1 && \
        yahs {params.yahs} -o results/Scaffolding/YAHS_Scaffolding/Hap2/{wildcards.sample}.yahs {input.fasta} {input.merged_bam} >> {log} 2>&1
        """

rule yahs_scaffolding_prim:
    input:
        fasta=lambda wildcards: "results/Decontamination/Contaminants/Solo_Asm_Primary/{sample}.p_ctg.clean.fasta" if config["gxdb"] else "results/Assembly/Contigging/Solo_Asm/{sample}.p_ctg.fa",
        merged_bam="results/Scaffolding/Initial_Contacts/Primary/{sample}.merged.bam",
    output:
        "results/Scaffolding/YAHS_Scaffolding/Primary/{sample}.yahs_scaffolds_final.fa"
    params:
        yahs=config['yahs']
    log:
        "logs/{sample}.yahs_scaffolding_prim.txt"
    benchmark:
        "benchmarks/{sample}.yahs_scaffolding_prim.txt"
    singularity:
       f"{config["sif_dir"]}/yahs.sif"
    shell:
        """
        samtools faidx {input.fasta} >> {log} 2>&1 && \
        yahs {params.yahs} -o results/Scaffolding/YAHS_Scaffolding/Primary/{wildcards.sample}.yahs {input.fasta} {input.merged_bam} >> {log} 2>&1
        """

rule yahs_scaffolding_alt:
    input:
        fasta=lambda wildcards: "results/Decontamination/Contaminants/Solo_Asm_Alt/{sample}.a_ctg.clean.fasta" if config["gxdb"] else "results/Assembly/Contigging/Solo_Asm/{sample}.a_ctg.fa",
        merged_bam="results/Scaffolding/Initial_Contacts/Alternative/{sample}.merged.bam",
    output:
        "results/Scaffolding/YAHS_Scaffolding/Alternative/{sample}.yahs_scaffolds_final.fa"
    params:
        yahs=config['yahs']
    log:
        "logs/{sample}.yahs_scaffolding_alt.txt"
    benchmark:
        "benchmarks/{sample}.yahs_scaffolding_alt.txt"
    singularity:
       f"{config["sif_dir"]}/yahs.sif"
    shell:
        """
        samtools faidx {input.fasta} >> {log} 2>&1 && \
        yahs {params.yahs} -o results/Scaffolding/YAHS_Scaffolding/Alternative/{wildcards.sample}.yahs {input.fasta} {input.merged_bam} >> {log} 2>&1
        """
