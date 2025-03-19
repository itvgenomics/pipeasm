if config["gxdb"]:
    index_out_hap1 = "results/Decontamination/Contaminants/Phased_Asm_Hap1/{sample}.hic.hap1.p_ctg.clean.fasta.bwt.2bit.64"
    index_out_hap2 = "results/Decontamination/Contaminants/Phased_Asm_Hap2/{sample}.hic.hap2.p_ctg.clean.fasta.bwt.2bit.64"
else:
    index_out_hap1 = "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa.bwt.2bit.64"
    index_out_hap2 = "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa.bwt.2bit.64"

rule initial_contacts_bwa_mem2_index_hap1:
    input:
        lambda wildcards: "results/Decontamination/Contaminants/Phased_Asm_Hap1/{sample}.hic.hap1.p_ctg.clean.fasta" if config["gxdb"] else "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa"
    output: index_out_hap1
    threads:
        config["software_threads"]["bwa_index"]
    log:
        "logs/{sample}.initial_contacts_bwa_mem2_index_hap1.txt"
    benchmark:
        "benchmarks/{sample}.initial_contacts_bwa_mem2_index_hap1.txt"
    singularity:
       f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        mkdir -p results/Scaffolding/Initial_Contacts/Hap1 && \
        bwa-mem2.avx index {input} >> {log} 2>&1
        """

rule initial_contacts_bwa_mem2_index_hap2:
    input:
        lambda wildcards: "results/Decontamination/Contaminants/Phased_Asm_Hap2/{sample}.hic.hap2.p_ctg.clean.fasta" if config["gxdb"] else "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa"
    output: index_out_hap2
    threads:
        config["software_threads"]["bwa_index"]
    log:
        "logs/{sample}.initial_contacts_bwa_mem2_index_hap2.txt"
    benchmark:
        "benchmarks/{sample}.initial_contacts_bwa_mem2_index_hap2.txt"
    singularity:
       f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        mkdir -p results/Scaffolding/Initial_Contacts/Hap2 && \
        bwa-mem2.avx index {input} >> {log} 2>&1
        """

rule initial_contacts_bwa_mem2_mapping_hap1_r1:
    input:
        fasta=lambda wildcards: "results/Decontamination/Contaminants/Phased_Asm_Hap1/{sample}.hic.hap1.p_ctg.clean.fasta" if config["gxdb"] else "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa",
        hic_r1="results/Trimming_QC/HiC/{sample}_R1.trimmed_paired.fastq.gz",
        index=lambda wildcards: "results/Decontamination/Contaminants/Phased_Asm_Hap1/{sample}.hic.hap1.p_ctg.clean.fasta.bwt.2bit.64" if config["gxdb"] else "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa.bwt.2bit.64"
    output:
        "results/Scaffolding/Initial_Contacts/Hap1/{sample}.R1.bam"
    threads:
        config["threads"]
    log:
        "logs/{sample}.initial_contacts_bwa_mem2_mapping_hap1_r1.txt"
    benchmark:
        "benchmarks/{sample}.initial_contacts_bwa_mem2_mapping_hap1_r1.txt"
    params:
        tmpdir="results/Scaffolding/Initial_Contacts/Hap1/tmp1"
    singularity:
       f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        mkdir -p {params.tmpdir} && \
        bwa-mem2.avx mem -t {threads} -v 1 {input.fasta} {input.hic_r1} | samtools sort -n -@ {threads} -T {params.tmpdir} -O bam -o {output} >> {log} 2>&1
        """

rule initial_contacts_bwa_mem2_mapping_hap1_r2:
    input:
        fasta=lambda wildcards: "results/Decontamination/Contaminants/Phased_Asm_Hap1/{sample}.hic.hap1.p_ctg.clean.fasta" if config["gxdb"] else "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa",
        hic_r2="results/Trimming_QC/HiC/{sample}_R2.trimmed_paired.fastq.gz",
        index=lambda wildcards: "results/Decontamination/Contaminants/Phased_Asm_Hap1/{sample}.hic.hap1.p_ctg.clean.fasta.bwt.2bit.64" if config["gxdb"] else "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa.bwt.2bit.64"
    output:
        "results/Scaffolding/Initial_Contacts/Hap1/{sample}.R2.bam"
    threads:
        config["threads"]
    log:
        "logs/{sample}.initial_contacts_bwa_mem2_mapping_hap1_r2.txt"
    benchmark:
        "benchmarks/{sample}.initial_contacts_bwa_mem2_mapping_hap1_r2.txt"
    params:
        tmpdir="results/Scaffolding/Initial_Contacts/Hap1/tmp2"
    singularity:
       f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        mkdir -p {params.tmpdir} && \
        bwa-mem2.avx mem -t {threads} -v 1 {input.fasta} {input.hic_r2} | samtools sort -n -@ {threads} -T {params.tmpdir} -O bam -o {output} >> {log} 2>&1
        """

rule initial_contacts_bwa_mem2_mapping_hap2_r1:
    input:
        fasta=lambda wildcards: "results/Decontamination/Contaminants/Phased_Asm_Hap2/{sample}.hic.hap2.p_ctg.clean.fasta" if config["gxdb"] else "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa",
        hic_r1="results/Trimming_QC/HiC/{sample}_R1.trimmed_paired.fastq.gz",
        index=lambda wildcards: "results/Decontamination/Contaminants/Phased_Asm_Hap2/{sample}.hic.hap2.p_ctg.clean.fasta.bwt.2bit.64" if config["gxdb"] else "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa.bwt.2bit.64"
    output:
        "results/Scaffolding/Initial_Contacts/Hap2/{sample}.R1.bam"
    threads:
        config["threads"]
    log:
        "logs/{sample}.initial_contacts_bwa_mem2_mapping_hap2_r1.txt"
    benchmark:
        "benchmarks/{sample}.initial_contacts_bwa_mem2_mapping_hap2_r1.txt"
    params:
        tmpdir="results/Scaffolding/Initial_Contacts/Hap2/tmp1"
    singularity:
       f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        mkdir -p {params.tmpdir} && \
        bwa-mem2.avx mem -t {threads} -v 1 {input.fasta} {input.hic_r1} | samtools sort -n -@ {threads} -T {params.tmpdir} -O bam -o {output} >> {log} 2>&1
        """

rule initial_contacts_bwa_mem2_mapping_hap2_r2:
    input:
        fasta=lambda wildcards: "results/Decontamination/Contaminants/Phased_Asm_Hap2/{sample}.hic.hap2.p_ctg.clean.fasta" if config["gxdb"] else "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa",
        hic_r2="results/Trimming_QC/HiC/{sample}_R2.trimmed_paired.fastq.gz",
        index=lambda wildcards: "results/Decontamination/Contaminants/Phased_Asm_Hap2/{sample}.hic.hap2.p_ctg.clean.fasta.bwt.2bit.64" if config["gxdb"] else "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa.bwt.2bit.64"
    output:
        "results/Scaffolding/Initial_Contacts/Hap2/{sample}.R2.bam"
    threads:
        config["threads"]
    log:
        "logs/{sample}.initial_contacts_bwa_mem2_mapping_hap2_r2.txt"
    benchmark:
        "benchmarks/{sample}.initial_contacts_bwa_mem2_mapping_hap2_r2.txt"
    params:
        tmpdir="results/Scaffolding/Initial_Contacts/Hap2/tmp2"
    singularity:
       f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        mkdir -p {params.tmpdir} && \
        bwa-mem2.avx mem -t {threads} -v 1 {input.fasta} {input.hic_r2} | samtools sort -n -@ {threads} -T {params.tmpdir} -O bam -o {output} >> {log} 2>&1
        """

rule initial_contacts_bellerophon_merge_hap1:
    input:
        bam_r1="results/Scaffolding/Initial_Contacts/Hap1/{sample}.R1.bam",
        bam_r2="results/Scaffolding/Initial_Contacts/Hap1/{sample}.R2.bam",
    output:
        "results/Scaffolding/Initial_Contacts/Hap1/{sample}.merged.bam"
    threads:
        config["threads"]
    params:
        quality= config['bellerophon'],
    log:
        "logs/{sample}.initial_contacts_bellerophon_merge_hap1.txt"
    benchmark:
        "benchmarks/{sample}.initial_contacts_bellerophon_merge_hap1.txt"
    singularity:
       f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        bellerophon --forward {input.bam_r1} --reverse {input.bam_r2} --output {output} {params} --threads {threads} >> {log} 2>&1
        """

rule initial_contacts_bellerophon_merge_hap2:
    input:
        bam_r1="results/Scaffolding/Initial_Contacts/Hap2/{sample}.R1.bam",
        bam_r2="results/Scaffolding/Initial_Contacts/Hap2/{sample}.R2.bam",
    output:
        "results/Scaffolding/Initial_Contacts/Hap2/{sample}.merged.bam"
    threads:
        config["threads"]
    params:
        quality= config['bellerophon'],
    log:
        "logs/{sample}.initial_contacts_bellerophon_merge_hap2.txt"
    benchmark:
        "benchmarks/{sample}.initial_contacts_bellerophon_merge_hap2.txt"
    singularity:
       f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        bellerophon --forward {input.bam_r1} --reverse {input.bam_r2} --output {output} {params} --threads {threads} >> {log} 2>&1
        """

rule initial_contacts_samtools_flagstats_hap1:
    input:
        merged_bam="results/Scaffolding/Initial_Contacts/Hap1/{sample}.merged.bam",
    output:
        flagstats="results/Scaffolding/Initial_Contacts/Hap1/{sample}.merged.bam.flagstats",
        hic_stats="results/Scaffolding/Initial_Contacts/Hap1/{sample}.merged.bam.hic_stats"
    threads:
        config["software_threads"]["samtools_flagstats"]
    log:
        "logs/{sample}.initial_contacts_samtools_flagstats_hap1.txt"
    benchmark:
        "benchmarks/{sample}.initial_contacts_samtools_flagstats_hap1.txt"
    singularity:
       f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        samtools flagstat -@ {threads} -O tsv {input.merged_bam} > {output.flagstats} 2>{log} && \
        perl workflow/scripts/get_stats.pl {input.merged_bam} > {output.hic_stats} 2>{log}
        """

rule initial_contacts_samtools_flagstats_hap2:
    input:
        merged_bam="results/Scaffolding/Initial_Contacts/Hap2/{sample}.merged.bam",
    output:
        flagstats="results/Scaffolding/Initial_Contacts/Hap2/{sample}.merged.bam.flagstats",
        hic_stats="results/Scaffolding/Initial_Contacts/Hap2/{sample}.merged.bam.hic_stats"
    threads:
        config["software_threads"]["samtools_flagstats"]
    log:
        "logs/{sample}.initial_contacts_samtools_flagstats_hap2.txt"
    benchmark:
        "benchmarks/{sample}.initial_contacts_samtools_flagstats_hap2.txt"
    singularity:
       f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        samtools flagstat -@ {threads} -O tsv {input.merged_bam} > {output.flagstats} 2>{log} && \
        perl workflow/scripts/get_stats.pl {input.merged_bam} > {output.hic_stats} 2>{log}
        """

rule initial_contacts_pretext_map_snapshot_hap1:
    input:
        merged_bam="results/Scaffolding/Initial_Contacts/Hap1/{sample}.merged.bam",
    output:
        "results/Scaffolding/Initial_Contacts/Hap1/{sample}.merged.bam.pretext"
    threads:
        config["software_threads"]["pretext_snapshot"]
    params:
        pretextmap= config['pretext']["map"],
        snapshot=  config['pretext']["snapshot"]
    log:
        "logs/{sample}.initial_contacts_pretext_map_snapshot_hap1.txt"
    benchmark:
        "benchmarks/{sample}.initial_contacts_pretext_map_snapshot_hap1.txt"
    singularity:
       f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        samtools view -h {input.merged_bam} | PretextMap -o {output} {params.pretextmap} >> {log} 2>&1 && \
        PretextSnapshot -m {output} {params.snapshot} -o results/Scaffolding/Initial_Contacts/Hap1/ --prefix {wildcards.sample}_Hap1 >> {log} 2>&1
        """

rule initial_contacts_pretext_map_snapshot_hap2:
    input:
        merged_bam="results/Scaffolding/Initial_Contacts/Hap2/{sample}.merged.bam",
    output:
        "results/Scaffolding/Initial_Contacts/Hap2/{sample}.merged.bam.pretext"
    threads:
        config["software_threads"]["pretext_snapshot"]
    params:
        pretextmap= config['pretext']["map"],
        snapshot=  config['pretext']["snapshot"]
    log:
        "logs/{sample}.initial_contacts_pretext_map_snapshot_hap2.txt"
    benchmark:
        "benchmarks/{sample}.initial_contacts_pretext_map_snapshot_hap2.txt"
    singularity:
       f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        samtools view -h {input.merged_bam} | PretextMap -o {output} {params.pretextmap} >> {log} 2>&1 && \
        PretextSnapshot -m {output} {params.snapshot} -o results/Scaffolding/Initial_Contacts/Hap2/ --prefix {wildcards.sample}_Hap2 >> {log} 2>&1
        """

rule initial_contacts_samtools_flagstats_bam_r1_hap1:
    input:
        bam_r1="results/Scaffolding/Initial_Contacts/Hap1/{sample}.R1.bam",
    output:
        flagstats_r1="results/Scaffolding/Initial_Contacts/Hap1/{sample}.R1.bam.flagstats"
    threads:
        config["software_threads"]["samtools_flagstats"]
    params:
        pretextmap= config['pretext']["map"],
        snapshot=  config['pretext']["snapshot"]
    log:
        "logs/{sample}.initial_contacts_samtools_flagstats_bam_r1_hap1.txt"
    benchmark:
        "benchmarks/{sample}.initial_contacts_samtools_flagstats_bam_r1_hap1.txt"
    singularity:
       f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        samtools flagstat -@ {threads} -O tsv {input.bam_r1} > {output.flagstats_r1} 2>{log}
        """

rule initial_contacts_samtools_flagstats_bam_r2_hap1:
    input:
        bam_r2="results/Scaffolding/Initial_Contacts/Hap1/{sample}.R2.bam",
    output:
        flagstats_r2="results/Scaffolding/Initial_Contacts/Hap1/{sample}.R2.bam.flagstats"
    threads:
        config["software_threads"]["samtools_flagstats"]
    params:
        pretextmap= config['pretext']["map"],
        snapshot=  config['pretext']["snapshot"]
    log:
        "logs/{sample}.initial_contacts_samtools_flagstats_bam_r2_hap1.txt"
    benchmark:
        "benchmarks/{sample}.initial_contacts_samtools_flagstats_bam_r2_hap1.txt"
    singularity:
       f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        samtools flagstat -@ {threads} -O tsv {input.bam_r2} > {output.flagstats_r2} 2>{log}
        """

rule initial_contacts_samtools_flagstats_bam_r1_hap2:
    input:
        bam_r1="results/Scaffolding/Initial_Contacts/Hap2/{sample}.R1.bam",
    output:
        flagstats_r1="results/Scaffolding/Initial_Contacts/Hap2/{sample}.R1.bam.flagstats"
    threads:
        config["software_threads"]["samtools_flagstats"]
    params:
        pretextmap= config['pretext']["map"],
        snapshot=  config['pretext']["snapshot"]
    log:
        "logs/{sample}.initial_contacts_samtools_flagstats_bam_r1_hap2.txt"
    benchmark:
        "benchmarks/{sample}.initial_contacts_samtools_flagstats_bam_r1_hap2.txt"
    singularity:
       f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        samtools flagstat -@ {threads} -O tsv {input.bam_r1} > {output.flagstats_r1} 2>{log}
        """

rule initial_contacts_samtools_flagstats_bam_r2_hap2:
    input:
        bam_r2="results/Scaffolding/Initial_Contacts/Hap2/{sample}.R2.bam",
    output:
        flagstats_r2="results/Scaffolding/Initial_Contacts/Hap2/{sample}.R2.bam.flagstats"
    threads:
        config["software_threads"]["samtools_flagstats"]
    params:
        pretextmap= config['pretext']["map"],
        snapshot=  config['pretext']["snapshot"]
    log:
        "logs/{sample}.initial_contacts_samtools_flagstats_bam_r2_hap2.txt"
    benchmark:
        "benchmarks/{sample}.initial_contacts_samtools_flagstats_bam_r2_hap2.txt"
    singularity:
       f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        samtools flagstat -@ {threads} -O tsv {input.bam_r2} > {output.flagstats_r2} 2>{log}
        """
