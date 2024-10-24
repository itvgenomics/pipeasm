rule final_contacts_bwa_mem2_index_hap1:
    input:
        fasta="results/Scaffolding/YAHS_Scaffolding/Hap1/{sample}.yahs_scaffolds_final.fa",
    output:
        "results/Scaffolding/YAHS_Scaffolding/Hap1/{sample}.yahs_scaffolds_final.fa.bwt.2bit.64"
    threads:
        config["software_threads"]["bwa_index"]
    log:
        "logs/{sample}.final_contacts_bwa_mem2_index_hap1.txt"
    benchmark:
        "benchmarks/{sample}.final_contacts_bwa_mem2_index_hap1.txt"
    singularity:
       "docker://itvdsbioinfo/hic_mapping:1.0"
    shell:
        """
        bwa-mem2.avx index {input.fasta} >> {log} 2>&1
        """

rule final_contacts_bwa_mem2_index_hap2:
    input:
        fasta="results/Scaffolding/YAHS_Scaffolding/Hap2/{sample}.yahs_scaffolds_final.fa",
    output:
        "results/Scaffolding/YAHS_Scaffolding/Hap2/{sample}.yahs_scaffolds_final.fa.bwt.2bit.64"
    threads:
        config["software_threads"]["bwa_index"]
    log:
        "logs/{sample}.final_contacts_bwa_mem2_index_hap2.txt"
    benchmark:
        "benchmarks/{sample}.final_contacts_bwa_mem2_index_hap2.txt"
    singularity:
       "docker://itvdsbioinfo/hic_mapping:1.0"
    shell:
        """
        bwa-mem2.avx index {input.fasta} >> {log} 2>&1
        """

rule final_contacts_bwa_mem2_mapping_hap1_r1:
    input:
        fasta="results/Scaffolding/YAHS_Scaffolding/Hap1/{sample}.yahs_scaffolds_final.fa",
        hic_r1="results/Trimming_QC/HiC/{sample}_R1.trimmed_paired.fastq.gz",
        index="results/Scaffolding/YAHS_Scaffolding/Hap1/{sample}.yahs_scaffolds_final.fa.bwt.2bit.64",
    output:
        "results/Scaffolding/Final_Contacts/Hap1/{sample}.R1.bam"
    threads:
        config["threads"]
    log:
        "logs/{sample}.final_contacts_bwa_mem2_mapping_hap1_r1.txt"
    benchmark:
        "benchmarks/{sample}.final_contacts_bwa_mem2_mapping_hap1_r1.txt"
    params:
        tmpdir="results/Scaffolding/Final_Contacts/Hap1/tmp1"
    singularity:
       "docker://itvdsbioinfo/hic_mapping:1.0"
    shell:
        """
        mkdir -p {params.tmpdir} && \
        bwa-mem2.avx mem -t {threads} -v 1 {input.fasta} {input.hic_r1} | samtools sort -n -@ {threads} -T {params.tmpdir} -O bam -o {output} >> {log} 2>&1
        """

rule final_contacts_bwa_mem2_mapping_hap1_r2:
    input:
        fasta="results/Scaffolding/YAHS_Scaffolding/Hap1/{sample}.yahs_scaffolds_final.fa",
        hic_r2="results/Trimming_QC/HiC/{sample}_R2.trimmed_paired.fastq.gz",
        index="results/Scaffolding/YAHS_Scaffolding/Hap1/{sample}.yahs_scaffolds_final.fa.bwt.2bit.64",
    output:
        "results/Scaffolding/Final_Contacts/Hap1/{sample}.R2.bam"
    threads:
        config["threads"]
    log:
        "logs/{sample}.final_contacts_bwa_mem2_mapping_hap1_r2.txt"
    benchmark:
        "benchmarks/{sample}.final_contacts_bwa_mem2_mapping_hap1_r2.txt"
    params:
        tmpdir="results/Scaffolding/Final_Contacts/Hap1/tmp2"
    singularity:
       "docker://itvdsbioinfo/hic_mapping:1.0"
    shell:
        """
        mkdir -p {params.tmpdir} && \
        bwa-mem2.avx mem -t {threads} -v 1 {input.fasta} {input.hic_r2} | samtools sort -n -@ {threads} -T {params.tmpdir} -O bam -o {output} >> {log} 2>&1
        """

rule final_contacts_bwa_mem2_mapping_hap2_r1:
    input:
        fasta="results/Scaffolding/YAHS_Scaffolding/Hap2/{sample}.yahs_scaffolds_final.fa",
        hic_r1="results/Trimming_QC/HiC/{sample}_R1.trimmed_paired.fastq.gz",
        index="results/Scaffolding/YAHS_Scaffolding/Hap2/{sample}.yahs_scaffolds_final.fa.bwt.2bit.64",
    output:
        "results/Scaffolding/Final_Contacts/Hap2/{sample}.R1.bam"
    threads:
        config["threads"]
    log:
        "logs/{sample}.final_contacts_bwa_mem2_mapping_hap2_r1.txt"
    benchmark:
        "benchmarks/{sample}.final_contacts_bwa_mem2_mapping_hap2_r1.txt"
    params:
        tmpdir="results/Scaffolding/Final_Contacts/Hap2/tmp1"
    singularity:
       "docker://itvdsbioinfo/hic_mapping:1.0"
    shell:
        """
        mkdir -p {params.tmpdir} && \
        bwa-mem2.avx mem -t {threads} -v 1 {input.fasta} {input.hic_r1} | samtools sort -n -@ {threads} -T {params.tmpdir} -O bam -o results/Scaffolding/Final_Contacts/Hap2/{wildcards.sample}.R1.bam >> {log} 2>&1
        """

rule final_contacts_bwa_mem2_mapping_hap2_r2:
    input:
        fasta="results/Scaffolding/YAHS_Scaffolding/Hap2/{sample}.yahs_scaffolds_final.fa",
        hic_r2="results/Trimming_QC/HiC/{sample}_R2.trimmed_paired.fastq.gz",
        index="results/Scaffolding/YAHS_Scaffolding/Hap2/{sample}.yahs_scaffolds_final.fa.bwt.2bit.64",
    output:
        "results/Scaffolding/Final_Contacts/Hap2/{sample}.R2.bam"
    threads:
        config["threads"]
    log:
        "logs/{sample}.final_contacts_bwa_mem2_mapping_hap2_r2.txt"
    benchmark:
        "benchmarks/{sample}.final_contacts_bwa_mem2_mapping_hap2_r2.txt"
    params:
        tmpdir="results/Scaffolding/Final_Contacts/Hap2/tmp2"
    singularity:
       "docker://itvdsbioinfo/hic_mapping:1.0"
    shell:
        """
        mkdir -p {params.tmpdir} && \
        bwa-mem2.avx mem -t {threads} -v 1 {input.fasta} {input.hic_r2} | samtools sort -n -@ {threads} -T {params.tmpdir} -O bam -o {output} >> {log} 2>&1
        """

rule final_contacts_bellerophon_merge_hap1:
    input:
        bam_r1="results/Scaffolding/Final_Contacts/Hap1/{sample}.R1.bam",
        bam_r2="results/Scaffolding/Final_Contacts/Hap1/{sample}.R2.bam",
    output:
        "results/Scaffolding/Final_Contacts/Hap1/{sample}.merged.bam"
    threads:
        config["threads"]
    params:
        quality= config['bellerophon'],
    log:
        "logs/{sample}.final_contacts_bellerophon_merge_hap1.txt"
    benchmark:
        "benchmarks/{sample}.final_contacts_bellerophon_merge_hap1.txt"
    singularity:
       "docker://itvdsbioinfo/hic_mapping:1.0"
    shell:
        """
        bellerophon --forward {input.bam_r1} --reverse {input.bam_r2} --output {output} {params} --threads {threads} >> {log} 2>&1
        """

rule final_contacts_bellerophon_merge_hap2:
    input:
        bam_r1="results/Scaffolding/Final_Contacts/Hap2/{sample}.R1.bam",
        bam_r2="results/Scaffolding/Final_Contacts/Hap2/{sample}.R2.bam",
    output:
        "results/Scaffolding/Final_Contacts/Hap2/{sample}.merged.bam"
    threads:
        config["threads"]
    params:
        quality= config['bellerophon'],
    log:
        "logs/{sample}.final_contacts_bellerophon_merge_hap2.txt"
    benchmark:
        "benchmarks/{sample}.final_contacts_bellerophon_merge_hap2.txt"
    singularity:
       "docker://itvdsbioinfo/hic_mapping:1.0"
    shell:
        """
        bellerophon --forward {input.bam_r1} --reverse {input.bam_r2} --output {output} {params} --threads {threads} >> {log} 2>&1
        """

rule final_contacts_samtools_flagstats_hap1:
    input:
        merged_bam="results/Scaffolding/Final_Contacts/Hap1/{sample}.merged.bam",
    output:
        flagstats="results/Scaffolding/Final_Contacts/Hap1/{sample}.merged.bam.flagstats",
        hic_stats="results/Scaffolding/Final_Contacts/Hap1/{sample}.merged.bam.hic_stats"
    threads:
        config["software_threads"]["samtools_flagstats"]
    log:
        "logs/{sample}.final_contacts_samtools_flagstats_hap1.txt"
    benchmark:
        "benchmarks/{sample}.final_contacts_samtools_flagstats_hap1.txt"
    singularity:
       "docker://itvdsbioinfo/hic_mapping:1.0"
    shell:
        """
        samtools flagstat -@ {threads} -O tsv {input.merged_bam} > {output.flagstats} 2>{log} && \
        perl workflow/scripts/get_stats.pl {input.merged_bam} > {output.hic_stats} 2>{log}
        """

rule final_contacts_samtools_flagstats_hap2:
    input:
        merged_bam="results/Scaffolding/Final_Contacts/Hap2/{sample}.merged.bam",
    output:
        flagstats="results/Scaffolding/Final_Contacts/Hap2/{sample}.merged.bam.flagstats",
        hic_stats="results/Scaffolding/Final_Contacts/Hap2/{sample}.merged.bam.hic_stats"
    threads:
        config["software_threads"]["samtools_flagstats"]
    log:
        "logs/{sample}.final_contacts_samtools_flagstats_hap2.txt"
    benchmark:
        "benchmarks/{sample}.final_contacts_samtools_flagstats_hap2.txt"
    singularity:
       "docker://itvdsbioinfo/hic_mapping:1.0"
    shell:
        """
        samtools flagstat -@ {threads} -O tsv {input.merged_bam} > {output.flagstats} 2>{log} && \
        perl workflow/scripts/get_stats.pl {input.merged_bam} > {output.hic_stats} 2>{log}
        """

rule final_contacts_pretext_map_snapshot_hap1:
    input:
        merged_bam="results/Scaffolding/Final_Contacts/Hap1/{sample}.merged.bam",
    output:
        "results/Scaffolding/Final_Contacts/Hap1/{sample}.merged.bam.pretext"
    threads:
        config["software_threads"]["pretext_snapshot"]
    params:
        pretextmap= config['pretext']["map"],
        snapshot=  config['pretext']["snapshot"]
    log:
        "logs/{sample}.final_contacts_pretext_map_snapshot_hap1.txt"
    benchmark:
        "benchmarks/{sample}.final_contacts_pretext_map_snapshot_hap1.txt"
    singularity:
       "docker://itvdsbioinfo/hic_mapping:1.0"
    shell:
        """
        samtools view -h {input.merged_bam} | PretextMap -o {output} {params.pretextmap} >> {log} 2>&1 && \
        PretextSnapshot -m {output} {params.snapshot} -o results/Scaffolding/Final_Contacts/Hap1/ --prefix {wildcards.sample}.final_contacts_Hap1 >> {log} 2>&1
        """

rule final_contacts_pretext_map_snapshot_hap2:
    input:
        merged_bam="results/Scaffolding/Final_Contacts/Hap2/{sample}.merged.bam",
    output:
        "results/Scaffolding/Final_Contacts/Hap2/{sample}.merged.bam.pretext"
    threads:
        config["software_threads"]["pretext_snapshot"]
    params:
        pretextmap= config['pretext']["map"],
        snapshot=  config['pretext']["snapshot"]
    log:
        "logs/{sample}.final_contacts_pretext_map_snapshot_hap2.txt"
    benchmark:
        "benchmarks/{sample}.final_contacts_pretext_map_snapshot_hap2.txt"
    singularity:
       "docker://itvdsbioinfo/hic_mapping:1.0"
    shell:
        """
        samtools view -h {input.merged_bam} | PretextMap -o {output} {params.pretextmap} >> {log} 2>&1 && \
        PretextSnapshot -m {output} {params.snapshot} -o results/Scaffolding/Final_Contacts/Hap2/ --prefix {wildcards.sample}.final_contacts_Hap2 >> {log} 2>&1
        """

rule final_contacts_samtools_flagstats_bam_r1_hap1:
    input:
        bam_r1="results/Scaffolding/Final_Contacts/Hap1/{sample}.R1.bam",
    output:
        flagstats_r1="results/Scaffolding/Final_Contacts/Hap1/{sample}.R1.bam.flagstats"
    threads:
        config["software_threads"]["samtools_flagstats"]
    params:
        pretextmap= config['pretext']["map"],
        snapshot=  config['pretext']["snapshot"]
    log:
        "logs/{sample}.final_contacts_samtools_flagstats_bam_r1_hap1.txt"
    benchmark:
        "benchmarks/{sample}.final_contacts_samtools_flagstats_bam_r1_hap1.txt"
    singularity:
       "docker://itvdsbioinfo/hic_mapping:1.0"
    shell:
        """
        samtools flagstat -@ {threads} -O tsv {input.bam_r1} > {output.flagstats_r1} 2>{log}
        """

rule final_contacts_samtools_flagstats_bam_r2_hap1:
    input:
        bam_r2="results/Scaffolding/Final_Contacts/Hap1/{sample}.R2.bam",
    output:
        flagstats_r2="results/Scaffolding/Final_Contacts/Hap1/{sample}.R2.bam.flagstats"
    threads:
        config["software_threads"]["samtools_flagstats"]
    params:
        pretextmap= config['pretext']["map"],
        snapshot=  config['pretext']["snapshot"]
    log:
        "logs/{sample}.final_contacts_samtools_flagstats_bam_r2_hap1.txt"
    benchmark:
        "benchmarks/{sample}.final_contacts_samtools_flagstats_bam_r2_hap1.txt"
    singularity:
       "docker://itvdsbioinfo/hic_mapping:1.0"
    shell:
        """
        samtools flagstat -@ {threads} -O tsv {input.bam_r2} > {output.flagstats_r2} 2>{log}
        """

rule final_contacts_samtools_flagstats_bam_r1_hap2:
    input:
        bam_r1="results/Scaffolding/Final_Contacts/Hap2/{sample}.R1.bam",
    output:
        flagstats_r1="results/Scaffolding/Final_Contacts/Hap2/{sample}.R1.bam.flagstats"
    threads:
        config["software_threads"]["samtools_flagstats"]
    params:
        pretextmap= config['pretext']["map"],
        snapshot=  config['pretext']["snapshot"]
    log:
        "logs/{sample}.final_contacts_samtools_flagstats_bam_r1_hap2.txt"
    benchmark:
        "benchmarks/{sample}.final_contacts_samtools_flagstats_bam_r1_hap2.txt"
    singularity:
       "docker://itvdsbioinfo/hic_mapping:1.0"
    shell:
        """
        samtools flagstat -@ {threads} -O tsv {input.bam_r1} > {output.flagstats_r1} 2>{log}
        """

rule final_contacts_samtools_flagstats_bam_r2_hap2:
    input:
        bam_r2="results/Scaffolding/Final_Contacts/Hap2/{sample}.R2.bam",
    output:
        flagstats_r2="results/Scaffolding/Final_Contacts/Hap2/{sample}.R2.bam.flagstats"
    threads:
        config["software_threads"]["samtools_flagstats"]
    params:
        pretextmap= config['pretext']["map"],
        snapshot=  config['pretext']["snapshot"]
    log:
        "logs/{sample}.final_contacts_samtools_flagstats_bam_r2_hap2.txt"
    benchmark:
        "benchmarks/{sample}.final_contacts_samtools_flagstats_bam_r2_hap2.txt"
    singularity:
       "docker://itvdsbioinfo/hic_mapping:1.0"
    shell:
        """
        samtools flagstat -@ {threads} -O tsv {input.bam_r2} > {output.flagstats_r2} 2>{log}
        """