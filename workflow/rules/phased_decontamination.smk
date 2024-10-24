rule phased_fcsadaptor_hap1:
    input:
        "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa"
    output:
        "results/Decontamination/FCS-Adaptor/Phased_Asm_Hap1/cleaned_sequences/{sample}.hic.hap1.p_ctg.fa"
    threads:
        config['threads']
    log:
        "logs/{sample}.phased_fcsadaptor_hap1.log"
    benchmark:
        "benchmarks/{sample}.phased_fcsadaptor_hap1.txt"
    shell:
        """
        bash -c 'export OMP_NUM_THREADS={threads} && \
        bash workflow/scripts/run_fcsadaptor.sh --fasta-input {input} --output-dir results/Decontamination/FCS-Adaptor/Phased_Asm_Hap1 --euk'
        """

rule phased_fcsadaptor_hap2:
    input:
        "results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa"
    output:
        "results/Decontamination/FCS-Adaptor/Phased_Asm_Hap2/cleaned_sequences/{sample}.hic.hap2.p_ctg.fa"
    threads:
        config['threads']
    log:
        "logs/{sample}.phased_fcsadaptor_hap2.log"
    benchmark:
        "benchmarks/{sample}.phased_fcsadaptor_hap2.txt"
    shell:
        """
        bash -c 'export OMP_NUM_THREADS={threads} && \
        bash workflow/scripts/run_fcsadaptor.sh --fasta-input {input} --output-dir results/Decontamination/FCS-Adaptor/Phased_Asm_Hap2 --euk'
        """

rule phased_fcsgx:
    input:
        hap1="results/Decontamination/FCS-Adaptor/Phased_Asm_Hap1/cleaned_sequences/{sample}.hic.hap1.p_ctg.fa",
        hap2="results/Decontamination/FCS-Adaptor/Phased_Asm_Hap2/cleaned_sequences/{sample}.hic.hap2.p_ctg.fa",
        sif="resources/fcs-gx.sif"
    output:
        hap1="results/Decontamination/Contaminants/Phased_Asm_Hap1/{sample}.hic.hap1.p_ctg.screen.check",
        hap2="results/Decontamination/Contaminants/Phased_Asm_Hap2/{sample}.hic.hap2.p_ctg.screen.check"
    threads:
        config['threads']
    params:
        taxid=config['taxid'],
        gxdb=config['gxdb']
    log:
        "logs/{sample}.phased_fcsgx.log"
    benchmark:
        "benchmarks/{sample}.phased_fcsgx.txt"
    shell:
        """
        bash -c 'echo "OMP_NUM_THREADS={threads}" > results/Decontamination/Contaminants/Phased_Asm_Hap1/gx_env.txt && \
        python workflow/scripts/fcs.py --env-file results/Decontamination/Contaminants/Phased_Asm_Hap1/gx_env.txt --image resources/fcs-gx.sif \
        screen genome \
        --fasta {input.hap1} \
        --out-dir results/Decontamination/Contaminants/Phased_Asm_Hap1 \
        --gx-db {params.gxdb} --tax-id {params.taxid} --generate-logfile True && \
        touch {output.hap1}'

        bash -c 'echo "OMP_NUM_THREADS={threads}" >>  results/Decontamination/Contaminants/Phased_Asm_Hap2/gx_env.txt && \
        python workflow/scripts/fcs.py --env-file results/Decontamination/Contaminants/Phased_Asm_Hap2/gx_env.txt --image resources/fcs-gx.sif \
        screen genome \
        --fasta {input.hap2} \
        --out-dir results/Decontamination/Contaminants/Phased_Asm_Hap2 \
        --gx-db {params.gxdb} --tax-id {params.taxid} --generate-logfile True && \
        touch {output.hap2}'
        """

rule phased_fcsgx_clean_hap1:
    input:
        hap1_fa="results/Decontamination/FCS-Adaptor/Phased_Asm_Hap1/cleaned_sequences/{sample}.hic.hap1.p_ctg.fa",
        hap1="results/Decontamination/Contaminants/Phased_Asm_Hap1/{sample}.hic.hap1.p_ctg.screen.check"
    output:
        hap1="results/Decontamination/Contaminants/Phased_Asm_Hap1/{sample}.hic.hap1.p_ctg.clean.fasta"
    threads:
        config['threads']
    params:
        config['taxid']
    log:
        "logs/{sample}.phased_fcsgx_clean_hap1.log"
    benchmark:
        "benchmarks/{sample}.phased_fcsgx_clean_hap1.txt"
    shell:
        """
        bash -c 'cat {input.hap1_fa} | python workflow/scripts/fcs.py \
        --image resources/fcs-gx.sif clean genome \
        --action-report results/Decontamination/Contaminants/Phased_Asm_Hap1/{wildcards.sample}.hic.hap1.p_ctg.{params}.fcs_gx_report.txt \
        --output {output.hap1} --contam-fasta-out results/Decontamination/Contaminants/Phased_Asm_Hap1/{wildcards.sample}.contam.fasta'
        """


rule phased_fcsgx_clean_hap2:
    input:
        hap2_fa="results/Decontamination/FCS-Adaptor/Phased_Asm_Hap2/cleaned_sequences/{sample}.hic.hap2.p_ctg.fa",
        hap2="results/Decontamination/Contaminants/Phased_Asm_Hap2/{sample}.hic.hap2.p_ctg.screen.check"
    output:
        hap2="results/Decontamination/Contaminants/Phased_Asm_Hap2/{sample}.hic.hap2.p_ctg.clean.fasta"
    threads:
        config['threads']
    params:
        config['taxid']
    log:
        "logs/{sample}.phased_fcsgx_clean_hap2.log"
    benchmark:
        "benchmarks/{sample}.phased_fcsgx_clean_hap2.txt"
    shell:
        """
        bash -c 'cat {input.hap2_fa} | python workflow/scripts/fcs.py \
        --image resources/fcs-gx.sif clean genome \
        --action-report results/Decontamination/Contaminants/Phased_Asm_Hap2/{wildcards.sample}.hic.hap2.p_ctg.{params}.fcs_gx_report.txt \
        --output {output.hap2} --contam-fasta-out results/Decontamination/Contaminants/Phased_Asm_Hap2/{wildcards.sample}.contam.fasta'
        """