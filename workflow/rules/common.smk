import sys

## Function to get all expected output files
def get_output(sampĺe_name):

    out = []

    ## Check if the input files exists
    out.extend(expand("results/Trimming_QC/HiFi/{sample}.fastq.gz", sample=sampĺe_name))

    if config['hic_r1'] and config['hic_r2']:
        out.extend(expand("results/Trimming_QC/HiC/{sample}_R{pair}.fastq.gz", sample=sampĺe_name, pair=["1","2"]))

    if config['ont_reads']:
        out.extend(expand("results/Trimming_QC/ONT/{sample}_ONT.fastq.gz", sample=sampĺe_name))

    ## Check the Trimming and QC output files
    if config['run_trimming_qc'].lower() == "yes":
        
        out.extend(expand("results/Trimming_QC/HiFi/{sample}.trimmed.fastq.gz", sample=sampĺe_name))

        out.extend(expand("results/Trimming_QC/QC/HiFi_FastQC/{sample}.trimmed_fastqc.{ext}", sample=config["sample"], ext=['zip','html']))

        nanoplot_out = [expand("results/Trimming_QC/QC/HiFi_NanoPlot/{sample}LengthvsQualityScatterPlot_dot.svg", sample=config["sample"]),
                        expand("results/Trimming_QC/QC/HiFi_NanoPlot/{sample}LengthvsQualityScatterPlot_kde.svg", sample=config["sample"]),
                        expand("results/Trimming_QC/QC/HiFi_NanoPlot/{sample}Non_weightedHistogramReadlength.svg", sample=config["sample"]),
                        expand("results/Trimming_QC/QC/HiFi_NanoPlot/{sample}Non_weightedLogTransformed_HistogramReadlength.svg", sample=config["sample"]),
                        expand("results/Trimming_QC/QC/HiFi_NanoPlot/{sample}WeightedHistogramReadlength.svg", sample=config["sample"]),
                        expand("results/Trimming_QC/QC/HiFi_NanoPlot/{sample}WeightedLogTransformed_HistogramReadlength.svg", sample=config["sample"]),
                        expand("results/Trimming_QC/QC/HiFi_NanoPlot/{sample}Yield_By_Length.svg", sample=config["sample"])]

        for file in nanoplot_out:
            out.append(file)

        if config["hic_r1"] and config["hic_r2"]:
            out.extend(expand("results/Trimming_QC/HiC/{sample}_R1.trimmed_paired.fastq.gz", sample=sampĺe_name))
            out.extend(expand("results/Trimming_QC/HiC/{sample}_R2.trimmed_paired.fastq.gz", sample=sampĺe_name))
            out.extend(expand("results/Trimming_QC/QC/HiC_FastQC/{sample}_R{pair}.trimmed_paired_fastqc.{ext}", sample=config["sample"], pair=['1','2'], ext=['zip','html']))
        
        if config['ont_reads']:
            
            out.extend(expand("results/Trimming_QC/ONT/{sample}_ONT.trimmed.fastq.gz", sample=sampĺe_name))

            ont_nanoplot_out = [expand("results/Trimming_QC/QC/ONT_NanoPlot/{sample}LengthvsQualityScatterPlot_dot.svg", sample=config["sample"]),
                        expand("results/Trimming_QC/QC/ONT_NanoPlot/{sample}LengthvsQualityScatterPlot_kde.svg", sample=config["sample"]),
                        expand("results/Trimming_QC/QC/ONT_NanoPlot/{sample}Non_weightedHistogramReadlength.svg", sample=config["sample"]),
                        expand("results/Trimming_QC/QC/ONT_NanoPlot/{sample}Non_weightedLogTransformed_HistogramReadlength.svg", sample=config["sample"]),
                        expand("results/Trimming_QC/QC/ONT_NanoPlot/{sample}WeightedHistogramReadlength.svg", sample=config["sample"]),
                        expand("results/Trimming_QC/QC/ONT_NanoPlot/{sample}WeightedLogTransformed_HistogramReadlength.svg", sample=config["sample"]),
                        expand("results/Trimming_QC/QC/ONT_NanoPlot/{sample}Yield_By_Length.svg", sample=config["sample"])]

            for file in ont_nanoplot_out:
                out.append(file)

    ## Check the HiFi reads kmer evaluation
    if config['run_kmer_evaluation'].lower() == "yes":
        out.append(expand("results/Trimming_QC/Meryl_DB/{sample}.trimmed.meryl.histo", sample=config["sample"]))

        genomescope2_out = [expand("results/Assembly/Genome_Stats/HiFi_GenomeScope2/{sample}_linear_plot.png", sample=config["sample"]),
                            expand("results/Assembly/Genome_Stats/HiFi_GenomeScope2/{sample}_log_plot.png", sample=config["sample"]),
                            expand("results/Assembly/Genome_Stats/HiFi_GenomeScope2/{sample}_transformed_linear_plot.png", sample=config["sample"]),
                            expand("results/Assembly/Genome_Stats/HiFi_GenomeScope2/{sample}_transformed_log_plot.png", sample=config["sample"]),
                            expand("results/Assembly/Genome_Stats/HiFi_GenomeScope2/{sample}_model.txt", sample=config["sample"]),]
        
        for file in genomescope2_out:
            out.append(file)

        ## Check SmudgePlot Files
        out.append(expand("results/Assembly/Genome_Stats/Smudgeplot/{sample}_smudgeplot_log10.pdf", sample=config["sample"]))
        out.append(expand("results/Assembly/Genome_Stats/Smudgeplot/{sample}_smudgeplot.pdf", sample=config["sample"]))

        ## Check KAT GCP Plot
        out.append(expand("results/Assembly/Genome_Stats/KAT/{sample}.mx.png", sample=config["sample"]))
        
    ## Check Assembly outputs
    if config['run_asm'].lower() == "yes":
        if config['solo_asm'].lower() == 'yes':
            out.append(expand("results/Assembly/Contigging/Solo_Asm/{sample}.p_ctg.gfa", sample=config["sample"]))
            out.append(expand("results/Assembly/Contigging/Solo_Asm/{sample}.a_ctg.gfa", sample=config["sample"]))

        if config["hic_r1"] and config["hic_r2"]:
            out.extend(expand("results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap{hap}.p_ctg.gfa", sample=config["sample"], hap=["1", "2"]))

        ## Check Assembly stats files
        if config['solo_asm'].lower() == 'yes':
            out.append(expand("results/Assembly/Contigging/Solo_Asm/{sample}.p_ctg.fa", sample=config["sample"]))
            out.append(expand("results/Assembly/Genome_Stats/GFAstats/{sample}.p_ctg.fa.stats", sample=config["sample"]))
            out.append(expand("results/Assembly/Contigging/Solo_Asm/{sample}.a_ctg.fa", sample=config["sample"]))
            out.append(expand("results/Assembly/Genome_Stats/GFAstats/{sample}.a_ctg.fa.stats", sample=config["sample"]))

        if config["hic_r1"] and config["hic_r2"]:
            out.extend(expand("results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap{hap}.p_ctg.fa", sample=config["sample"], hap=["1", "2"]))
            out.extend(expand("results/Assembly/Genome_Stats/GFAstats/{sample}.hic.hap{hap}.p_ctg.fa.stats", sample=config["sample"], hap=["1", "2"]))

        ## Check Assembly evaluation
        if config['solo_asm'].lower() == 'yes':
         out.extend(expand("results/Assembly/Genome_Stats/Merqury/Solo_Asm/{sample}.completeness.stats", sample=config["sample"]))

        if config["hic_r1"] and config["hic_r2"]:
            out.extend(expand("results/Assembly/Genome_Stats/Merqury/Phased_Asm/{sample}.completeness.stats", sample=config["sample"]))

        # ## Check completeness
        # out.extend(expand("logs/{sample}_buscodb.check", sample=config["sample"]))

        if config['solo_asm'].lower() == 'yes':
            out.extend(expand("results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Primary/{sample}.summary.txt", sample=config["sample"]))
            out.extend(expand("results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Alt/{sample}.summary.txt", sample=config["sample"]))

        if config["hic_r1"] and config["hic_r2"]:
            out.extend(expand("results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap1/{sample}.summary.txt", sample=config["sample"]))
            out.extend(expand("results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap2/{sample}.summary.txt", sample=config["sample"]))

        ## Check SnailPlot files
        if config['solo_asm'].lower() == 'yes':
            out.append(expand("results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Primary/{sample}.full_table_busco_format_edit.tsv", sample=config["sample"]))
            out.append(expand("results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/00-Solo-Hap1/{sample}.bloobtools.create.check", sample=config["sample"]))
            out.append(expand("results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/00-Solo-Hap1/{sample}_Solo_Hap1.snail.png", sample=config["sample"]))

            out.append(expand("results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Alt/{sample}.full_table_busco_format_edit.tsv", sample=config["sample"]))
            out.append(expand("results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/01-Solo-Hap2/{sample}.bloobtools.create.check", sample=config["sample"]))
            out.append(expand("results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/01-Solo-Hap2/{sample}_Solo_Hap2.snail.png", sample=config["sample"]))

        if config["hic_r1"] and config["hic_r2"]:
            out.append(expand("results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap1/{sample}.full_table_busco_format_edit.tsv", sample=config["sample"]))
            out.append(expand("results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap1/{sample}.bloobtools.create.check", sample=config["sample"]))
            out.append(expand("results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap1/{sample}_Phased_Hap1.snail.png", sample=config["sample"]))

            out.append(expand("results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap2/{sample}.full_table_busco_format_edit.tsv", sample=config["sample"]))
            out.append(expand("results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap2/{sample}.bloobtools.create.check", sample=config["sample"]))
            out.append(expand("results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap2/{sample}_Phased_Hap2.snail.png", sample=config["sample"]))
    
        ## Check MitoHifi files
        out.append("resources/mitohifi.sif")

        ## Check Mitogenomes references
        out.append(expand("resources/{sample}.reference.fasta", sample=sampĺe_name))
        out.append(expand("resources/{sample}.reference.gb", sample=sampĺe_name))

        if config['solo_asm'].lower() == 'yes':
            out.append(expand("results/Assembly/Mitogenome/Solo_Asm/{sample}.contigs_stats.tsv", sample=config["sample"]))

        if config["hic_r1"] and config["hic_r2"]:
            out.append(expand("results/Assembly/Mitogenome/Phased_Asm/{sample}.contigs_stats.tsv", sample=config["sample"]))
    
    
    if config['gxdb'] and config['run_asm'].lower() == "yes":
        out.append("resources/fcs-gx.sif")

        ## Check FCS-Adaptors output
        if config['solo_asm'].lower() == 'yes':
            out.append(expand("results/Decontamination/FCS-Adaptor/Solo_Asm_Primary/cleaned_sequences/{sample}.p_ctg.fa", sample=config["sample"]))
            out.append(expand("results/Decontamination/FCS-Adaptor/Solo_Asm_Alt/cleaned_sequences/{sample}.a_ctg.fa", sample=config["sample"]))
        
        if config["hic_r1"] and config["hic_r2"]:
            out.append(expand("results/Decontamination/FCS-Adaptor/Phased_Asm_Hap1/cleaned_sequences/{sample}.hic.hap1.p_ctg.fa", sample=config["sample"]))
            out.append(expand("results/Decontamination/FCS-Adaptor/Phased_Asm_Hap2/cleaned_sequences/{sample}.hic.hap2.p_ctg.fa", sample=config["sample"]))
        
        ## Check FCS-GX output
        if config['solo_asm'].lower() == 'yes':
            out.append(expand("results/Decontamination/Contaminants/Solo_Asm_Primary/{sample}.p_ctg.screen.check", sample=config["sample"]))
            out.append(expand("results/Decontamination/Contaminants/Solo_Asm_Primary/{sample}.p_ctg.clean.fasta", sample=config["sample"]))
            out.append(expand("results/Decontamination/Contaminants/Solo_Asm_Alt/{sample}.a_ctg.screen.check", sample=config["sample"]))
            out.append(expand("results/Decontamination/Contaminants/Solo_Asm_Alt/{sample}.a_ctg.clean.fasta", sample=config["sample"]))
            
        if config["hic_r1"] and config["hic_r2"]:
            out.append(expand("results/Decontamination/Contaminants/Phased_Asm_Hap1/{sample}.hic.hap1.p_ctg.screen.check", sample=config["sample"]))
            out.append(expand("results/Decontamination/Contaminants/Phased_Asm_Hap2/{sample}.hic.hap2.p_ctg.screen.check", sample=config["sample"]))
            out.append(expand("results/Decontamination/Contaminants/Phased_Asm_Hap1/{sample}.hic.hap1.p_ctg.clean.fasta", sample=config["sample"]))
            out.append(expand("results/Decontamination/Contaminants/Phased_Asm_Hap2/{sample}.hic.hap2.p_ctg.clean.fasta", sample=config["sample"]))

    
    if config['run_auto_scaffolding'].lower() == 'yes':
        if config["hic_r1"] and config["hic_r2"]:

            if config['gxdb']:
                out.append(expand("results/Decontamination/Contaminants/Phased_Asm_Hap1/{sample}.hic.hap1.p_ctg.clean.fasta.bwt.2bit.64", sample=config["sample"]))
                out.append(expand("results/Decontamination/Contaminants/Phased_Asm_Hap2/{sample}.hic.hap2.p_ctg.clean.fasta.bwt.2bit.64", sample=config["sample"]))
        
            else:
                out.append(expand("results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa.bwt.2bit.64", sample=config["sample"]))
                out.append(expand("results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa.bwt.2bit.64", sample=config["sample"]))

            ## Check HiC initial contact maps
            out.append(expand("results/Scaffolding/Initial_Contacts/Hap1/{sample}.R1.bam", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Initial_Contacts/Hap1/{sample}.R2.bam", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Initial_Contacts/Hap2/{sample}.R1.bam", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Initial_Contacts/Hap2/{sample}.R2.bam", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Initial_Contacts/Hap1/{sample}.R1.bam.flagstats", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Initial_Contacts/Hap1/{sample}.R2.bam.flagstats", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Initial_Contacts/Hap2/{sample}.R1.bam.flagstats", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Initial_Contacts/Hap2/{sample}.R2.bam.flagstats", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Initial_Contacts/Hap1/{sample}.merged.bam", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Initial_Contacts/Hap2/{sample}.merged.bam", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Initial_Contacts/Hap1/{sample}.merged.bam.flagstats", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Initial_Contacts/Hap1/{sample}.merged.bam.hic_stats", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Initial_Contacts/Hap2/{sample}.merged.bam.flagstats", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Initial_Contacts/Hap2/{sample}.merged.bam.hic_stats", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Initial_Contacts/Hap1/{sample}.merged.bam.pretext", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Initial_Contacts/Hap2/{sample}.merged.bam.pretext", sample=config["sample"]))

            ## Check YAHS output
            out.append(expand("results/Scaffolding/YAHS_Scaffolding/Hap1/{sample}.yahs_scaffolds_final.fa", sample=config["sample"]))
            out.append(expand("results/Scaffolding/YAHS_Scaffolding/Hap2/{sample}.yahs_scaffolds_final.fa", sample=config["sample"]))
            
            ## Check HiC final contact maps
            out.append(expand("results/Scaffolding/Final_Contacts/Hap1/{sample}.R1.bam", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Final_Contacts/Hap1/{sample}.R2.bam", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Final_Contacts/Hap2/{sample}.R1.bam", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Final_Contacts/Hap2/{sample}.R2.bam", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Final_Contacts/Hap1/{sample}.R1.bam.flagstats", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Final_Contacts/Hap1/{sample}.R2.bam.flagstats", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Final_Contacts/Hap2/{sample}.R1.bam.flagstats", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Final_Contacts/Hap2/{sample}.R2.bam.flagstats", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Final_Contacts/Hap1/{sample}.merged.bam", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Final_Contacts/Hap2/{sample}.merged.bam", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Final_Contacts/Hap1/{sample}.merged.bam.flagstats", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Final_Contacts/Hap1/{sample}.merged.bam.hic_stats", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Final_Contacts/Hap2/{sample}.merged.bam.flagstats", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Final_Contacts/Hap2/{sample}.merged.bam.hic_stats", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Final_Contacts/Hap1/{sample}.merged.bam.pretext", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Final_Contacts/Hap2/{sample}.merged.bam.pretext", sample=config["sample"]))

        ## Check Scaffolding stats files
        if config["hic_r1"] and config["hic_r2"]:
            out.extend(expand("results/Scaffolding/Scaffolding_stats/GFAstats/{sample}.yahs_scaffolds_hap{hap}.fa.stats", sample=config["sample"], hap=["1", "2"]))
    
        ## Check Scaffolding completeness
        if config["hic_r1"] and config["hic_r2"]:
            out.extend(expand("results/Scaffolding/Scaffolding_stats/Compleasm/Hap1/{sample}.summary.txt", sample=config["sample"]))
            out.extend(expand("results/Scaffolding/Scaffolding_stats/Compleasm/Hap2/{sample}.summary.txt", sample=config["sample"]))

        ## Check Scaffolding SnailPlot files
        if config["hic_r1"] and config["hic_r2"]:
            out.append(expand("results/Scaffolding/Scaffolding_stats/Compleasm/Hap1/{sample}.full_table_busco_format_edit.tsv", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Scaffolding_stats/SnailPlot/Hap1/{sample}.bloobtools.create.check", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Scaffolding_stats/SnailPlot/Hap1/{sample}_Scaffolding_Hap1.snail.png", sample=config["sample"]))

            out.append(expand("results/Scaffolding/Scaffolding_stats/Compleasm/Hap2/{sample}.full_table_busco_format_edit.tsv", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Scaffolding_stats/SnailPlot/Hap2/{sample}.bloobtools.create.check", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Scaffolding_stats/SnailPlot/Hap2/{sample}_Scaffolding_Hap2.snail.png", sample=config["sample"]))
        
        ## Check Scaffolding Evaluation
        if config["hic_r1"] and config["hic_r2"]:
            out.append(expand("results/Scaffolding/Scaffolding_stats/Merqury/{sample}.completeness.stats", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Scaffolding_stats/Merqury/{sample}.{sample}.yahs_hap1.spectra-cn.st.png", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Scaffolding_stats/Merqury/{sample}.{sample}.yahs_hap2.spectra-cn.st.png", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Scaffolding_stats/Merqury/{sample}.spectra-asm.st.png", sample=config["sample"]))
            out.append(expand("results/Scaffolding/Scaffolding_stats/Merqury/{sample}.spectra-cn.st.png", sample=config["sample"]))

    return out

# ### Function to create a list of files for checking before remove buscodb folder
# def files_remove_buscodb():
    
#     files = []

#     if config['hic_r1'] and config['hic_r2']:
#         files.append(expand("results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap1/{sample}.summary.txt", sample=config["sample"]))
#         files.append(expand("results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap2/{sample}.summary.txt", sample=config["sample"]))
    
#         if config['run_auto_scaffolding'].lower() == "yes":
#             files.append(expand("results/Scaffolding/Scaffolding_stats/Compleasm/Hap1/{sample}.summary.txt", sample=config["sample"]))
#             files.append(expand("results/Scaffolding/Scaffolding_stats/Compleasm/Hap2/{sample}.summary.txt", sample=config["sample"]))
    
#     elif config['solo_asm'].lower() == 'yes':
#         files.append(expand("results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Primary/{sample}.summary.txt", sample=config["sample"]))
#         files.append(expand("results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Alt/{sample}.summary.txt", sample=config["sample"]))
        

#     return files