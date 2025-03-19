rule trimming_pacbio:
	input:
		"results/Trimming_QC/HiFi/{sample}.fastq.gz"
	output:
		"results/Trimming_QC/HiFi/{sample}.trimmed.fastq.gz"
	threads:
		config["threads"]
	params:
		adaptors=config['cutadapt']['adaptors'],
		filtering=config['cutadapt']['filtering'],
		adap_find=config['cutadapt']['adap_find']
	log:
		"logs/{sample}.trimming_pacbio.log"
	benchmark:
		"benchmarks/{sample}.trimming_pacbio.txt"
	singularity:
		f"{config["sif_dir"]}/cutadapt.sif"
	shell:
		"cutadapt -j {threads} {params} -o {output} {input} >> {log} 2>&1"


rule qc_trim_pacbio:
	input:
		"results/Trimming_QC/HiFi/{sample}.trimmed.fastq.gz"
	output:
		expand("results/Trimming_QC/QC/HiFi_FastQC/{{sample}}.trimmed_fastqc.{ext}", ext=['zip','html'])
	threads:
		config["software_threads"]["fastqc"]
	log:
		"logs/{sample}.qc_trim_pacbio.log"
	benchmark:
		"benchmarks/{sample}.qc_trim_pacbio.txt"
	singularity:
		f"{config["sif_dir"]}/fastqc.sif"
	params:
		adapt_txt="resources/pacbio_adapters.txt",
		outdir="results/Trimming_QC/QC/HiFi_FastQC/"
	shell:
		"""
		mkdir -p tmp && \
		(fastqc -t {threads} --dir tmp -a {params.adapt_txt} --svg \
		-o {params.outdir} {input} >> {log} 2>&1 \
		&& rm -r .cache .java)
		"""

rule qc_nanoplot_pacbio:
	input:
		"results/Trimming_QC/HiFi/{sample}.trimmed.fastq.gz"
	output:
		"results/Trimming_QC/QC/HiFi_NanoPlot/{sample}LengthvsQualityScatterPlot_dot.svg",
		"results/Trimming_QC/QC/HiFi_NanoPlot/{sample}LengthvsQualityScatterPlot_kde.svg",
		"results/Trimming_QC/QC/HiFi_NanoPlot/{sample}Non_weightedHistogramReadlength.svg",
		"results/Trimming_QC/QC/HiFi_NanoPlot/{sample}Non_weightedLogTransformed_HistogramReadlength.svg",
		"results/Trimming_QC/QC/HiFi_NanoPlot/{sample}WeightedHistogramReadlength.svg",
		"results/Trimming_QC/QC/HiFi_NanoPlot/{sample}WeightedLogTransformed_HistogramReadlength.svg",
		"results/Trimming_QC/QC/HiFi_NanoPlot/{sample}Yield_By_Length.svg",
		"results/Trimming_QC/QC/HiFi_NanoPlot/{sample}NanoPlot-report.html"
	threads:
		config["software_threads"]["nanoplot"]
	params:
		args= config['nanoplot']['args'],
		plots= config['nanoplot']['plots'],
		outdir="results/Trimming_QC/QC/HiFi_NanoPlot"
	log:
		"logs/{sample}.qc_nanoplot_pacbio.log"
	benchmark:
		"benchmarks/{sample}.qc_nanoplot_pacbio.txt"
	singularity:
		f"{config["sif_dir"]}/nanoplot.sif"
	shell:
		"""
		NanoPlot -t {threads} \
		-o {params.outdir} -p {wildcards.sample} \
		--fastq {input} {params.args} {params.plots} >> {log} 2>&1
		"""
