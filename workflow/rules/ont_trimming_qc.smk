rule trimming_ont:
	input:
		"results/Trimming_QC/ONT/{sample}_ONT.fastq.gz"
	output:
		"results/Trimming_QC/ONT/{sample}_ONT.trimmed.fastq.gz"
	threads:
		config["threads"]
	log:
		"logs/{sample}.trimming_ont.log"
	benchmark:
		"benchmarks/{sample}.trimming_ont.txt"
	singularity:
		f"{config["sif_dir"]}/dorado.sif"
	shell:
		"""
        dorado trim {input} --emit-fastq -t {threads} > results/Trimming_QC/ONT/{wildcards.sample}_ONT.trimmed.fastq 2> {log} && \
        gzip results/Trimming_QC/ONT/{wildcards.sample}_ONT.trimmed.fastq
        """

rule qc_nanoplot_ont:
	input:
		"results/Trimming_QC/ONT/{sample}_ONT.trimmed.fastq.gz"
	output:
		"results/Trimming_QC/QC/ONT_NanoPlot/{sample}LengthvsQualityScatterPlot_dot.svg",
		"results/Trimming_QC/QC/ONT_NanoPlot/{sample}LengthvsQualityScatterPlot_kde.svg",
		"results/Trimming_QC/QC/ONT_NanoPlot/{sample}Non_weightedHistogramReadlength.svg",
		"results/Trimming_QC/QC/ONT_NanoPlot/{sample}Non_weightedLogTransformed_HistogramReadlength.svg",
		"results/Trimming_QC/QC/ONT_NanoPlot/{sample}WeightedHistogramReadlength.svg",
		"results/Trimming_QC/QC/ONT_NanoPlot/{sample}WeightedLogTransformed_HistogramReadlength.svg",
		"results/Trimming_QC/QC/ONT_NanoPlot/{sample}Yield_By_Length.svg",
		"results/Trimming_QC/QC/ONT_NanoPlot/{sample}NanoPlot-report.html"
	threads:
		config["software_threads"]["nanoplot"]
	params:
		args= config['nanoplot']['args'],
		plots= config['nanoplot']['plots'],
		outdir="results/Trimming_QC/QC/ONT_NanoPlot"
	log:
		"logs/{sample}.qc_nanoplot_ont.log"
	benchmark:
		"benchmarks/{sample}.qc_nanoplot_ont.txt"
	singularity:
		f"{config["sif_dir"]}/nanoplot.sif"
	shell:
		"""
		NanoPlot -t {threads} \
		-o {params.outdir} -p {wildcards.sample} \
		--fastq {input} {params.args} {params.plots} >> {log} 2>&1
		"""
