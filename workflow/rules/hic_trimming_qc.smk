rule trimming_hic:
	input:
		r1="results/Trimming_QC/HiC/{sample}_R1.fastq.gz",
		r2="results/Trimming_QC/HiC/{sample}_R2.fastq.gz"
	output:
		r1="results/Trimming_QC/HiC/{sample}_R1.trimmed_paired.fastq.gz",
		r2="results/Trimming_QC/HiC/{sample}_R2.trimmed_paired.fastq.gz"
	threads:
		config["software_threads"]["fastp"]
	params:
		flags= config['fastp']['flags'],
		min_len= config['fastp']['min_len'],
		crop= config['fastp']['crop'],
		quality= config['fastp']['quality']
	log:
		"logs/{sample}.trimming_hic.log"
	benchmark:
		"benchmarks/{sample}.trimming_hic.txt"
	singularity:
		f"{config["sif_dir"]}/fastp.sif"
	shell:
		"""
		fastp --in1 {input.r1} --in2 {input.r2} \
		{params} --thread {threads} --html results/Trimming_QC/HiC/fastp.html --json results/Trimming_QC/HiC/fastp.json \
		--out1 {output.r1} --out2 {output.r2} >> {log} 2>&1
		"""

rule qc_trim_hic:
	input:
		expand("results/Trimming_QC/HiC/{{sample}}_R{pair}.trimmed_paired.fastq.gz", pair=['1','2'])
	output:
		expand("results/Trimming_QC/QC/HiC_FastQC/{{sample}}_R{pair}.trimmed_paired_fastqc.{ext}", pair=['1','2'], ext=['zip','html'])
	threads:
		config["software_threads"]["fastqc"]
	log:
		"logs/{sample}.qc_trim_hic.log"
	benchmark:
		"benchmarks/{sample}.qc_trim_hic.txt"
	singularity:
		f"{config["sif_dir"]}/fastqc.sif"
	params:
		outdir="results/Trimming_QC/QC/HiC_FastQC/"
	shell:
		"""
		mkdir -p tmp && \
		(fastqc -t {threads} --dir tmp --svg \
		-o {params.outdir} {input} >> {log} 2>&1 \
		&& rm -r .cache .java)
		"""
