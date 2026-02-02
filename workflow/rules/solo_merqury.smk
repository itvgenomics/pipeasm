rule solo_merqury:
	input:
		hap1="results/Assembly/Contigging/Solo_Asm/{sample}.p_ctg.fa",
		hap2="results/Assembly/Contigging/Solo_Asm/{sample}.a_ctg.fa",
		fastk="results/Trimming_QC/FastK/{sample}.ktab"
	output:
		"results/Assembly/Genome_Stats/MerquryFK/Solo_Asm/{sample}.completeness.stats",
		"results/Assembly/Genome_Stats/MerquryFK/Solo_Asm/{sample}.{sample}.p_ctg.spectra-cn.st.png",
		"results/Assembly/Genome_Stats/MerquryFK/Solo_Asm/{sample}.{sample}.a_ctg.spectra-cn.st.png",
		"results/Assembly/Genome_Stats/MerquryFK/Solo_Asm/{sample}.spectra-asm.st.png",
		"results/Assembly/Genome_Stats/MerquryFK/Solo_Asm/{sample}.spectra-cn.st.png"
	log:
		"logs/{sample}.solo_merqury.log"
	benchmark:
		"benchmarks/{sample}.solo_merqury.txt"
	params:
		outdir="results/Assembly/Genome_Stats/MerquryFK/Solo_Asm"
	singularity:
		f"{config["sif_dir"]}/fkutils.sif"
	shell:
		"""
		mkdir -p {params.outdir} && \
		MerquryFK -v -X250 -Ptmp -lfs -T{threads} {input.fastk} {input.hap1} {input.hap2} {params.outdir}/{wildcards.sample} >> {log} 2>&1 && \
		mv {params.outdir}/{wildcards.sample}.hap1.spectra-cn.st.png {params.outdir}/{wildcards.sample}.{wildcards.sample}.p_ctg.spectra-cn.st.png && \
		mv {params.outdir}/{wildcards.sample}.hap2.spectra-cn.st.png {params.outdir}/{wildcards.sample}.{wildcards.sample}.a_ctg.spectra-cn.st.png
		"""
