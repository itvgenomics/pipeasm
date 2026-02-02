rule phased_merqury:
	input:
		hap1="results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa",
		hap2="results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa",
		fastk="results/Trimming_QC/FastK/{sample}.ktab"
	output:
		"results/Assembly/Genome_Stats/MerquryFK/Phased_Asm/{sample}.completeness.stats",
		"results/Assembly/Genome_Stats/MerquryFK/Phased_Asm/{sample}.{sample}.hic.hap1.p_ctg.spectra-cn.st.png",
		"results/Assembly/Genome_Stats/MerquryFK/Phased_Asm/{sample}.{sample}.hic.hap2.p_ctg.spectra-cn.st.png",
		"results/Assembly/Genome_Stats/MerquryFK/Phased_Asm/{sample}.spectra-asm.st.png",
		"results/Assembly/Genome_Stats/MerquryFK/Phased_Asm/{sample}.spectra-cn.st.png"
	log:
		"logs/{sample}.phased_merqury.log"
	benchmark:
		"benchmarks/{sample}.phased_merqury.txt"
	params:
		outdir="results/Assembly/Genome_Stats/MerquryFK/Phased_Asm"
	singularity:
		f"{config["sif_dir"]}/fkutils.sif"
	shell:
		"""
		mkdir -p {params.outdir} && \
		MerquryFK -v -X250 -Ptmp -lfs -T{threads} {input.fastk} {input.hap1} {input.hap2} {params.outdir}/{wildcards.sample} >> {log} 2>&1
		"""
