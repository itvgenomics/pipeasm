rule scaffolding_merqury:
	input:
		hap1="results/Scaffolding/YAHS_Scaffolding/Hap1/{sample}.yahs_scaffolds_final.fa",
		hap2="results/Scaffolding/YAHS_Scaffolding/Hap2/{sample}.yahs_scaffolds_final.fa",
		fastk="results/Trimming_QC/FastK/{sample}.ktab"
	output:
		"results/Scaffolding/Scaffolding_stats/MerquryFK/{sample}.completeness.stats",
		"results/Scaffolding/Scaffolding_stats/MerquryFK/{sample}.{sample}.yahs_hap1.spectra-cn.st.png",
		"results/Scaffolding/Scaffolding_stats/MerquryFK/{sample}.{sample}.yahs_hap2.spectra-cn.st.png",
		"results/Scaffolding/Scaffolding_stats/MerquryFK/{sample}.spectra-asm.st.png",
		"results/Scaffolding/Scaffolding_stats/MerquryFK/{sample}.spectra-cn.st.png"
	log:
		"logs/{sample}.scaffolding_merqury.log"
	benchmark:
		"benchmarks/{sample}.scaffolding_merqury.txt"
	params:
		outdir="results/Scaffolding/Scaffolding_stats/MerquryFK"
	singularity:
		f"{config["sif_dir"]}/fkutils.sif"
	shell:
		"""
		mkdir -p {params.outdir} && \
		MerquryFK -v -X250 -Ptmp -T{threads} {input.fastk} {input.hap1} {input.hap2} {params.outdir}/{wildcards.sample} >> {log} 2>&1 && \
		mv {params.outdir}/{wildcards.sample}.hap1.spectra-cn.st.png {params.outdir}/{wildcards.sample}.{wildcards.sample}.yahs_hap1.spectra-cn.st.png && \
		mv {params.outdir}/{wildcards.sample}.hap2.spectra-cn.st.png {params.outdir}/{wildcards.sample}.{wildcards.sample}.yahs_hap2.spectra-cn.st.png
		"""
