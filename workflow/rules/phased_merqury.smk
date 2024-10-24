rule phased_merqury:
	input:
		meryl="results/Trimming_QC/Meryl_DB/{sample}.trimmed.union.meryl",
		hap1="results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap1.p_ctg.fa",
		hap2="results/Assembly/Contigging/Phased_Asm/{sample}.hic.hap2.p_ctg.fa",
	output:
		"results/Assembly/Genome_Stats/Merqury/Phased_Asm/{sample}.completeness.stats",
		"results/Assembly/Genome_Stats/Merqury/Phased_Asm/{sample}.{sample}.hic.hap1.p_ctg.spectra-cn.st.png",
		"results/Assembly/Genome_Stats/Merqury/Phased_Asm/{sample}.{sample}.hic.hap2.p_ctg.spectra-cn.st.png",
		"results/Assembly/Genome_Stats/Merqury/Phased_Asm/{sample}.spectra-asm.st.png",
		"results/Assembly/Genome_Stats/Merqury/Phased_Asm/{sample}.spectra-cn.st.png"
	log:
		"logs/{sample}.phased_merqury.log"
	benchmark:
		"benchmarks/{sample}.phased_merqury.txt"
	threads:
		config["threads"]
	params:
		"results/Assembly/Genome_Stats/Merqury/Phased_Asm"
	singularity:
		"docker://itvdsbioinfo/merqury:1.3"
	shell:
		"""
		bash -c 'export OMP_NUM_THREADS={threads} && \
		export MERQURY=/opt/conda/share/merqury && \
		mkdir -p {params} && \
		cd {params} && \
		merqury.sh ../../../../../{input.meryl} ../../../../../{input.hap1} ../../../../../{input.hap2} {wildcards.sample} >> ../../../../../{log} 2>&1'
		"""