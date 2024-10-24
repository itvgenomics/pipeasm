rule solo_merqury:
	input:
		meryl="results/Trimming_QC/Meryl_DB/{sample}.trimmed.union.meryl",
		hap1="results/Assembly/Contigging/Solo_Asm/{sample}.p_ctg.fa",
		hap2="results/Assembly/Contigging/Solo_Asm/{sample}.a_ctg.fa"
	output:
		"results/Assembly/Genome_Stats/Merqury/Solo_Asm/{sample}.completeness.stats",
		"results/Assembly/Genome_Stats/Merqury/Solo_Asm/{sample}.{sample}.p_ctg.spectra-cn.st.png",
		"results/Assembly/Genome_Stats/Merqury/Solo_Asm/{sample}.{sample}.a_ctg.spectra-cn.st.png",
		"results/Assembly/Genome_Stats/Merqury/Solo_Asm/{sample}.spectra-asm.st.png",
		"results/Assembly/Genome_Stats/Merqury/Solo_Asm/{sample}.spectra-cn.st.png"
	log:
		"logs/{sample}.solo_merqury.log"
	benchmark:
		"benchmarks/{sample}.solo_merqury.txt"
	threads:
		config["threads"]
	params:
		"results/Assembly/Genome_Stats/Merqury/Solo_Asm"
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