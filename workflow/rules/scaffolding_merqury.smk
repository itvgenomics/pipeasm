rule scaffolding_merqury:
	input:
		meryl="results/Trimming_QC/Meryl_DB/{sample}.trimmed.union.meryl",
		hap1="results/Scaffolding/YAHS_Scaffolding/Hap1/{sample}.yahs_scaffolds_final.fa",
		hap2="results/Scaffolding/YAHS_Scaffolding/Hap2/{sample}.yahs_scaffolds_final.fa",
	output:
		"results/Scaffolding/Scaffolding_stats/Merqury/{sample}.completeness.stats",
		"results/Scaffolding/Scaffolding_stats/Merqury/{sample}.{sample}.yahs_hap1.spectra-cn.st.png",
		"results/Scaffolding/Scaffolding_stats/Merqury/{sample}.{sample}.yahs_hap2.spectra-cn.st.png",
		"results/Scaffolding/Scaffolding_stats/Merqury/{sample}.spectra-asm.st.png",
		"results/Scaffolding/Scaffolding_stats/Merqury/{sample}.spectra-cn.st.png"
	log:
		"logs/{sample}.scaffolding_merqury.log"
	benchmark:
		"benchmarks/{sample}.scaffolding_merqury.txt"
	threads:
		config["threads"]
	params:
		"results/Scaffolding/Scaffolding_stats/Merqury"
	singularity:
		"docker://itvdsbioinfo/merqury:1.3"
	shell:
		"""
		bash -c 'export OMP_NUM_THREADS={threads} && \
		export MERQURY=/opt/conda/share/merqury && \
		mkdir -p {params} && \
        cp {input.hap1} results/Scaffolding/Scaffolding_stats/Merqury/{wildcards.sample}.yahs_hap1.fa && \
        cp {input.hap2} results/Scaffolding/Scaffolding_stats/Merqury/{wildcards.sample}.yahs_hap2.fa && \
		cd {params} && \
		merqury.sh ../../../../{input.meryl} {wildcards.sample}.yahs_hap1.fa {wildcards.sample}.yahs_hap2.fa {wildcards.sample} >> ../../../../{log} 2>&1 && \
        rm {wildcards.sample}.yahs_hap1.fa {wildcards.sample}.yahs_hap2.fa'
		"""