# Pipeasm: a tool for automated large chromosome-scale genome assembly and evaluation

[![DOI](https://img.shields.io/badge/DOI-10.1093%2Fbioadv%2Fvbaf326-blue)](https://doi.org/10.1093/bioadv/vbaf326)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Snakemake](https://img.shields.io/badge/snakemake-8.5.1-brightgreen.svg)](https://snakemake.readthedocs.io)
[![Release](https://img.shields.io/github/v/tag/itvgenomics/pipeasm?label=release)](https://github.com/itvgenomics/pipeasm/tags)

<p align="center">
  <img src="Pipeasm.png" alt="Pipeasm Logo"/>
</p>

Pipeasm streamlines the VGP and DToL assembly pipelines into a single automated workflow. It is written in [Snakemake](https://snakemake.readthedocs.io) and runs every tool inside Singularity containers, so the only software you need to install is Snakemake and Singularity.

## Table of Contents

- [Quick Start](#quick-start)
- [How Pipeasm Works](#how-pipeasm-works)
- [Installation](#installation)
- [Hardware Requirements](#hardware-requirements)
- [Inputs and Configuration](#inputs-and-configuration)
- [Running Pipeasm](#running-pipeasm)
- [Outputs](#outputs)
- [Test Dataset](#test-dataset)
- [Troubleshooting](#troubleshooting)
- [How to Cite](#how-to-cite)
- [Changelog](#changelog)

## Quick Start

```bash
# 1. Create an environment with Snakemake and Singularity
conda create -n pipeasm -c conda-forge -c bioconda python=3.11 singularity=3.8.5 snakemake=8.5.1
conda activate pipeasm

# 2. Clone the repository
git clone https://github.com/itvgenomics/pipeasm.git
cd pipeasm

# 3. Fill in config/config.yaml (sample, reads, taxid, buscodb, genome_size, sif_dir, ...)

# 4. Dry run: checks inputs, pulls the container images and builds the job graph
bash Pipeasm.sh -d . -t 32 -np

# 5. Run
bash Pipeasm.sh -d . -t 32
```

## How Pipeasm Works

Pipeasm integrates a suite of tools to process and assemble genomic data, producing comprehensive reports and statistics. The workflow is organized into six main steps:

<p align="center">
  <img src="workflow.png" alt="Workflow"/>
</p>

1. **Trimming and Quality Control (QC)**
   - Accepts **PacBio HiFi (CCS) long reads**, **Oxford Nanopore (ONT) long reads** and **Hi-C Illumina paired-end short reads**.
   - **Trimming**: `Cutadapt` for HiFi reads, `fastp` for Hi-C reads and `Dorado` for ONT reads, to remove adapters and primers.
   - **Quality assessment**: `FastQC` and `NanoPlot`.

2. **K-mer Profiling** (HiFi reads)
   - **K-mer counting**: `FastK`.
   - **Genome statistics**: `GeneScopeFK` estimates genome size, heterozygosity and repeat content.
   - **Ploidy**: `Smudgeplot` infers ploidy from k-mer pairs.
   - **K-mer spectra**: `KatGC` and `KatComp` (from `MerquryFK`) plot k-mer spectra to help identify genomic features and contaminants.

3. **Assembly**
   - `Hifiasm` produces haplotype-resolved (Hi-C phased) or primary/alternate assemblies. ONT reads, when provided, are integrated as ultra-long reads (`--ul`).
   - **Mitogenome** (optional): `MitoHiFi` recovers the mitochondrial genome.

4. **Assembly Statistics**
   - **Summary statistics**: `gfastats` converts the assembly graph to FASTA and reports contig counts, N50, etc.
   - **Completeness**: `Compleasm` assesses gene-space completeness against BUSCO lineages.
   - **K-mer QV and completeness**: `MerquryFK` compares assembled k-mers with the read k-mers.
   - **Visualization**: Snail plots from `BlobToolkit`.

5. **Decontamination** (optional)
   - `FCS-adaptor` removes adapters and `FCS-GX` removes other contaminants.

6. **Hi-C Mapping and Scaffolding** (requires Hi-C reads)
   - **Mapping**: `BWA-MEM2` indexes the assembly and maps the Hi-C reads.
   - **Filtering**: reads are sorted with `SAMtools`, then filtered and merged with `Bellerophon`.
   - **Scaffolding**: `YaHS` links contigs into chromosome-scale scaffolds.
   - **Contact maps**: `PretextMap` (ultra-resolution, `--ultraRes`) and `PretextSnapshot` produce Hi-C contact maps before and after scaffolding.

### Software and Containers

The images and their versions are listed in [`workflow/containers.yaml`](workflow/containers.yaml). Before each run, `Pipeasm.sh` pulls the missing ones from Docker Hub into `sif_dir` as `<name>_<tag>.sif` (e.g., `hic_mapping_1.1.sif`), so different versions can share the same `sif_dir` and an image is pulled again only when its version changes. The Dockerfiles for the `itvdsbioinfo/*` images are in [`container_files/dockerfiles`](container_files/dockerfiles).

If the machine running Pipeasm has no internet access, `Pipeasm.sh` stops and prints the `singularity pull` command for each missing image. Run them on another machine and copy the `.sif` files to `sif_dir`.

| Tool(s) | Version | Image |
| --- | --- | --- |
| Cutadapt | 4.4 | `pipecraft/cutadapt:4.4` |
| fastp | 0.23.4 | `staphb/fastp:0.23.4` |
| FastQC | 0.12.1 | `staphb/fastqc:0.12.1` |
| NanoPlot | 1.41.6 | `staphb/nanoplot:1.41.6` |
| Dorado | pinned build | `nanoporetech/dorado:shae36d1b49fe470a60e006afad90bedd2fc2774a89` |
| FastK / MerquryFK / GeneScopeFK / Smudgeplot | 1.2 / 1.2 / 1.0.0 / 0.5.3 | `itvdsbioinfo/fkutils:1.2` |
| Hifiasm | 0.25.0 | `itvdsbioinfo/hifiasm:0.25.0` |
| MitoHiFi | master | `ghcr.io/marcelauliano/mitohifi:master` |
| gfastats | 1.3.6 | `staphb/gfastats:1.3.6` |
| Compleasm | 0.2.8 | `itvdsbioinfo/compleasm:0.2.8` |
| BlobToolkit | 4.3.5 | `genomehubs/blobtoolkit:4.3.5` |
| FCS-adaptor / FCS-GX | latest | `ncbi/fcs-adaptor:latest`, `ncbi/fcs-gx:latest` |
| BWA-MEM2 / SAMtools / Bellerophon / PretextMap / PretextSnapshot | 2.2.1 / 1.19 / 1.0 / 0.2.4 / 0.1.1 | `itvdsbioinfo/hic_mapping:1.1` |
| YaHS | 1.2a.2 | `itvdsbioinfo/yahs:1.2a.2` |

## Installation

### Prerequisites
- A Linux operating system (e.g., Ubuntu, CentOS, Fedora).
- [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html).
- Python 3.11 or later (required by Snakemake 8). It is installed in the environment below.

### Snakemake and Singularity
Create a dedicated environment (recommended):

```bash
conda create -n pipeasm -c conda-forge -c bioconda python=3.11 singularity=3.8.5 snakemake=8.5.1
conda activate pipeasm
```

Or install them in an existing environment:

```bash
mamba install -c conda-forge -c bioconda singularity=3.8.5 snakemake=8.5.1
```

**Docker is optional.** Pipeasm pulls Docker Hub images through Singularity (`singularity pull docker://...`), so Docker is only needed if you want to rebuild the images from [`container_files/dockerfiles`](container_files/dockerfiles) (tested with Docker 24.0.7).

### Clone the Repository

```bash
git clone https://github.com/itvgenomics/pipeasm.git
cd pipeasm
```

### Additional Resources
- [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/)
- [Snakemake GitHub repository](https://github.com/snakemake/snakemake)

## Hardware Requirements

Pipeasm targets HPC nodes or large workstations. By default, the heavy rules use **31 threads**, so a node with at least **32 cores** is recommended.

The memory requested per rule in SLURM mode depends on `genome_size` (see `config/resources/{small,medium,large}.yaml`). The table shows the largest requests:

| `genome_size` | Genome size | Hifiasm assembly | FCS-GX (if enabled) | BWA-MEM2 index |
| --- | --- | --- | --- | --- |
| `small` | < 2 Gb | 200 GB | 480 GB | 40 GB |
| `medium` | 2–4 Gb | 400 GB | 480 GB | 120 GB |
| `large` | > 4 Gb | 800 GB | 480 GB | 240 GB |

- **Without decontamination** (`gxdb` empty), the peak memory is set by the Hifiasm assembly step.
- **With decontamination**, FCS-GX needs ~480–512 GB of RAM to hold its database in memory (see [FCS-GX Database](#fcs-gx-database)).
- **Hi-C contact maps** are built with `PretextMap --ultraRes`, which needs ~40 GB of RAM per map (50 GB is requested per job). Replace it with `--highRes` (~16 GB) or remove it in `pretext.map` (`config/config.yaml`) to save memory.
- **Disk**: plan for several times the size of the raw reads (trimmed reads, FastK tables, BAM files and intermediate assemblies), plus the container images in `sif_dir` and, if used, ~470 GB for the FCS-GX database.

**Reference runtimes** (from the [Pipeasm paper](https://doi.org/10.1093/bioadv/vbaf326); Intel Xeon Gold 6252, 3 TB RAM, 64 threads with 32 per rule):

| Species | Genome size | HiFi | Hi-C | Phased assembly + scaffolding |
| --- | --- | --- | --- | --- |
| *Taeniopygia guttata* (bTaeGut) | 1.1 Gb | 39.5 Gb | 99.0 Gb | ~14 h |
| *Pseudophryne corroboree* (aPseCor) | 8.9 Gb | 230.3 Gb | 381.0 Gb | ~80 h |

In general, an assembly ready for manual curation takes about 1–4 days.

## Inputs and Configuration

### Supported Input Combinations

HiFi reads are always required. Hi-C and ONT reads are optional and change which assemblies are produced:

| Inputs | `solo_asm` | `diff_species_hic` | Assemblies produced | Scaffolding |
| --- | --- | --- | --- | --- |
| HiFi | `Yes` | `No` | Primary + alternate (Hifiasm solo) | No |
| HiFi | `No` | `No` | None (only trimming/QC and k-mer profiling) | No |
| HiFi + Hi-C (same species) | `No` | `No` | Hap1 + Hap2 (Hi-C phased) | Hap1 and Hap2 |
| HiFi + Hi-C (same species) | `Yes` | `No` | Hap1 + Hap2 **and** primary + alternate | Hap1 and Hap2 |
| HiFi + Hi-C (different species) | any | `Yes` | Primary + alternate (`--primary`) | Primary and alternate |
| Any of the above + ONT | — | — | Same as above, with ONT reads added as ultra-long reads (`--ul`) | Same as above |

ONT-only and Hi-C-only runs are not supported.

### Raw Data
Each read type must be merged into a single gzip-compressed FASTQ file:

```
/path/to/hifi_reads/hifi_reads.fastq.gz
/path/to/hic_reads/r1_reads.fastq.gz
/path/to/hic_reads/r2_reads.fastq.gz
/path/to/ont_reads/ont_reads.fastq.gz
```

> **Note:** raw reads copied into the working directory are marked as temporary and deleted at the end of the run. Only the trimmed reads are kept.

### `config/config.yaml`

The config file is split into sections (sample and reads, assembly options, completeness and contaminants, mitogenome, containers, and advanced tool parameters). Each field has a comment explaining its accepted values. Use absolute paths, and write `"Yes"`/`"No"` in quotes.

| Field | Example | Description | Requirement |
| --- | --- | --- | --- |
| sample | `"bHypStr1"` | Unique sample identifier (letters, numbers and underscores), e.g., DToL ID. | Required |
| hifi_reads | `"/data/raw_data/hifi/bHypStr1.fastq.gz"` | Absolute path to the HiFi reads, merged into a single `.fastq.gz` file. | Required |
| hic_r1 / hic_r2 | `"/data/raw_data/hic/bHypStr1.R1.fastq.gz"` | Absolute paths to the Hi-C forward and reverse reads (`.fastq.gz`). Set both or neither. | Optional |
| ont_reads | `"/data/raw_data/ont/bHypStr1.ONT.fastq.gz"` | Absolute path to the Oxford Nanopore reads (`.fastq.gz`). | Optional |
| genome_size | `"medium"` | `small` (< 2 Gb), `medium` (2–4 Gb) or `large` (> 4 Gb). Selects the SLURM resource file and must be set in **both local and SLURM** modes. | Required |
| solo_asm | `"No"` | `"Yes"` for both solo and phased assemblies, `"No"` for phased only. Set `"Yes"` when there is no Hi-C. | Required |
| diff_species_hic | `"No"` | `"Yes"` if the Hi-C data comes from a different species. Hifiasm then runs in solo primary mode and YaHS scaffolds the primary and alternate assemblies. | Required |
| buscodb | `"aves"` | [BUSCO lineage](https://busco-data.ezlab.org/v5/data/lineages/) of the closest taxonomic group, without the `_odb10` suffix. | Required |
| taxid | `1202453` | [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) ID. Use a similar species if unavailable. Used by the decontamination step. | Required |
| gxdb | `"/databases/fcs-gx/gxdb"` | Path to the FCS-GX database. Leave empty to skip decontamination. | Optional |
| run_mitohifi | `"No"` | `"Yes"` to assemble the mitogenome with MitoHiFi, `"No"` to skip it. | Required |
| species | `"Hypocnemis striata"` | Scientific name, used to search for a mitochondrial reference. Must exist in the NCBI Taxonomy database. | Required if `run_mitohifi: "Yes"` |
| geneticcode | `2` | [Genetic code](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for mitogenome annotation, e.g., `2` for vertebrates. | Required if `run_mitohifi: "Yes"` |
| sif_dir | `"/data/singularity_images/pipeasm"` | Directory where the Singularity images are stored. | Required |

The **advanced tool parameters** at the end of the file (fastp, Cutadapt, NanoPlot, FastK, GeneScopeFK, Smudgeplot, Hifiasm, gfastats, Bellerophon, Pretext, YaHS) follow the VGP/DToL defaults. Each one has a comment explaining what it does. For example, change the YaHS `-e` enzyme motifs to match your Hi-C kit, or use `-l 0` in Hifiasm for highly homozygous genomes.

**Configuration check.** Before every run, `Pipeasm.sh` validates the config file with `workflow/scripts/check_config.py`. It reports all problems at once (missing required fields, files not found, relative paths, invalid `"Yes"`/`"No"` values, Hi-C R1 without R2, missing `species`/`geneticcode` for MitoHiFi, mismatched k-mer sizes, etc.) and stops before anything runs. The config file is never modified by Pipeasm.

### Threads and Memory

Threads and resources are kept in a few small files under `config/`, so each value is set in one place only:

| File | Used in | Contents |
| --- | --- | --- |
| `config/threads.yaml` | Local and SLURM | Threads per rule (rules not listed use 1 thread). |
| `config/local_params.yaml` | Local | General Snakemake settings for local runs. |
| `config/slurm_params.yaml` | SLURM | General SLURM settings and the default resources for every rule (2 GB, 72 h, the `-partition` queue, and an optional `slurm_account`). |
| `config/resources/{small,medium,large}.yaml` | SLURM | Memory per rule for each `genome_size`. Only rules that need more than the defaults are listed. |

`Pipeasm.sh` joins the files for the selected mode, fills in the placeholders in braces (`{ WORKDIR }`, `{ THREADS }`, `{ PARTITION }`; do not edit them) and writes `profiles/local/config.yaml` or `profiles/slurm/config.yaml`. Those generated files are overwritten on every run, so edit the files in `config/` instead.

Examples:

```yaml
# config/threads.yaml - give Hifiasm 64 threads
  phased_assembly: 64

# config/resources/large.yaml - more memory, another partition or a longer time limit for one rule
  phased_assembly: {mem_mb: 1000000}
  phased_fcsgx: {mem_mb: 480000, slurm_partition: "highmem"}
  yahs_scaffolding_hap1: {mem_mb: 16000, runtime: 7200}
```

In local mode, a job never uses more threads than `-t`, and Snakemake does not limit memory. Make sure the machine has enough RAM for the heavy steps (see [Hardware Requirements](#hardware-requirements)).

### FCS-GX Database

Decontamination runs only when `gxdb` is set. The FCS-GX database takes **~470 GB of disk**, and FCS-GX needs **~512 GB of RAM** to run at full speed, so store it on fast local storage when possible.

Download it with the `fcs.py` wrapper shipped with Pipeasm. First build the FCS-GX image (Pipeasm also builds it in `resources/` when `gxdb` is set, so you can reuse it later):

```bash
singularity build resources/fcs-gx.sif docker://ncbi/fcs-gx:latest

python3 workflow/scripts/fcs.py --image resources/fcs-gx.sif db get \
    --mft "https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest" \
    --dir /absolute/path/to/gxdb
```

Then set `gxdb: "/absolute/path/to/gxdb"` in `config/config.yaml`. See the [FCS-GX wiki](https://github.com/ncbi/fcs/wiki/FCS-GX) for details.

## Running Pipeasm

### Command-line Options

```bash
bash Pipeasm.sh -d </path/to/work/dir> [-c <config.yaml>] [-s <Snakefile>] (-t <threads> | -slurm -partition <name> -j <jobs>) [step] [-np] [-unlock]
```

| Flag | Required | Description |
| --- | --- | --- |
| `-d <path>` | Yes | Working directory (the cloned repository). |
| `-t <int>` | Yes, in local mode | Number of threads (cores) available to Snakemake. Not allowed with `-slurm`. |
| `-j <int>` | Yes, with `-slurm` | Maximum number of jobs submitted to the SLURM queue at the same time. |
| `-slurm` | No | Submit each rule as a SLURM job, using `config/slurm_params.yaml` and `config/resources/{genome_size}.yaml`. |
| `-partition <name>` | Yes, with `-slurm` | SLURM partition (queue) for the jobs. |
| `-c <path>` | No | Alternative config file (default: `config/config.yaml`). |
| `-s <path>` | No | Alternative Snakefile (default: `workflow/Snakefile`). |
| `-np` | No | Dry run: build the job graph without running any rule. |
| `-unlock` | No | Unlock the working directory after an interrupted run. |

**Steps.** By default, all steps run. To run only some of them, add one or more of:

| Flag | Step |
| --- | --- |
| `-trimming_qc` | Trimming and quality control |
| `-kmer_eval` | K-mer profiling statistics and plots |
| `-assembly` | Assembly, statistics and decontamination |
| `-scaffolding` | YaHS scaffolding and Hi-C contact maps |

Steps depend on the outputs of the previous ones. Run them in order, e.g. `-trimming_qc` before `-kmer_eval`.

### Local Mode

```bash
bash Pipeasm.sh -d /path/to/pipeasm -t 64
```

### SLURM Mode

```bash
bash Pipeasm.sh -d /path/to/pipeasm -slurm -partition <partition> -j 20
```

A submission script template is provided in [`pipeasm.slurm`](pipeasm.slurm). The Snakemake controller job itself needs few resources (e.g., 4 CPUs and 2 GB of RAM) but must run for the whole pipeline.

### Temporary Files

Pipeasm creates `singularity/` (Singularity cache) and `tmp/` (temporary files) in the working directory. Both can grow large, so the working directory should be on a file system with enough free space.

## Outputs

All results are written to `results/`. The main files are also collected into `workflow/report/`.

### Report (`workflow/report/`)

- **Summary**
  - `{sample}.summary.txt`: a plain-text summary of the main metrics for each assembly.
- **Assembly metrics**
  - `{sample}.gfastats.csv`: contig and scaffold counts, N50, L50, etc.
  - `{sample}.merqury.csv`: k-mer QV and completeness.
  - `{sample}.compleasm.csv`: BUSCO gene completeness.
  - `{sample}.genomescope2_model.csv` and `{sample}.genomescope2_summary.csv`: genome size, heterozygosity and repeat content estimates.
  - `{sample}.fcsadaptor.csv` and `{sample}.fcsgx.csv`: adapters and contaminants found.
- **Read and mapping metrics**
  - `{sample}.flagstats.csv`: Hi-C mapping statistics.
  - `{sample}.hifi_nanoplot.csv` and `{sample}.ont_nanoplot.csv`: read length and quality metrics.
  - `{sample}.fastqc_basic_stats.csv` and `{sample}.fastqc_overrepresented_sequences.csv`: FastQC summaries.
- **`images/`**
  - Hi-C contact maps (`*FullMap.png`), before and after scaffolding.
  - MerquryFK spectra plots (`*.spectra-cn.*.png`, `*.spectra-asm.*.png`).
  - Snail plots (`*_snail.png`).
  - KatGC and KatComp plots.
  - Read length and quality plots (`HiFi_*`, `ONT_*`).
  - Smudgeplots (`{sample}_smudgeplot*.png`).
- **`files/`**
  - FASTA files from Hifiasm and YaHS.
  - Mitogenome files (when `run_mitohifi: 'Yes'`).
  - FastQC HTML reports (`{sample}.trimmed_fastqc.html`, `{sample}_R1.trimmed_paired_fastqc.html`, `{sample}_R2.trimmed_paired_fastqc.html`).
  - Detailed statistics, e.g. `{sample}.completeness.stats` and `{sample}.contigs_stats.tsv`.

### Benchmark Report

At the end of every run (but not on dry runs), Pipeasm writes an HTML report of runtime, CPU and memory usage to `workflow/report/benchmark.html`. The report is also written when some jobs fail, so it always covers the jobs finished so far. It is built from the Snakemake benchmark files in `benchmarks/`, so it only shows the rules that actually ran for your inputs and options. It includes:

- Summary tiles: cumulative runtime, CPU time, peak memory, longest job and data read/written.
- The run configuration (inputs and options that decide which rules run).
- Runtime by stage (Trimming and QC, K-mer profiling, Assembly, Assembly statistics, Decontamination, Hi-C mapping and scaffolding, Scaffolding statistics).
- Resources by rule: runtime, CPU time, average cores used vs. the threads set in `config/threads.yaml`, and peak memory vs. the memory requested in `config/resources/{genome_size}.yaml`. Rules using 90% or more of the requested memory are highlighted. Use this table to tune the SLURM resources for your next genomes.
- The rules with the highest peak memory.

If the report cannot be generated, the run is not affected and the reason is written to `logs/benchmark_report.log`. To rebuild it manually, run from the Pipeasm working directory:

```bash
python workflow/scripts/benchmark_report.py -c config/config.yaml
```

See `--help` for the other options (e.g., `--all-samples` when the `benchmarks/` directory holds runs of more than one sample).

### Results Directory Structure

Only the directories that apply to your inputs and settings are created.

```
results/
├── Trimming_QC
│   ├── HiFi
│   ├── HiC
│   ├── ONT
│   ├── FastK
│   └── QC
│       ├── HiFi_FastQC
│       ├── HiFi_NanoPlot
│       ├── HiC_FastQC
│       └── ONT_NanoPlot
├── Assembly
│   ├── Contigging
│   │   ├── Phased_Asm
│   │   └── Solo_Asm
│   ├── Genome_Stats
│   │   ├── Compleasm
│   │   │   ├── Phased_Asm_Hap1
│   │   │   ├── Phased_Asm_Hap2
│   │   │   ├── Solo_Asm_Alt
│   │   │   └── Solo_Asm_Primary
│   │   ├── GFAstats
│   │   ├── GeneScopeFK
│   │   ├── KatComp
│   │   ├── KatGC
│   │   ├── MerquryFK
│   │   │   ├── Phased_Asm
│   │   │   └── Solo_Asm
│   │   ├── Smudgeplot
│   │   └── SnailPlot
│   │       ├── Phased_Asm
│   │       │   ├── Phased-Hap1
│   │       │   └── Phased-Hap2
│   │       └── Solo_Asm
│   │           ├── 00-Solo-Hap1
│   │           └── 01-Solo-Hap2
│   └── Mitogenome
│       ├── Phased_Asm
│       └── Solo_Asm
├── Decontamination
│   ├── Contaminants
│   │   ├── Phased_Asm_Hap1
│   │   ├── Phased_Asm_Hap2
│   │   ├── Solo_Asm_Alt
│   │   └── Solo_Asm_Primary
│   └── FCS-Adaptor
│       ├── Phased_Asm_Hap1
│       ├── Phased_Asm_Hap2
│       ├── Solo_Asm_Alt
│       └── Solo_Asm_Primary
└── Scaffolding
    ├── Initial_Contacts
    │   └── {Hap1,Hap2,Primary,Alternative}
    ├── YAHS_Scaffolding
    │   └── {Hap1,Hap2,Primary,Alternative}
    ├── Final_Contacts
    │   └── {Hap1,Hap2,Primary,Alternative}
    └── Scaffolding_stats
        ├── Compleasm
        │   └── {Hap1,Hap2,Primary,Alternative}
        ├── GFAstats
        ├── MerquryFK
        └── SnailPlot
            └── {Hap1,Hap2,Primary,Alternative}
```

## Test Dataset

We provide a test dataset, a ready-to-use configuration file and all expected output files on [Zenodo](https://zenodo.org/records/17243106?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImVjYjM1NWJjLTg1YzYtNDljOC1hMTYxLWFhNDEwNzM5ZDkyZSIsImRhdGEiOnt9LCJyYW5kb20iOiJhMzk0NDQ1NjdhYWY2MDBhZjM3ODExNjllODc0OWQ0OSJ9.b65MGFEuokHENviKgaDa9ITTpS9dOPRyCHzDSUAiikQDkW0CQDI5nUyAbWaHntvjqqNU847jdc9B5_7uEbSBcw), for *Panthera onca* (mPanOnc):

1. Example HiFi and Hi-C read files.
2. A pre-configured `config.yaml` for the test dataset.
3. The complete outputs and reports of a successful run, to compare with your results.

Use it to check that Snakemake, Singularity and all containers work before running your own data.

## Troubleshooting

- **`LockException` / "Directory cannot be locked"**: a previous run was interrupted. Run the same command with `-unlock`, then run it again without it.
- **Incomplete outputs after a crash**: just rerun the same command. The profiles enable `rerun-incomplete` and `keep-going`, so Snakemake resumes from where it stopped.
- **Where to look when a rule fails**: each rule writes its log to `logs/` and its resource usage to `benchmarks/` in the working directory. Snakemake also prints the failed rule's log at the end of the run, and `workflow/report/benchmark.html` shows the runtime and memory of the jobs that finished.
- **"No space left on device"**: `singularity/` and `tmp/` in the working directory and `sif_dir` can grow large. Move the working directory to a larger file system or clean old runs. Old image versions in `sif_dir` that are no longer listed in `workflow/containers.yaml` can be deleted.
- **`ERROR: Invalid genome_size`**: `genome_size` in `config.yaml` must be `small`, `medium` or `large`, even in local mode.
- **`-t is not allowed with -slurm`**: in SLURM mode, use `-j` for the number of queued jobs and `-partition` for the queue.
- **`ERROR: Found N problem(s) in ...config.yaml`**: the configuration check found invalid fields. Fix each item listed and run again.
- **A SLURM job runs out of memory or time**: raise `mem_mb` or `runtime` for that rule in `config/resources/{genome_size}.yaml` (see [Threads and Memory](#threads-and-memory)) and run the same command again.

## How to Cite

> Bruno Marques Silva, Fernanda de Jesus Trindade, Lucas Eduardo Costa Canesin, Giordano Souza, Alexandre Aleixo, Gisele Nunes, Renato Renison Moreira-Oliveira. Pipeasm: a tool for automated large chromosome-scale genome assembly and evaluation. *Bioinformatics Advances*, Volume 6, Issue 1, 2026, vbaf326. https://doi.org/10.1093/bioadv/vbaf326

```bibtex
@article{silva2026pipeasm,
  author  = {Silva, Bruno Marques and Trindade, Fernanda de Jesus and Canesin, Lucas Eduardo Costa and Souza, Giordano and Aleixo, Alexandre and Nunes, Gisele and Moreira-Oliveira, Renato Renison},
  title   = {Pipeasm: a tool for automated large chromosome-scale genome assembly and evaluation},
  journal = {Bioinformatics Advances},
  volume  = {6},
  number  = {1},
  pages   = {vbaf326},
  year    = {2026},
  doi     = {10.1093/bioadv/vbaf326}
}
```

## Changelog

See [CHANGELOG.md](CHANGELOG.md).
