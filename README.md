# Pipeasm: a tool for automated large chromosome-scale genome assembly and evaluation

<p align="center">
  <img src="https://github.com/itvgenomics/pipeasm/blob/main/Pipeasm.png" alt="Pipeasm Logo"/>
</p>

Pipeasm is the streamlining of the VGP and DToL assembly pipelines. The pipeline is written in Snakemake, a powerful workflow management system for creating and executing data analysis pipelines, and both container daemons Docker and Singularity. This document will guide you through the installation process of these software on Linux.

## How Pipeasm works?

Pipeasm is a Snakemake-based workflow that integrates a suite of tools to process and assemble genomic data, producing comprehensive reports and statistics. The pipeline is organized into six main steps, as shown in the schematic diagram.

<p align="center">
  <img src="https://github.com/itvgenomics/pipeasm/blob/main/workflow.png" alt="Workflow"/>
</p>

1. **Trimming and Quality Control (QC)**
- Pipeasm can use **PacBio high-fidelity (HiFi-CCS) long reads**, **Oxford Nanopore Long-Reads (ONT)**, and **chromatin-contact (Hi-C) Illumina short reads**.
- **Trimming**: It uses `Cutadapt` for HiFi reads, `Fastp` for Hi-C reads and `Dorado` for ONT reads, to remove adapters and primers.
- **Quality Assessment**: `FastQC` and `Nanoplot` perform quality control on the reads.

2. **K-mer Profiling**
- After trimming, the pipeline performs a series of analyses on the HiFi long reads to gather detailed genomic statistics.
- **K-mer Counting**: `Meryl` counts and evaluates K-mers.
- **Genome Statistics**: `GenomeScope2` estimates genome size, heterozygosity, and repeat content.
- **Ploidy Identification**: `SmudgePlot` identifies ploidy levels through K-mer distribution analysis.
- **Graphical Representation**: `KAT-GCP` provides graphical representations of K-mer spectra to help identify genomic features and contaminants.

3. **Assembly**
- This is the core of the pipeline where the genome is assembled.
- **Hifiasm** produces haplotype-resolved genomes, which is crucial for creating accurate and contiguous assemblies.
- **Mitogenome Assembly**: `MitoHiFi` recovers the mitochondrial genome.

4. **Assembly Statistics**
- After assembly, the pipeline generates a fasta file and collects detailed statistics on the assembled genome.
- **Summaries and Statistics**: `GFAstats` generates the fasta file and provides statistics like contig lengths, N50 values, and scaffold counts.
- **Completeness Evaluation**: `Compleasm` assesses gene-space completeness by comparing the assembled genome against BUSCO databases.
- **Quality Assessment**: `Merqury` counts assembled K-mers and compares them to the Meryl database.
- **Visualization**: `Snailplot`, which is a part of BlobToolkit , provides a visual overview of the assembled genome.

5. **Decontamination**
- This step eliminates remaining contaminants from the assembled genome.
- `FCS Adaptor` removes adapters, and `FCS-GX` targets other contaminants.

6. **Hi-C Mapping and Scaffolding**
- This final stage links contigs into longer, chromosome-scale scaffolds.
- **Mapping**: `BWA-MEM2` is used for genome indexing and mapping Hi-C reads.
- **Sorting and Filtering**: Aligned reads are sorted with `SAMtools`, and low-quality and duplicated reads are filtered and merged using `Bellerophon`.
- **Scaffolding**: `YAHS` links the contigs into longer scaffolds, enhancing the assembly.
- **Visualization**: `Pretext` tools (Map and Snapshot) convert the mapping results into visual representations of Hi-C contact maps

## How to Install Pipeasm Enviroment on Linux

### Prerequisites
Before installing the required software, make sure you have the following:
- A Linux-based operating system (e.g., Ubuntu, CentOS, Fedora)
- Python (version 3.5 or later) installed on your system

### Snakemake/Singularity/Docker Instalation
1. - **Install Conda or Mamba**: The pipeline requires a conda/mamba installation with snakemake and singularity installed in the base environment. You can find the installation steps [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). Then activate the base environment and install:
  ```
  conda install -c bioconda -c conda-forge singularity=3.8.5 snakemake=8.5.1
  ```
  or
  ```
  mamba install -c bioconda -c conda-forge singularity=3.8.5 snakemake=8.5.1
  ```

  - You can also create a conda new enviroment:
  ```
  conda create -n snakemake -c bioconda -c conda-forge singularity=3.8.5 snakemake=8.5.1

  conda activate snakemake
  ```

**Make sure you are installing Singularity, Snakemake and Docker compatible versions.**

2. - **Install Docker**: Docker installation can be found [here](https://docs.docker.com/engine/install/). If you instal the above recommended versions, use the Docker 24.0.7 version.

### Additional Resources
- [Snakemake Documentation](https://snakemake.readthedocs.io/en/stable/)
- [Snakemake GitHub Repository](https://github.com/snakemake/snakemake)

## How to run the Pipeasm Workflow

### 1. Clone the Snakemake Repository from GitHub
- Open your terminal.
- Navigate to the directory where you want to clone the repository.
- Run the following command to clone the repository:
  ```bash
  git clone https://github.com/itvgenomics/pipeasm.git
  ```

### 2. Navigate to the Cloned Repository
- Change directory to the cloned repository:
  ```bash
  cd <repository_name>
  ```
  Replace `<repository_name>` with the name of the cloned repository directory.

### 4. Check Raw Data
- Check if your HiFi/Hi-C/ONT reads are available and joint in its respective files and compressed with gzip:
  - /path/to/hifi_reads/hifi_reads.fasq.gz
  - /path/to/hic_reads/r1_reads.fasq.gz
  - /path/to/hic_reads/r2_reads.fasq.gz
  - /path/to/ont_reads/ont_reads.fasq.gz

### 3. Configure the Pipeline
- Edit the configuration file according to your requirements at the `config` directory.
  - `config.yaml`: General configuration settings.
    ```bash
    vim config/config.yaml
    ```

The fields to be edited are the following:

| Field    | Example | Description | Field Requirement |
| -------- | ------- | ------- | ------- |
| species | 'Hypocnemis striata' | Name of the species to be assembled. This will be used to search a complete mitochondrial genome to be used as reference in MitoHiFi. | Required |
| sample | 'bHypStr1.1' | This is the DToL_ID naming the assembly for each species. You can check the documentation and get the corresponding ID for your species [here](https://id.tol.sanger.ac.uk) | Required |
| hifi_reads | '/absolute/path/to/raw_data/hifi/bHypStr1.fastq.gz' | Full path to your HiFI raw reads. Make sure to get all your reads into a single `fastq.gz` file | Required |
| hic_r1 | '/absolute/path/to/raw_data/hic/bHypStr1.R1.fastq.gz' | Full path to your HiC forward raw reads. Make sure to get all your forward reads into a single `fastq.gz` file | Optional |
| hic_r2 | '/absolute/path/to/raw_data/hic/bHypStr1.R2.fastq.gz' | Full path to your HiC reverse raw reads. Make sure to get all your reverse reads into a single `fastq.gz` file | Optional |
| ont_reads | '/absolute/path/to/raw_data/hic/bHypStr1.ONT.fastq.gz' | Full path to your ONT demultiplexed reads. | Optional |
| geneticcode | 2 | Genetic Code to be used at mitochondrial genome annotation. This can be left as it is. Only basal eukaryotes or special cases will require to change this parameter | Required |
| taxid | 1202453 | You can check the Tax_ID for your species in NCBI's Taxonomy database. If your species does not have a TaxID, you can set a similar organism. This will be used at the decontamination step. | Required |
| buscodb | 'aves' | Pick the closest taxonomic group available for your species. You can check the available datasets [here](https://busco-data.ezlab.org/v5/data/lineages/) | Required |
| solo_asm | 'Yes' | If "Yes", run both solo and phased contig assembly. If "No", run only phased assembly. | Required |
| threads | 32 | We recommend 32 to 64 threads (which may require ~150Gb of RAM for a 1Gb genome assembly, and up to 500Gb for decontamination and kmer analysis) | Required |
| sif_dir | '/absolute/path/to/singularity/images/' | Directory to build all singularity image files used in the pipeline. If the path already contains the images, they will not be pulled. | Required |
| gxdb | '/absolute/path/to/fcs-gx_DB/gxdb' | FCS-GX database to perform the contamination step. To download, see this [guide](https://github.com/ncbi/fcs/wiki/FCS-GX-quickstart#download-the-fcs-gx-database). Leave it blank if you dont want to run this step. | Optional |

- We have set default parameters to some software threads and flags. You can change them at the `config/config.yaml` file.

### 4. Run the Pipeline
- Navigate back to the root directory of the cloned repository.
- Run the pipeline with default parameters:
  ``` bash
  bash Pipeasm.sh -d /path/to/work/dir -t <n_threads>

 - Run with all available parameters:
  ``` bash
  bash Pipeasm.sh -d /path/to/work/dir -c </path/to/config.yaml> -s </path/to/snakefile> -t <n_threads> <step>
  ```
- **Extra Flags**:
  - **-d** </path/to/work/dir> (Required) = Path to your working directory where all the workflow file are
  - **-c** </path/to/config.yaml> (Optional) = Overwrite the default configuration file with all nedded parameters (config/config.yaml)
  - **-s** </path/to/snakefile> (Optional) = Overwrite the default snakefile path (workflow/snakefile)
  - **-t** {int} (Required) = Number of threads to use. If running in **SLURM** mode, this value determines how many jobs will be submitted to the queue.
  - **-slurm** (Optional) = Use the `config/slurm_params.yaml` file to run the workflow with SLURM job submission using Snakemake’s profile system. This enables use of SLURM-specific resource configuration, submission rules, and cluster-specific options. If you want to change any default SLURM settings, such as the partition: Edit `config/slurm_params.yaml` and set the appropriate value for the `slurm_partition` variable
  - **-partition** {string} (Required when **-slurm** is used) — Specifies the SLURM partition (queue) to which the jobs will be submitted.
  - **-np** (Optional) = Perform a dry run to see what jobs will be executed without actually running them.
  - **-unlock** (Optional) = Unlock the working directory if Snakemake has somehow locked it.

  - Only the **-d** and **-t** flags are required if you want to run Pipeasm with default parameters. If you do not want to run the decontamination step, set **gxdb**, in the **config.yaml** file, **blank**.
  - You can perform a dry-run (build only the DAG, no rule will be run) with `-np` and unlock the directory, if needed, with `-unlock`
  - If you want to change any default settings, such as the threads number and memory usage: Edit `config/local_params.yaml`. **DO NOT CHANGE THE {WORKDIR}/{THREADS}/{PARTITION} VARIABLES**. Change the `{PARTITION}` only if you want to set multiple slurm partitions.

- You can choose a Pipeasm step with:
  -  `-trimming_qc` (for Trimming and Quality Control);
  -  `-kmer_eval` (for k-mer Evaluation stats/plots);
  -  `-assembly` (for all Assembly and Decontamination steps);
  -  `-scaffolding` (to run YAHS auto-scaffolding and create the Hi-C Maps);

- Just make sure to run previous steps before running the next one.
  - Run
  `bash Pipeasm.sh -d /path/to/work/dir -t 64 -trimming_qc`
  before
  `bash Pipeasm.sh -d /path/to/work/dir -t 64 -kmer_eval`


- Pipeasm will create a temporary directory at you working dir to store all singularity images and Snakemake temp files.

### 5. Monitor Progress and View Results
- Snakemake will display progress information and any errors encountered during the execution.
- Once the pipeline completes successfully, you can find the output files in the `results` output directory as specified.
- Pipeasm provides a extensive report data located at `workflow/report` directory.

## How to Test the Pipeline

We provide a test dataset, a ready-to-use configuration file, and all expected output files to help you quickly test and validate your Pipeasm installation.
These resources are available on [Zenodo](https://zenodo.org/records/17243106?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImVjYjM1NWJjLTg1YzYtNDljOC1hMTYxLWFhNDEwNzM5ZDkyZSIsImRhdGEiOnt9LCJyYW5kb20iOiJhMzk0NDQ1NjdhYWY2MDBhZjM3ODExNjllODc0OWQ0OSJ9.b65MGFEuokHENviKgaDa9ITTpS9dOPRyCHzDSUAiikQDkW0CQDI5nUyAbWaHntvjqqNU847jdc9B5_7uEbSBcw).

This Zenodo repository contains, for the *Panthera onca* (mPanOnc):

1. Example HiFi and Hi-C read files for testing.
2. A pre-configured config.yaml file tailored for the test dataset.
3. Complete expected output files and reports from a successful run — useful for verifying your results match the reference outputs.

You can use this dataset to ensure that your installation of Snakemake, Singularity, and all dependent containers are working properly before running your own genomic data.


#### Summary Files

The main output files include:

- **Assembly Metrics**
  - `{sample}.gfastats.csv`: General assembly statistics including contig and scaffold counts, N50, and L50.
  - `{sample}.merqury.csv`: K-mer-based assembly quality metrics for evaluating base-level accuracy.
  - `{sample}.compleasm.csv`: Completeness assessment using conserved single-copy orthologs.
  - `{sample}.genomescope2_model.csv` and `{sample}.genomescope2_summary.csv`: Genome characterization and heterozygosity estimates.
  - `{sample}.fcsadaptor.csv` and `{sample}.fcsgx.csv`: Adapters and contamination statistics.

- **Metrics**
  - `{sample}.flagstats.csv`: Alignment statistics from mapped reads.
  - `{sample}.hifi_nanoplot.csv`: HiFi read quality metrics and length distributions.
  - `Files/{sample}.completeness.stats` and `Files/{sample}.contigs_stats.tsv`: Detailed statistics for contigs, scaffolds, and completeness.

#### Visualizations

The `Images` folder contains a range of figures to aid interpretation:

- **Contact Maps**:
  - `{sample}.final_contacts_Hap1FullMap.png` and `{sample}.final_contacts_Hap2FullMap.png`: Hi-C contact maps for haplotype assemblies.

- **Assembly Spectra**:
  - Merqury phased and solo assemblies visualized via spectra plots (`*.spectra-cn.st.png`, `*.spectra-asm.st.png`).

- **Scaffolding and Snail Plots**:
  - `*_snail.png` files provide circular visualizations of scaffolding quality and haplotype structure.

- **Read Quality and Length Distribution**:
  - HiFi read length and quality plots (`HiFi_*` files) including linear, log-transformed, and weighted histograms.

- **Smudgeplots**:
  - `{sample}_smudgeplot.pdf` and `{sample}_smudgeplot_log10.pdf` visualize k-mer multiplicity and genome ploidy/heterozygosity patterns.

#### HTML Reports

Interactive HTML reports are included in the `Files` folder for detailed examination of read quality:

- `{sample}.trimmed_fastqc.html`
- `{sample}_R1.trimmed_paired_fastqc.html`
- `{sample}_R2.trimmed_paired_fastqc.html`

These reports allow in-depth inspection of per-base quality, GC content, duplication rates, and other QC metrics.

### Additional Notes
- Make sure to review and modify the Snakefile if necessary to fit your specific pipeline requirements.
- Make sure conda is installed with singularity and snakemake in the active environment.
- We are providing a [PBS](https://www.openpbs.org/) example file to run Pipeasm in a cloud/HPC environment. See `pipeasm.pbs`.

### Example of the Results Structure
```
├── Assembly
│ ├── Contigging
│ │ ├── Phased_Asm
│ │ └── Solo_Asm
│ ├── Genome_Stats
│ │ ├── Compleasm
│ │ │ ├── Phased_Asm_Hap1
│ │ │ ├── Phased_Asm_Hap2
│ │ │ ├── Solo_Asm_Alt
│ │ │ └── Solo_Asm_Primary
│ │ ├── GFAstats
│ │ ├── HiFi_GenomeScope2
│ │ ├── KAT
│ │ ├── Merqury
│ │ │ ├── Phased_Asm
│ │ │ └── Solo_Asm
│ │ ├── Smudgeplot
│ │ └── SnailPlot
│ │     ├── Phased_Asm
│ │     │ ├── Phased-Hap1
│ │     │ └── Phased-Hap2
│ │     └── Solo_Asm
│ │         ├── 00-Solo-Hap1
│ │         └── 01-Solo-Hap2
│ └── Mitogenome
│     ├── Phased_Asm
│     └── Solo_Asm
├── Decontamination
│ ├── Contaminants
│ │ ├── Phased_Asm_Hap1
│ │ ├── Phased_Asm_Hap2
│ │ ├── Solo_Asm_Alt
│ │ └── Solo_Asm_Primary
│ └── FCS-Adaptor
│     ├── Phased_Asm_Hap1
│     ├── Phased_Asm_Hap2
│     ├── Solo_Asm_Alt
│     └── Solo_Asm_Primary
├── Scaffolding
│ ├── Final_Contacts
│ │ ├── Hap1
│ │ └── Hap2
│ ├── Initial_Contacts
│ │ ├── Hap1
│ │ └── Hap2
│ ├── Scaffolding_stats
│ │ ├── Compleasm
│ │ │ ├── Hap1
│ │ │ └── Hap2
│ │ ├── GFAstats
│ │ ├── Merqury
│ │ └── SnailPlot
│ │     ├── Hap1
│ │     └── Hap2
│ └── YAHS_Scaffolding
│     ├── Hap1
│     └── Hap2
└── Trimming_QC
    ├── HiC
    ├── HiFi
    ├── Meryl_DB
    └── QC
        ├── HiC_FastQC
        ├── HiFi_FastQC
        └── HiFi_NanoPlot
```

# Changelog

## Pipeasm 1.0.1 – Changelog

### New Features
- Added Slurm support:
  - To modify default parameters, edit `profiles/slurm/config.yaml`.
  - Set your partition using the `slurm_partition` variable.

### Improvements
- Enhanced handling when `mitoHiFi` fails to assemble the mitogenome.
- Added `Singularity` temporary directory (`tmpdir`) handling in the `.sh` scripts.
- Included YAHS FASTA outputs in the reports.
- Updated `hifiasm` to version 0.25.0.

### Bug Fixes
- Fixed abstract text formatting.
- Corrected file checking and `fcs-adaptor` Singularity image issues.
- Fixed `fcs` stderr output.
