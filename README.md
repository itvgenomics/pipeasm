

# Pipeasm: a tool for automated large chromosome-scale genome assembly and evaluation

Pipeasm is the streamlining of the VGP and DToL assembly pipelines. The pipeline is written in Snakemake, a powerful workflow management system for creating and executing data analysis pipelines, and both container daemons Docker and Singularity. This document will guide you through the installation process of these software on Linux.

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
| solo_asm | 'Yes' | Run both solo and phased contig assembly. Fill with "Yes" or "No" | Required |
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
  - **-t** {int} (Required) = Number of threads to use

  - Only the **-d** and **-t** flags are required if you want to run Pipeasm with default parameters. If you do not want to run the decontamination step, set **gxdb**, in the **config.yaml** file, **blank**.
  - You can perform a dry-run (build only the DAG, no rule will be run) with -np and unlock the directory, if needed, with --unlock

- You can choose a Pipeasm step with:
  -  `--trimming_qc` (for Trimming and Quality Control);
  -  `--kmer_eval` (for k-mer Evaluation stats/plots);
  -  `--assembly` (for all Assembly and Decontamination steps);
  -  `--scaffolding` (to run YAHS auto-scaffolding and create the Hi-C Maps);

- Just make sure to run previous steps before running the next one.
  - Run
  `bash Pipeasm.sh -d /path/to/work/dir -t 64 --trimming_qc`
  before
  `bash Pipeasm.sh -d /path/to/work/dir -t 64 --kmer_eval`


- Pipeasm will create a temporary directory at you working dir to store all singularity images and Snakemake temp files.

### 5. Monitor Progress and View Results
- Snakemake will display progress information and any errors encountered during the execution.
- Once the pipeline completes successfully, you can find the output files in the `results` output directory as specified.
- Pipeasm provides a extensive report data located at `workflow/report` directory.

### Additional Notes
- Make sure to review and modify the Snakefile if necessary to fit your specific pipeline requirements.
- Make sure conda is installed with singularity and snakemake in the active environment.
- We are providing a [PBS](https://www.openpbs.org/) example file to run Pipeasm in a cloud/HPC environment. See `pipeasm.pbs`.

### Example Repository Structure
```
├── config
│   └── config.yaml
├── resources
│   ├── pacbio_adapters.txt
│   └── TruSeq3-PE.fa
├── workflow
│   ├── scripts
│   ├── rules
│   ├── report
│   ├── singularity
│   └── snakefile
├── results
│   ├── Trimming_QC
│   ├── Assembly
│   ├── Decontamination
│   └── Scaffolding
├── README.md
├── pipeasm.pbs
├── benchmarks
└── logs
```
