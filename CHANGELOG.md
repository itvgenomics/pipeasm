# Changelog

## Pipeasm 1.1.3

### Improvements
- Added a `--memory` limit to FastQC for HiFi reads.
- Rewrote the README: quick start, software versions, hardware requirements, supported input combinations, FCS-GX database setup and troubleshooting. Moved the changelog to `CHANGELOG.md`.

### Bug Fixes
- Fixed the `{ WORKDIR }` placeholder in `config/local_params.yaml`, which was not replaced in local mode, leaving the Singularity bind path as the literal `{WORKDIR}`.
- Replaced hard-coded SLURM partitions in `config/{small,medium,large}_genome_params.yaml` with the `{ PARTITION }` placeholder set by `-partition`.
- Removed leftover `run_*` flags from `config/config.yaml`.
- Fixed the `genome_size` description.

## Pipeasm 1.1.2

### New Features
- Added small, medium and large genome SLURM configuration files.
- Added the `genome_size` option to `config.yaml` to select the genome-size-specific SLURM configuration.

### Improvements
- Renamed `{sample}.abstract.txt` to `{sample}.summary.txt` in `workflow/report`.

### Bug Fixes
- Fixed parsing of Merqury QV values in the summary.
- Fixed the Cutadapt comment in the configuration file.
- Fixed the Merqury parser.

## Pipeasm 1.1.1

### Improvements
- Added mitogenome files to `workflow/report`.
- Updated `Compleasm` to fix a filename parsing failure.

## Pipeasm 1.1.0

### New Features
- Replaced `meryl` as the k-mer counter with:
  - `FastK`
  - `MerquryFK`
  - `GeneScopeFK`
- Added support for Hi-C data from a different species:
  - Performs a primary-only assembly with the `--primary` flag.
  - Applies `YaHS` scaffolding to the resulting contigs.
- Added the `-j` flag for SLURM mode to control the number of jobs submitted to the queue.

### Improvements
- Updated `Compleasm` to version 0.2.7.
- Reduced RAM usage in assembly rules from 700 GB to 500 GB.
- Marked raw reads as temporary files:
  - Raw reads are now automatically deleted after pipeline execution.
  - Only trimmed reads are kept.

### Changes
- Made `MitoHiFi` optional:
  - Controlled by the `run_mitohifi` option in the configuration file.

## Pipeasm 1.0.1

### New Features
- Added SLURM support (`-slurm`). Resources are now set in `config/{small,medium,large}_genome_params.yaml` (see 1.1.2).

### Improvements
- Improved handling of `MitoHiFi` runs that fail to assemble the mitogenome.
- Added a Singularity temporary directory (`tmpdir`) to the `.sh` scripts.
- Included YaHS FASTA outputs in the reports.
- Updated `Hifiasm` to version 0.25.0.

### Bug Fixes
- Fixed the abstract text formatting.
- Fixed file checking and the `fcs-adaptor` Singularity image.
- Fixed the `fcs` stderr output.
