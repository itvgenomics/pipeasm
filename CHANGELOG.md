# Changelog

## Pipeasm 1.1.3

### New Features
- Added a benchmark report (`workflow/report/benchmark.html`), written at the end of every run by `workflow/scripts/benchmark_report.py`. It reads whatever benchmark files exist, so it adapts to the rules that ran, and shows runtime, CPU and memory by stage and by rule, compared with the threads and SLURM memory set in `config/`.
- Added `workflow/scripts/check_config.py`, which validates the config file before every run and lists all problems at once (missing required fields, files not found, relative paths, invalid "Yes"/"No" values, incomplete Hi-C pairs, missing MitoHiFi fields, mismatched k-mer sizes).

### Changes
- Split the Snakemake profile templates so each value is set in one place:
  - `config/threads.yaml`: threads per rule, shared by local and SLURM modes.
  - `config/slurm_params.yaml`: general SLURM settings and default resources (2 GB, 72 h, partition, optional account).
  - `config/resources/{small,medium,large}.yaml`: memory per rule for each `genome_size`.
  - Removed `config/{small,medium,large}_genome_params.yaml`. Only rules that differ from the defaults are listed (about 900 to 140 lines per file).
- `Pipeasm.sh` no longer modifies the config file: the selected steps are passed to Snakemake with `--config`.
- Reorganized `config/config.yaml` into sections, with accepted values, examples and links for every field, and corrected descriptions of the advanced tool parameters.
- Set `max-jobs-per-second` to 10 in SLURM mode (it was equal to `-j`).

### Improvements
- Added a `--memory` limit to FastQC for HiFi reads.
- Rewrote the README: quick start, software versions, hardware requirements, supported input combinations, FCS-GX database setup and troubleshooting. Moved the changelog to `CHANGELOG.md`.

### Bug Fixes
- `phased_run_mitohifi` and `solo_run_mitohifi` wrote the same benchmark file (`{sample}.run_mitohifi.txt`); they now write `{sample}.phased_run_mitohifi.txt` and `{sample}.solo_run_mitohifi.txt`.
- Fixed the SLURM resources of `scaffolding_merqury_prim` (listed as `scaffolding_merqury_alt`), which ran with 1 thread and Snakemake's default memory when `diff_species_hic` was "Yes".
- Removed stale or misspelled rule names from the thread settings (`remove_buscodb`, `caffolding_gfastats_stats_alt`, `scaffolding_edit_busco_table_alt`).
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
