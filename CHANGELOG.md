# Changelog

## Pipeasm 1.2.0

### Upgrading from 1.1.x
- Custom threads or memory set in `config/{small,medium,large}_genome_params.yaml` must be moved to `config/threads.yaml` and `config/resources/{small,medium,large}.yaml`.
- Images in `sif_dir` are now named `<name>_<tag>.sif`, so all images are pulled again on the first run. The old `.sif` files without a version (e.g., `hic_mapping.sif`) can be deleted.
- Write `"Yes"`/`"No"` values in quotes in `config/config.yaml`; the configuration check rejects unquoted values.

### New Features
- Added a benchmark report (`workflow/report/benchmark.html`), written at the end of every run by `workflow/scripts/benchmark_report.py`. It reads whatever benchmark files exist, so it adapts to the rules that ran, and shows runtime, CPU and memory by stage and by rule, compared with the threads and SLURM memory set in `config/`.
- Added `workflow/scripts/check_config.py`, which validates the config file before every run and lists all problems at once (missing required fields, files not found, relative paths, invalid "Yes"/"No" values, incomplete Hi-C pairs, missing MitoHiFi fields, mismatched k-mer sizes).

### Changes
- Split the Snakemake profile templates so each value is set in one place:
  - `config/threads.yaml`: threads per rule, shared by local and SLURM modes.
  - `config/slurm_params.yaml`: general SLURM settings and default resources.
  - `config/resources/{small,medium,large}.yaml`: memory per rule for each `genome_size`.
  - Removed `config/{small,medium,large}_genome_params.yaml`. Only rules that differ from the defaults are listed.
- `Pipeasm.sh` no longer modifies the config file: the selected steps are passed to Snakemake with `--config`.
- Reorganized `config/config.yaml` into sections, with accepted values, examples and links for every field, and corrected descriptions of the advanced tool parameters.
- Set `max-jobs-per-second` to 10 in SLURM mode (it was equal to `-j`).
- Updated `PretextMap` to version 0.2.4 and `PretextSnapshot` to version 0.1.1 in the new `itvdsbioinfo/hic_mapping:1.1` image. Removed the unused `PretextGraph` and `pretext-suite`.
- Hi-C contact maps are now built with `PretextMap --ultraRes`. The Pretext rules now request 50 GB of RAM and 4 threads.
- Container images are now listed in `workflow/containers.yaml` and saved in `sif_dir` as `<name>_<tag>.sif`. The rules use `sif("<name>")`, so updating a tag makes Pipeasm pull and use the new image, while older versions stay available.
- The `/tmp` of every container is now mounted on `$WORKDIR/tmp`:
  - Rules with containers: `-B { WORKDIR }/tmp:/tmp` in `singularity-args` (`config/local_params.yaml` and `config/slurm_params.yaml`).
  - FCS-adaptor and FCS-GX: new `--tmp-dir <dir>` option in `workflow/scripts/run_fcsadaptor.sh` and `workflow/scripts/fcs.py`, used with `--tmp-dir tmp` by the decontamination rules.

### Improvements
- Added the Hi-C mapping statistics (`{sample}.merged.bam.hic_stats` and `{sample}.merged.bam.flagstats`) of the initial and final contact maps (Hap1/Hap2, or Primary/Alternative when `diff_species_hic` is "Yes") to `workflow/report/files`, named `Initial_Contacts_<assembly>_*` and `Final_Contacts_<assembly>_*`.
- Added a `--memory` limit to FastQC for HiFi reads.
- Rewrote the README: quick start, software versions, hardware requirements, supported input combinations, FCS-GX database setup and troubleshooting. Moved the changelog to `CHANGELOG.md`.

### Bug Fixes
- An updated image was never used when a `.sif` with the same name already existed in `sif_dir`.
- An interrupted image pull left a partial `.sif` that was reused in later runs. Pulls now go to a temporary file that is renamed only when complete.
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
- Added SLURM support (`-slurm`). Resources are now set in `config/{small,medium,large}_genome_params.yaml` (see 1.1.2; replaced in 1.2.0).

### Improvements
- Improved handling of `MitoHiFi` runs that fail to assemble the mitogenome.
- Added a Singularity temporary directory (`tmpdir`) to the `.sh` scripts.
- Included YaHS FASTA outputs in the reports.
- Updated `Hifiasm` to version 0.25.0.

### Bug Fixes
- Fixed the abstract text formatting.
- Fixed file checking and the `fcs-adaptor` Singularity image.
- Fixed the `fcs` stderr output.
