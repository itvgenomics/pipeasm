"""Check the Pipeasm configuration file before running the workflow.

Usage: python check_config.py --config config/config.yaml

Prints every problem found and exits with status 1 if there is any error.
Warnings are printed but do not stop the run.
"""

import argparse
import os
import re
import sys

import yaml

YES_NO_FIELDS = ["solo_asm", "diff_species_hic", "run_mitohifi"]
GENOME_SIZES = ["small", "medium", "large"]
TOOL_SECTIONS = [
    "fastp",
    "cutadapt",
    "nanoplot",
    "fastk",
    "genescopeFK",
    "smudgeplot",
    "hifiasm",
    "gfastats",
    "bellerophon",
    "pretext",
    "yahs",
]


def is_empty(value):
    return value is None or str(value).strip() == ""


def as_text(value):
    return "" if value is None else str(value).strip()


def check_reads(config, field, errors, required=False):
    value = as_text(config.get(field))
    if not value:
        if required:
            errors.append(f"'{field}' is required: set the absolute path to the reads (.fastq.gz).")
        return
    if not os.path.isabs(value):
        errors.append(f"'{field}' must be an absolute path: '{value}'.")
    elif not os.path.isfile(value):
        errors.append(f"'{field}' file not found: '{value}'.")
    if not value.endswith((".fastq.gz", ".fq.gz")):
        errors.append(f"'{field}' must be a gzip-compressed FASTQ file (.fastq.gz): '{value}'.")


def check_config(config):
    errors, warnings = [], []

    # Sample name (must match the `sample` wildcard constraint in workflow/rules/common.smk)
    sample = as_text(config.get("sample"))
    if not sample:
        errors.append("'sample' is required.")
    elif not re.fullmatch(r"[A-Za-z0-9_]+", sample):
        errors.append(f"'sample' must contain only letters, numbers and underscores: '{sample}'.")

    # Reads
    check_reads(config, "hifi_reads", errors, required=True)
    check_reads(config, "hic_r1", errors)
    check_reads(config, "hic_r2", errors)
    check_reads(config, "ont_reads", errors)
    has_hic = not is_empty(config.get("hic_r1")) and not is_empty(config.get("hic_r2"))
    if is_empty(config.get("hic_r1")) != is_empty(config.get("hic_r2")):
        errors.append("'hic_r1' and 'hic_r2' must be set together (or both left empty).")

    # Yes/No options
    for field in YES_NO_FIELDS:
        value = as_text(config.get(field))
        if isinstance(config.get(field), bool):
            # YAML reads unquoted Yes/No as booleans
            errors.append(f"'{field}' must be quoted: \"Yes\" or \"No\".")
        elif value.lower() not in ("yes", "no"):
            errors.append(f"'{field}' must be \"Yes\" or \"No\", got: '{value}'.")

    genome_size = as_text(config.get("genome_size"))
    if genome_size not in GENOME_SIZES:
        errors.append(f"'genome_size' must be one of {GENOME_SIZES}, got: '{genome_size}'.")

    solo_asm = as_text(config.get("solo_asm")).lower()
    diff_species_hic = as_text(config.get("diff_species_hic")).lower()
    if diff_species_hic == "yes" and not has_hic:
        errors.append("'diff_species_hic' is \"Yes\" but no Hi-C reads were set in 'hic_r1'/'hic_r2'.")
    if not has_hic and solo_asm == "no":
        warnings.append(
            "No Hi-C reads and 'solo_asm' is \"No\": no assembly will be produced "
            "(only trimming/QC and k-mer profiling). Set 'solo_asm' to \"Yes\" to assemble the HiFi reads."
        )

    # Completeness and contaminants
    buscodb = as_text(config.get("buscodb"))
    if not buscodb:
        errors.append("'buscodb' is required (e.g., \"aves\").")
    elif buscodb.endswith("_odb10") or buscodb.endswith("_odb12"):
        errors.append(f"'buscodb' must not include the '_odb' suffix: '{buscodb}'.")

    taxid = as_text(config.get("taxid"))
    if not taxid:
        errors.append("'taxid' is required (NCBI Taxonomy ID, e.g., 1202453).")
    elif not taxid.isdigit():
        errors.append(f"'taxid' must be a number: '{taxid}'.")

    gxdb = as_text(config.get("gxdb"))
    if gxdb:
        if not os.path.isabs(gxdb):
            errors.append(f"'gxdb' must be an absolute path: '{gxdb}'.")
        elif not os.path.isdir(gxdb):
            errors.append(f"'gxdb' directory not found: '{gxdb}'.")

    # Mitogenome
    if as_text(config.get("run_mitohifi")).lower() == "yes":
        if is_empty(config.get("species")):
            errors.append("'species' is required when 'run_mitohifi' is \"Yes\".")
        geneticcode = as_text(config.get("geneticcode"))
        if not geneticcode:
            errors.append("'geneticcode' is required when 'run_mitohifi' is \"Yes\".")
        elif not geneticcode.isdigit():
            errors.append(f"'geneticcode' must be a number: '{geneticcode}'.")

    # Containers
    sif_dir = as_text(config.get("sif_dir"))
    if not sif_dir:
        errors.append("'sif_dir' is required (directory to store the Singularity images).")
    elif not os.path.isabs(sif_dir):
        errors.append(f"'sif_dir' must be an absolute path: '{sif_dir}'.")

    # Advanced tool parameters
    for section in TOOL_SECTIONS:
        if section not in config:
            errors.append(f"Tool parameter section '{section}' is missing. Copy it from the default config/config.yaml.")

    fastk_k = re.search(r"\d+", as_text((config.get("fastk") or {}).get("kmer")))
    genescope_k = re.search(r"\d+", as_text((config.get("genescopeFK") or {}).get("kmer")))
    if fastk_k and genescope_k and fastk_k.group() != genescope_k.group():
        errors.append(
            f"'fastk.kmer' (k={fastk_k.group()}) and 'genescopeFK.kmer' (k={genescope_k.group()}) must use the same k-mer size."
        )

    return errors, warnings


def main():
    parser = argparse.ArgumentParser(description="Check the Pipeasm configuration file.")
    parser.add_argument("--config", required=True, help="Path to the configuration YAML file")
    args = parser.parse_args()

    try:
        with open(args.config) as handle:
            config = yaml.safe_load(handle) or {}
    except (OSError, yaml.YAMLError) as error:
        print(f"ERROR: Could not read the configuration file '{args.config}':\n{error}")
        sys.exit(1)

    errors, warnings = check_config(config)

    for warning in warnings:
        print(f"WARNING: {warning}")
    if errors:
        print(f"ERROR: Found {len(errors)} problem(s) in '{args.config}':")
        for error in errors:
            print(f"  - {error}")
        sys.exit(1)

    print("INFO: Configuration file OK.")


if __name__ == "__main__":
    main()
