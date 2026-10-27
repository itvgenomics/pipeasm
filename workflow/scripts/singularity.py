"""Pull the Singularity images listed in workflow/containers.yaml into sif_dir.

Each image is saved as <name>_<tag>.sif (see sif_name), the same path the
rules get from sif("<name>") in workflow/rules/common.smk. Images that already
exist are skipped. Each pull goes to a temporary file that is renamed only when
the pull succeeds, so an interrupted pull never leaves a partial .sif behind.
"""

import argparse
import os
import re
import subprocess
import sys

import yaml

parser = argparse.ArgumentParser(description="Pull the Pipeasm Singularity images.")
parser.add_argument(
    "--config", type=str, required=True, help="Path to the configuration YAML file"
)
parser.add_argument(
    "--containers",
    type=str,
    default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "containers.yaml"),
    help="Container registry (default: workflow/containers.yaml)",
)


def sif_name(name, image):
    """<name>_<tag>.sif; must match sif() in workflow/rules/common.smk."""
    tag = image.rsplit(":", 1)[1] if ":" in image.split("/")[-1] else "latest"
    return f"{name}_{re.sub(r'[^A-Za-z0-9._-]', '_', tag)}.sif"


def main():
    args = parser.parse_args()

    with open(args.config) as handle:
        sif_dir = str((yaml.safe_load(handle) or {}).get("sif_dir") or "").strip()
    if not sif_dir:
        print("ERROR: 'sif_dir' is not set in the configuration file.")
        sys.exit(1)
    with open(args.containers) as handle:
        containers = yaml.safe_load(handle) or {}

    sif_dir = os.path.abspath(sif_dir)
    os.makedirs(sif_dir, exist_ok=True)

    failed = []
    for name, image in containers.items():
        output_file = os.path.join(sif_dir, sif_name(name, image))

        if os.path.exists(output_file):
            print(f"INFO: {output_file} already exists. Skipping pull for {image}.")
            continue

        print(f"INFO: Fetching {image} to {output_file}.")
        partial_file = output_file + ".partial"
        if os.path.exists(partial_file):
            os.remove(partial_file)

        result = subprocess.run(["singularity", "pull", partial_file, f"docker://{image}"])
        if result.returncode != 0 or not os.path.exists(partial_file):
            if os.path.exists(partial_file):
                os.remove(partial_file)
            failed.append((image, output_file))
            continue

        os.replace(partial_file, output_file)
        print(f"INFO: Successfully pulled {image} to {output_file}")

    if failed:
        print(f"ERROR: Could not pull {len(failed)} image(s). If this machine has no internet access, "
              "pull them on another machine and copy them to sif_dir:")
        for image, output_file in failed:
            print(f"  singularity pull {os.path.basename(output_file)} docker://{image}   # -> {output_file}")
        sys.exit(1)


if __name__ == "__main__":
    main()
