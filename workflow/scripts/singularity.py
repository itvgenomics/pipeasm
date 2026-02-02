import os
import subprocess
import argparse

parser = argparse.ArgumentParser(description="Process some Docker images.")
parser.add_argument(
    "--config", type=str, required=True, help="Path to the configuration YAML file"
)

args = parser.parse_args()


def read_yaml(file_path):
    sif_dir = None
    with open(file_path, "r") as file:
        for line in file:
            # Strip whitespace and check for the sif_dir key
            line = line.strip()
            if line.startswith("sif_dir:"):
                # Extract the value after the colon and strip whitespace
                sif_dir = line.split(":", 1)[1].strip()
                break
    return sif_dir


docker_images = [
    "itvdsbioinfo/hic_mapping:1.0",
    "itvdsbioinfo/yahs:1.2a.2",
    "staphb/fastp:0.23.4",
    "staphb/fastqc:0.12.1",
    "itvdsbioinfo/fkutils:1.2",
    "pipecraft/cutadapt:4.4",
    "itvdsbioinfo/hifiasm:0.25.0",
    "nanoporetech/dorado:shae36d1b49fe470a60e006afad90bedd2fc2774a89",
    "staphb/nanoplot:1.41.6",
    "genomehubs/blobtoolkit:4.3.5",
    "huangnengcsu/compleasm:v0.2.7",
    "ghcr.io/marcelauliano/mitohifi:master",
    "staphb/gfastats:1.3.6",
]

sif_dir = read_yaml(args.config)
sif_dir = sif_dir.replace('"', "").replace("'", "")

os.makedirs(os.path.abspath(sif_dir), exist_ok=True)

for image in docker_images:
    print(f"INFO: Fetching {image} to {sif_dir}.")

    # Construct the singularity pull command with output directory
    output_file = os.path.join(
        os.path.abspath(sif_dir), f"{image.split('/')[-1].split(':')[0]}.sif"
    )

    # Check if the SIF file already exists
    if os.path.exists(output_file):
        print(f"INFO: {output_file} already exists. Skipping pull for {image}.")
        continue

    command = f"singularity pull {output_file} docker://{image}"

    # Execute the command
    subprocess.run(command, shell=True, check=True)
    print(f"INFO: Successfully pulled {image} to {output_file}")
