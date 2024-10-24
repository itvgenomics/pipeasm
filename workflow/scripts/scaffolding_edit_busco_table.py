import re

def parse_scaffolding_summary(sample, hap):

    if hap == 1:
        # Read the text file
        with open(f"results/Scaffolding/Scaffolding_stats/Compleasm/Hap1/{sample}.summary.txt", "r") as file:
            file_content = file.read()
    elif hap == 2:
        with open(f"results/Scaffolding/Scaffolding_stats/Compleasm/Hap2/{sample}.summary.txt", "r") as file:
            file_content = file.read()

    # Regular expression patterns
    lineage_pattern = r"## lineage:\s*(\w+)"
    n_number_pattern = r"N:(\d+)"

    # Find lineage
    lineage_match = re.search(lineage_pattern, file_content)
    if lineage_match:
        lineage = lineage_match.group(1)

    # Find N number
    n_number_match = re.search(n_number_pattern, file_content)
    if n_number_match:
        n_number = int(n_number_match.group(1))

    return lineage, n_number

def edit_scaffolding_table(buscodb, n_busco, sample, hap):
    line1 = "# BUSCO version is: 5.0.0"
    line2 = f"# The lineage dataset is: {buscodb} (Creation date: 2021-02-19, number of species: 4085, number of BUSCOs: {n_busco})"

    if hap == 1:
        # Open the file and read existing content
        with open(
            f"results/Scaffolding/Scaffolding_stats/Compleasm/Hap1/{buscodb}/full_table_busco_format.tsv",
            "r",
        ) as file:
            existing_content = file.readlines()

        # Prepend the new lines to the existing content
        updated_content = [line1 + "\n", line2 + "\n"] + existing_content

        # Write the updated content back to the file
        with open(
            f"results/Scaffolding/Scaffolding_stats/Compleasm/Hap1/{sample}.full_table_busco_format_edit.tsv",
            "w",
        ) as outfile:
            outfile.writelines(updated_content)

    elif hap == 2:
        # Open the file and read existing content
        with open(
            f"results/Scaffolding/Scaffolding_stats/Compleasm/Hap2/{buscodb}/full_table_busco_format.tsv",
            "r",
        ) as file:
            existing_content = file.readlines()

        # Prepend the new lines to the existing content
        updated_content = [line1 + "\n", line2 + "\n"] + existing_content

        # Write the updated content back to the file
        with open(
            f"results/Scaffolding/Scaffolding_stats/Compleasm/Hap2/{sample}.full_table_busco_format_edit.tsv",
            "w",
        ) as outfile:
            outfile.writelines(updated_content)

# Hap1
buscodb, n_busco = parse_scaffolding_summary(snakemake.config['sample'], 1)
edit_scaffolding_table(buscodb, n_busco, snakemake.config['sample'], 1)

# Hap2
buscodb, n_busco = parse_scaffolding_summary(snakemake.config['sample'], 2)
edit_scaffolding_table(buscodb, n_busco, snakemake.config['sample'], 2)