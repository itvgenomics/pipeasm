import shutil
import os
import re
import pandas as pd
import zipfile
from io import StringIO


def copy_files(sample):

    files = [
        f"results/Trimming_QC/QC/HiFi_FastQC/{sample}.trimmed_fastqc.html",
        f"results/Trimming_QC/QC/HiFi_FastQC/{sample}.trimmed_fastqc.zip",
        f"results/Trimming_QC/QC/HiFi_NanoPlot/{sample}LengthvsQualityScatterPlot_kde.svg",
        f"results/Trimming_QC/QC/HiFi_NanoPlot/{sample}Non_weightedHistogramReadlength.svg",
        f"results/Trimming_QC/QC/HiFi_NanoPlot/{sample}Non_weightedLogTransformed_HistogramReadlength.svg",
        f"results/Trimming_QC/QC/HiFi_NanoPlot/{sample}WeightedHistogramReadlength.svg",
        f"results/Trimming_QC/QC/HiFi_NanoPlot/{sample}WeightedLogTransformed_HistogramReadlength.svg",
        f"results/Trimming_QC/QC/HiFi_NanoPlot/{sample}Yield_By_Length.svg",
        f"results/Trimming_QC/QC/HiC_FastQC/{sample}_R1.trimmed_paired_fastqc.html",
        f"results/Trimming_QC/QC/HiC_FastQC/{sample}_R1.trimmed_paired_fastqc.zip",
        f"results/Trimming_QC/QC/HiC_FastQC/{sample}_R2.trimmed_paired_fastqc.html",
        f"results/Trimming_QC/QC/HiC_FastQC/{sample}_R2.trimmed_paired_fastqc.zip",
        f"results/Assembly/Genome_Stats/HiFi_GenomeScope2/{sample}_linear_plot.png",
        f"results/Assembly/Genome_Stats/HiFi_GenomeScope2/{sample}_log_plot.png",
        f"results/Assembly/Genome_Stats/HiFi_GenomeScope2/{sample}_transformed_linear_plot.png",
        f"results/Assembly/Genome_Stats/HiFi_GenomeScope2/{sample}_transformed_log_plot.png",
        f"results/Assembly/Genome_Stats/Merqury/Solo_Asm/{sample}.{sample}.p_ctg.spectra-cn.st.png",
        f"results/Assembly/Genome_Stats/Merqury/Solo_Asm/{sample}.{sample}.a_ctg.spectra-cn.st.png",
        f"results/Assembly/Genome_Stats/Merqury/Solo_Asm/{sample}.spectra-asm.st.png",
        f"results/Assembly/Genome_Stats/Merqury/Solo_Asm/{sample}.spectra-cn.st.png",
        f"results/Assembly/Genome_Stats/Merqury/Phased_Asm/{sample}.{sample}.hic.hap1.p_ctg.spectra-cn.st.png",
        f"results/Assembly/Genome_Stats/Merqury/Phased_Asm/{sample}.{sample}.hic.hap2.p_ctg.spectra-cn.st.png",
        f"results/Assembly/Genome_Stats/Merqury/Phased_Asm/{sample}.spectra-asm.st.png",
        f"results/Assembly/Genome_Stats/Merqury/Phased_Asm/{sample}.spectra-cn.st.png",
        f"results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/00-Solo-Hap1/{sample}_Solo_Hap1.snail.png",
        f"results/Assembly/Genome_Stats/SnailPlot/Solo_Asm/01-Solo-Hap2/{sample}_Solo_Hap2.snail.png",
        f"results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap1/{sample}_Phased_Hap1.snail.png",
        f"results/Assembly/Genome_Stats/SnailPlot/Phased_Asm/Phased-Hap2/{sample}_Phased_Hap2.snail.png",
        f"results/Assembly/Mitogenome/Solo_Asm/{sample}.contigs_stats.tsv",
        f"results/Assembly/Mitogenome/Phased_Asm/{sample}.contigs_stats.tsv",
        f"results/Assembly/Genome_Stats/Smudgeplot/{sample}_smudgeplot_log10.pdf",
        f"results/Assembly/Genome_Stats/Smudgeplot/{sample}_smudgeplot.pdf",
        f"results/Assembly/Genome_Stats/KAT/{sample}.mx.png",
        f"results/Scaffolding/Initial_Contacts/Hap1/{sample}_Hap1FullMap.png",
        f"results/Scaffolding/Initial_Contacts/Hap2/{sample}_Hap2FullMap.png",
        f"results/Scaffolding/Final_Contacts/Hap1/{sample}.final_contacts_Hap1FullMap.png",
        f"results/Scaffolding/Final_Contacts/Hap2/{sample}.final_contacts_Hap2FullMap.png",
        f"results/Scaffolding/Scaffolding_stats/SnailPlot/Hap1/{sample}_Scaffolding_Hap1.snail.png",
        f"results/Scaffolding/Scaffolding_stats/SnailPlot/Hap2/{sample}_Scaffolding_Hap2.snail.png",
        f"results/Trimming_QC/QC/ONT_NanoPlot/{sample}LengthvsQualityScatterPlot_dot.svg",
        f"results/Trimming_QC/QC/ONT_NanoPlot/{sample}LengthvsQualityScatterPlot_kde.svg",
        f"results/Trimming_QC/QC/ONT_NanoPlot/{sample}Non_weightedHistogramReadlength.svg",
        f"results/Trimming_QC/QC/ONT_NanoPlot/{sample}Non_weightedLogTransformed_HistogramReadlength.svg",
        f"results/Trimming_QC/QC/ONT_NanoPlot/{sample}WeightedHistogramReadlength.svg",
        f"results/Trimming_QC/QC/ONT_NanoPlot/{sample}WeightedLogTransformed_HistogramReadlength.svg",
        f"results/Trimming_QC/QC/ONT_NanoPlot/{sample}Yield_By_Length.svg",
        f"results/Scaffolding/Scaffolding_stats/Merqury/{sample}.completeness.stats",
		f"results/Scaffolding/Scaffolding_stats/Merqury/{sample}.{sample}.yahs_hap1.spectra-cn.st.png",
		f"results/Scaffolding/Scaffolding_stats/Merqury/{sample}.{sample}.yahs_hap2.spectra-cn.st.png",
		f"results/Scaffolding/Scaffolding_stats/Merqury/{sample}.spectra-asm.st.png",
		f"results/Scaffolding/Scaffolding_stats/Merqury/{sample}.spectra-cn.st.png"
    ]

    directories = [
        "workflow/report/Images",
        "workflow/report/Files",
    ]

    for dir in directories:
        if not os.path.exists(dir):
            print(f"Creating Directory: {dir}")
            os.makedirs(dir)

    for file in files:
        if os.path.exists(file):
            print(f"Copying file: {file}")
            if ".png" in file or ".svg" in file or ".pdf" in file:
                file_name = os.path.basename(file)
                if "HiFi" in file:
                    shutil.copy(file, f"workflow/report/Images/HiFi_{file_name}")
                elif "HiC" in file:
                    shutil.copy(file, f"workflow/report/Images/HiC_{file_name}")
                elif "ONT" in file:
                    shutil.copy(file, f"workflow/report/Images/ONT_{file_name}")
                elif "Solo_Asm" in file:
                    shutil.copy(file, f"workflow/report/Images/Solo_{file_name}")
                elif "Phased_Asm" in file:
                    shutil.copy(file, f"workflow/report/Images/Phased_{file_name}")
                elif "Scaffolding_stats/Merqury" in file:
                    shutil.copy(file, f"workflow/report/Images/Scaffolding_{file_name}")
                else:
                    shutil.copy(file, "workflow/report/Images")
            else:
                shutil.copy(file, "workflow/report/Files")


def parse_fastqc(sample):
    fastqc_files = [
        f"results/Trimming_QC/QC/HiFi_FastQC/{sample}.trimmed_fastqc.zip",
        f"results/Trimming_QC/QC/HiC_FastQC/{sample}_R1.trimmed_paired_fastqc.zip",
        f"results/Trimming_QC/QC/HiC_FastQC/{sample}_R2.trimmed_paired_fastqc.zip",
    ]

    basic_statistics = []
    overrepresented_sequences = []

    for fastqc_file in fastqc_files:
        if os.path.exists(fastqc_file):
            # Open the zip file
            with zipfile.ZipFile(fastqc_file, "r") as zip_ref:
                file_name = fastqc_file.split("/")[-1].split(".zip")[0]
                # Read the contents of a specific file in the zip archive
                with zip_ref.open(f"{file_name}/fastqc_data.txt") as file:
                    # Read the content of the file
                    content = file.read().decode("utf-8")

                    # Define the pattern to match everything between ">>Basic Statistics" and ">>END_MODULE"
                    bs_pattern = re.compile(
                        r">>Basic Statistics(.*?)>>END_MODULE", re.DOTALL
                    )

                    # Search for the pattern in the text
                    bs_match = bs_pattern.search(content)

                    if bs_match:
                        # Extract the matched substring
                        extracted_data = bs_match.group(1).strip()

                    bs_data = {}

                    for line in extracted_data.split("\n"):
                        if (
                            "pass" in line
                            or "#" in line
                            or "fail" in line
                            or "warn" in line
                        ):
                            pass
                        else:
                            key, value = line.strip().split("\t")
                            bs_data[key] = value

                    basic_statistics.append(bs_data)

                    if "paired" in fastqc_file:
                        # Define the pattern to match everything between ">>Overrepresented sequences" and ">>END_MODULE"
                        os_pattern = re.compile(
                            r">>Overrepresented sequences(.*?)>>END_MODULE", re.DOTALL
                        )

                        # Search for the pattern in the text
                        os_match = os_pattern.search(content)

                        if os_match:
                            # Extract the matched substring
                            extracted_data = os_match.group(1).strip()

                        for line in extracted_data.split("\n"):
                            os_data = {}
                            if "pass" in line:
                                break
                            elif "#" in line or "fail" in line or "warn" in line:
                                pass
                            else:
                                Sequence, Count, Percentage, Possible_Source = (
                                    line.strip().split("\t")
                                )
                                os_data["Sequence"] = Sequence
                                os_data["Count"] = Count
                                os_data["Percentage"] = Percentage
                                os_data["Possible Source"] = Possible_Source

                                overrepresented_sequences.append(os_data)

    df_bs = pd.DataFrame(basic_statistics)
    df_bs.to_csv(f"workflow/report/{sample}.fastqc_basic_stats.csv", index=False)

    df_os = pd.DataFrame(overrepresented_sequences)
    df_os.to_csv(f"workflow/report/{sample}.fastqc_overrepresented_sequences.csv", index=False)


def parse_nanoplot(sample):
    hifi_nanoplot_file = f"results/Trimming_QC/QC/HiFi_NanoPlot/{sample}NanoStats.txt"
    if os.path.exists(hifi_nanoplot_file):
        # Read the text file
        with open(hifi_nanoplot_file, "r") as file:
            lines = file.readlines()

        # Initialize empty lists to store data
        columns = []
        data = []

        # Process each line in the file
        for line in lines:
            # Split the line by tabs
            items = line.strip().split("\t")

            # If it's the first line, it contains column names
            if not columns:
                columns = items
            else:
                # Append the data to the list
                data.append(items)

        # Create a DataFrame from the data
        df = pd.DataFrame(data, columns=columns)
        df.to_csv(f"workflow/report/{sample}.hifi_nanoplot.csv", index=False)

    ont_nanoplot_file = f"results/Trimming_QC/QC/ONT_NanoPlot/{sample}NanoStats.txt"
    if os.path.exists(ont_nanoplot_file):
        # Read the text file
        with open(ont_nanoplot_file, "r") as file:
            lines = file.readlines()

        # Initialize empty lists to store data
        columns = []
        data = []

        # Process each line in the file
        for line in lines:
            # Split the line by tabs
            items = line.strip().split("\t")

            # If it's the first line, it contains column names
            if not columns:
                columns = items
            else:
                # Append the data to the list
                data.append(items)

        # Create a DataFrame from the data
        df = pd.DataFrame(data, columns=columns)
        df.to_csv(f"workflow/report/{sample}.ont_nanoplot.csv", index=False)

def parse_genomescope2(sample):
    genomescope2_file = f"results/Assembly/Genome_Stats/HiFi_GenomeScope2/{sample}_model.txt"

    if os.path.exists(genomescope2_file):
        # Read the text from the file
        with open(genomescope2_file, "r") as file:
            text = file.read()

        # Define a regular expression pattern to match parameter lines
        pattern = r"^\s*(\w+)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+(<\d*e[-+]?\d+|\d*\.\d+|\d+)"

        # Initialize empty lists to store data
        parameters = []

        # Find all matches in the text
        matches = re.findall(pattern, text, re.MULTILINE)

        # Process each match and append to parameters list
        for match in matches:
            parameters.append(
                {
                    "Parameter": match[0],
                    "Estimate": float(match[1]),
                    "Std. Error": float(match[2]),
                    "t value": float(match[3]),
                    "Pr(>|t|)": str(match[4]),
                }
            )

        # Create a DataFrame from the parameters
        df = pd.DataFrame(parameters)

        # Display the DataFrame
        df.to_csv(f"workflow/report/{sample}.genomescope2_model.csv", index=False)

        # Read the text from the file
        with open(
            f"results/Assembly/Genome_Stats/HiFi_GenomeScope2/{sample}_summary.txt",
            "r",
        ) as file:
            content = file.read()

        # Replace "(aa)" and "(ab)" with empty strings
        content = re.sub(r"\(\w+\)", "", content)

        # Define the pattern to match the lines below "property, min, and max"
        pattern = (
            r"(\w+(?:\s\w+)*)\s+(\d[\d,.]*%?|\d[\d,.]*\s\w+)\s+(\d[\d,.]*%?|\d[\d,.]*\s\w+)"
        )

        # Find all matches using the pattern
        matches = re.findall(pattern, content)

        df = pd.DataFrame(matches, columns=["property", "min", "max"])

        df.to_csv(f"workflow/report/{sample}.genomescope2_summary.csv", index=False)


def parse_merqury(sample):
    # Initialize empty DataFrames for completeness and QV data
    df_solo_completeness = pd.DataFrame()
    df_phased_completeness = pd.DataFrame()
    df_solo_qv = pd.DataFrame()
    df_phased_qv = pd.DataFrame()

    # Define common column names for completeness and QV files
    completeness_columns = ["asm", "kmer_set", "kmer_in_asm", "kmer_total", "completeness(%)"]
    qv_columns = ["asm", "kmer_asm_uniq", "kmer_both", "qv", "error_rate"]

    # File paths
    solo_completeness_path = f"results/Assembly/Genome_Stats/Merqury/Solo_Asm/{sample}.completeness.stats"
    phased_completeness_path = f"results/Assembly/Genome_Stats/Merqury/Phased_Asm/{sample}.completeness.stats"
    solo_qv_path = f"results/Assembly/Genome_Stats/Merqury/Solo_Asm/{sample}.qv"
    phased_qv_path = f"results/Assembly/Genome_Stats/Merqury/Phased_Asm/{sample}.qv"

    # Read HiFi completeness data if the file exists
    if os.path.exists(solo_completeness_path):
        df_solo_completeness = pd.read_csv(
            solo_completeness_path, sep="\t", header=None, names=completeness_columns
        )
        df_solo_completeness["asm"] = df_solo_completeness["asm"].replace("both", "both_solo")

    # Read HiC completeness data if the file exists
    if os.path.exists(phased_completeness_path):
        df_phased_completeness = pd.read_csv(
            phased_completeness_path, sep="\t", header=None, names=completeness_columns
        )
        df_phased_completeness["asm"] = df_phased_completeness["asm"].replace("both", "both_phased")

    # Concatenate completeness DataFrames
    df_completeness = pd.concat([df_solo_completeness, df_phased_completeness], ignore_index=True)

    # Read HiFi QV data if the file exists
    if os.path.exists(solo_qv_path):
        df_solo_qv = pd.read_csv(
            solo_qv_path, sep="\t", header=None, names=qv_columns
        )
        df_solo_qv["asm"] = df_solo_qv["asm"].replace("Both", "both_solo")

    # Read HiC QV data if the file exists
    if os.path.exists(phased_qv_path):
        df_phased_qv = pd.read_csv(
            phased_qv_path, sep="\t", header=None, names=qv_columns
        )
        df_phased_qv["asm"] = df_phased_qv["asm"].replace("Both", "both_phased")

    # Concatenate QV DataFrames
    df_qv = pd.concat([df_solo_qv, df_phased_qv], ignore_index=True)

    # Merge completeness and QV DataFrames
    if not df_completeness.empty and not df_qv.empty:
        df_completo = pd.merge(df_completeness, df_qv, on="asm", how="left")
    else:
        df_completo = df_completeness if not df_completeness.empty else df_qv

    # Save the final merged DataFrame to CSV
    df_completo.to_csv(f"workflow/report/{sample}.merqury.csv", index=False)

def parse_compleasm(sample):
    all_data = []

    if os.path.exists(
        f"results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Primary/{sample}.summary.txt"
    ):
        # Open the text file and read line by line
        with open(
            f"results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Primary/{sample}.summary.txt",
            "r",
        ) as file:

            # Initialize empty lists to store data
            s_values = []
            d_values = []
            f_values = []
            m_values = []
            n_values = []

            for line in file:
                if line.startswith("## lineage:"):
                    # Extract file name
                    asm = "Solo_hap1"
                elif line.startswith(("S:", "D:", "F:", "M:", "N:")):
                    # Extract values
                    # parts = line.strip().split(",")
                    value = line.split(":")[1].replace("\n", "")
                    if line.startswith("S:"):
                        s_values.append(value)
                    elif line.startswith("D:"):
                        d_values.append(value)
                    elif line.startswith("F:"):
                        f_values.append(value)
                    elif line.startswith("M:"):
                        m_values.append(value)
                    elif line.startswith("N:"):
                        n_values.append(value)

        # Create pandas DataFrame
        data = {
            "asm": asm,
            "Complete and Single Copy": s_values,
            "Complete and Duplicated:": d_values,
            "Fragmented:": f_values,
            "Missing:": m_values,
            "Total Genes:": n_values,
        }

        all_data.append(data)

    if os.path.exists(
        f"results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Alt/{sample}.summary.txt"
    ):
        # Open the text file and read line by line
        with open(
            f"results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Alt/{sample}.summary.txt",
            "r",
        ) as file:

            # Initialize empty lists to store data
            s_values = []
            d_values = []
            f_values = []
            m_values = []
            n_values = []

            for line in file:
                if line.startswith("## lineage:"):
                    # Extract file name
                    asm = "Solo_hap2"
                elif line.startswith(("S:", "D:", "F:", "M:", "N:")):
                    # Extract values
                    # parts = line.strip().split(",")
                    value = line.split(":")[1].replace("\n", "")
                    if line.startswith("S:"):
                        s_values.append(value)
                    elif line.startswith("D:"):
                        d_values.append(value)
                    elif line.startswith("F:"):
                        f_values.append(value)
                    elif line.startswith("M:"):
                        m_values.append(value)
                    elif line.startswith("N:"):
                        n_values.append(value)

        # Create pandas DataFrame
        data = {
            "asm": asm,
            "Complete and Single Copy": s_values,
            "Complete and Duplicated:": d_values,
            "Fragmented:": f_values,
            "Missing:": m_values,
            "Total Genes:": n_values,
        }

        all_data.append(data)

    if os.path.exists(
        f"results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap1/{sample}.summary.txt"
    ):
        # Open the text file and read line by line
        with open(
            f"results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap1/{sample}.summary.txt",
            "r",
        ) as file:

            # Initialize empty lists to store data
            s_values = []
            d_values = []
            f_values = []
            m_values = []
            n_values = []

            for line in file:
                if line.startswith("## lineage:"):
                    # Extract file name
                    asm = "Phased_hap1"
                elif line.startswith(("S:", "D:", "F:", "M:", "N:")):
                    # Extract values
                    # parts = line.strip().split(",")
                    value = line.split(":")[1].replace("\n", "")
                    if line.startswith("S:"):
                        s_values.append(value)
                    elif line.startswith("D:"):
                        d_values.append(value)
                    elif line.startswith("F:"):
                        f_values.append(value)
                    elif line.startswith("M:"):
                        m_values.append(value)
                    elif line.startswith("N:"):
                        n_values.append(value)

        # Create pandas DataFrame
        data = {
            "asm": asm,
            "Complete and Single Copy": s_values,
            "Complete and Duplicated:": d_values,
            "Fragmented:": f_values,
            "Missing:": m_values,
            "Total Genes:": n_values,
        }

        all_data.append(data)

    if os.path.exists(
        f"results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap2/{sample}.summary.txt"
    ):
        # Open the text file and read line by line
        with open(
            f"results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap2/{sample}.summary.txt",
            "r",
        ) as file:

            # Initialize empty lists to store data
            s_values = []
            d_values = []
            f_values = []
            m_values = []
            n_values = []

            for line in file:
                if line.startswith("## lineage:"):
                    # Extract file name
                    asm = "Phased_hap2"
                elif line.startswith(("S:", "D:", "F:", "M:", "N:")):
                    # Extract values
                    # parts = line.strip().split(",")
                    value = line.split(":")[1].replace("\n", "")
                    if line.startswith("S:"):
                        s_values.append(value)
                    elif line.startswith("D:"):
                        d_values.append(value)
                    elif line.startswith("F:"):
                        f_values.append(value)
                    elif line.startswith("M:"):
                        m_values.append(value)
                    elif line.startswith("N:"):
                        n_values.append(value)

        # Create pandas DataFrame
        data = {
            "asm": asm,
            "Complete and Single Copy": s_values,
            "Complete and Duplicated:": d_values,
            "Fragmented:": f_values,
            "Missing:": m_values,
            "Total Genes:": n_values,
        }

        all_data.append(data)

    if os.path.exists(
        f"results/Scaffolding/Scaffolding_stats/Compleasm/Hap1/{sample}.summary.txt"
    ):
        # Open the text file and read line by line
        with open(
            f"results/Scaffolding/Scaffolding_stats/Compleasm/Hap1/{sample}.summary.txt",
            "r",
        ) as file:

            # Initialize empty lists to store data
            s_values = []
            d_values = []
            f_values = []
            m_values = []
            n_values = []

            for line in file:
                if line.startswith("## lineage:"):
                    # Extract file name
                    asm = "Scaffolding_hap1"
                elif line.startswith(("S:", "D:", "F:", "M:", "N:")):
                    # Extract values
                    # parts = line.strip().split(",")
                    value = line.split(":")[1].replace("\n", "")
                    if line.startswith("S:"):
                        s_values.append(value)
                    elif line.startswith("D:"):
                        d_values.append(value)
                    elif line.startswith("F:"):
                        f_values.append(value)
                    elif line.startswith("M:"):
                        m_values.append(value)
                    elif line.startswith("N:"):
                        n_values.append(value)

        # Create pandas DataFrame
        data = {
            "asm": asm,
            "Complete and Single Copy": s_values,
            "Complete and Duplicated:": d_values,
            "Fragmented:": f_values,
            "Missing:": m_values,
            "Total Genes:": n_values,
        }

        all_data.append(data)

    if os.path.exists(
        f"results/Scaffolding/Scaffolding_stats/Compleasm/Hap2/{sample}.summary.txt"
    ):
        # Open the text file and read line by line
        with open(
            f"results/Scaffolding/Scaffolding_stats/Compleasm/Hap2/{sample}.summary.txt",
            "r",
        ) as file:

            # Initialize empty lists to store data
            s_values = []
            d_values = []
            f_values = []
            m_values = []
            n_values = []

            for line in file:
                if line.startswith("## lineage:"):
                    # Extract file name
                    asm = "Scaffolding_hap2"
                elif line.startswith(("S:", "D:", "F:", "M:", "N:")):
                    # Extract values
                    # parts = line.strip().split(",")
                    value = line.split(":")[1].replace("\n", "")
                    if line.startswith("S:"):
                        s_values.append(value)
                    elif line.startswith("D:"):
                        d_values.append(value)
                    elif line.startswith("F:"):
                        f_values.append(value)
                    elif line.startswith("M:"):
                        m_values.append(value)
                    elif line.startswith("N:"):
                        n_values.append(value)

        # Create pandas DataFrame
        data = {
            "asm": asm,
            "Complete and Single Copy": s_values,
            "Complete and Duplicated:": d_values,
            "Fragmented:": f_values,
            "Missing:": m_values,
            "Total Genes:": n_values,
        }

        all_data.append(data)

    df = pd.DataFrame(all_data)

    df.to_csv(f"workflow/report/{sample}.compleasm.csv", index=False)


def abstract(sample):
    abstract_text = ""

    if os.path.exists(f"results/Assembly/Genome_Stats/GFAstats/{sample}.p_ctg.fa.stats") and os.path.exists(f"results/Assembly/Genome_Stats/GFAstats/{sample}.a_ctg.fa.stats"):
        data_hap1 = []
        with open(
            f"results/Assembly/Genome_Stats/GFAstats/{sample}.p_ctg.fa.stats",
            "r",
        ) as file:
            # Iterate through each line
            for line in file:
                # Use regular expression to search for lines containing patterns
                match = re.search(r"# scaffolds: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap1.append(number)
                match = re.search(r"Scaffold N50: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap1.append(number)
                match = re.search(r"Total scaffold length: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap1.append(number)
                match = re.search(r"# contigs: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap1.append(number)
                match = re.search(r"Largest scaffold: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap1.append(number)

        data_hap2 = []
        with open(
            f"results/Assembly/Genome_Stats/GFAstats/{sample}.a_ctg.fa.stats",
            "r",
        ) as file:
            # Iterate through each line
            for line in file:
                # Use regular expression to search for lines containing patterns
                match = re.search(r"# scaffolds: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap2.append(number)
                match = re.search(r"Scaffold N50: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap2.append(number)
                match = re.search(r"Total scaffold length: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap2.append(number)
                match = re.search(r"# contigs: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap2.append(number)
                match = re.search(r"Largest scaffold: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap2.append(number)

        data_busco = []
        with open(
            f"results/Assembly/Genome_Stats/Compleasm/Solo_Asm_Primary/{sample}.summary.txt",
            "r",
        ) as file:

            # Iterate through each line
            for line in file:
                # Use regular expression to search for lines containing "S:", "D:", "F:", "M:", "N:"
                match = re.search(r"(S|D|F|M|N):(\d+\.\d+)%", line)
                if match:
                    # Extract the percentage and the letter (S, D, F, M, or N)
                    percentage = float(match.group(2))
                    # Append the percentage to the percentages list
                    data_busco.append(percentage)

                if "N:" in line:
                    n_number = line.split("N:")[1].replace("\n", "")
                    data_busco.append(n_number)

        data_qv = []
        with open(
            f"results/Assembly/Genome_Stats/Merqury/Solo_Asm/{sample}.qv", "r"
        ) as file:
            for line in file:
                data_qv = line.split("\t")
                break

        abstract_text = f"""--- Solo Assembly ---
N Scaffolds -> Primary: {data_hap1[0]} Alternate: {data_hap2[0]}
Assembled Bases (bp) -> Primary: {data_hap1[1]} Alternate: {data_hap2[1]}
N50 -> Primary: {data_hap1[2]} Alternate: {data_hap2[2]}
N Contigs -> Primary: {data_hap1[3]} Alternate: {data_hap2[3]}
Largest scaffold (bp) -> Primary: {data_hap1[4]} Alternate: {data_hap2[4]}
Busco (Prim. only) -> C(%): {data_busco[0]} D(%): {data_busco[1]} F(%): {data_busco[2]} M(%): {data_busco[3]} Genes:  {data_busco[4]}
QV (Prim. only) -> {data_qv[3]}
"""
    with open(f"workflow/report/{sample}.abstract.txt", "w") as output:
        output.write(abstract_text)

    if os.path.exists(
        f"results/Assembly/Genome_Stats/GFAstats/{sample}.hic.hap1.p_ctg.fa.stats"
    ):
        data_hap1 = []
        with open(
            f"results/Assembly/Genome_Stats/GFAstats/{sample}.hic.hap1.p_ctg.fa.stats",
            "r",
        ) as file:
            # Iterate through each line
            for line in file:
                # Use regular expression to search for lines containing patterns
                match = re.search(r"# scaffolds: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap1.append(number)
                match = re.search(r"Scaffold N50: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap1.append(number)
                match = re.search(r"Total scaffold length: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap1.append(number)
                match = re.search(r"# contigs: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap1.append(number)
                match = re.search(r"Largest scaffold: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap1.append(number)

        data_hap2 = []
        with open(
            f"results/Assembly/Genome_Stats/GFAstats/{sample}.hic.hap2.p_ctg.fa.stats",
            "r",
        ) as file:
            # Iterate through each line
            for line in file:
                # Use regular expression to search for lines containing patterns
                match = re.search(r"# scaffolds: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap2.append(number)
                match = re.search(r"Scaffold N50: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap2.append(number)
                match = re.search(r"Total scaffold length: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap2.append(number)
                match = re.search(r"# contigs: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap2.append(number)
                match = re.search(r"Largest scaffold: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap2.append(number)

        data_busco = []
        with open(
            f"results/Assembly/Genome_Stats/Compleasm/Phased_Asm_Hap1/{sample}.summary.txt",
            "r",
        ) as file:

            # Iterate through each line
            for line in file:
                # Use regular expression to search for lines containing "S:", "D:", "F:", "M:", "N:"
                match = re.search(r"(S|D|F|M|N):(\d+\.\d+)%", line)
                if match:
                    # Extract the percentage and the letter (S, D, F, M, or N)
                    percentage = float(match.group(2))
                    # Append the percentage to the percentages list
                    data_busco.append(percentage)

                if "N:" in line:
                    n_number = line.split("N:")[1].replace("\n", "")
                    data_busco.append(n_number)

        data_qv = []
        with open(
            f"results/Assembly/Genome_Stats/Merqury/Phased_Asm/{sample}.qv", "r"
        ) as file:
            for line in file:
                data_qv = line.split("\t")
                break

        abstract_text = f"""--- Phased Assembly ---
N Scaffolds -> Hap1: {data_hap1[0]} Hap2: {data_hap2[0]}
Assembled Bases (bp) -> Hap1: {data_hap1[1]} Hap2: {data_hap2[1]}
N50 -> Hap1: {data_hap1[2]} Hap2: {data_hap2[2]}
N Contigs -> Hap1: {data_hap1[3]} Hap2: {data_hap2[3]}
Largest scaffold (bp) -> Hap1: {data_hap1[4]} Hap2: {data_hap2[4]}
Busco (Hap1. only) -> C(%): {data_busco[0]} D(%): {data_busco[1]} F(%): {data_busco[2]} M(%): {data_busco[3]} Genes:  {data_busco[4]}
QV (Hap1. only) -> {data_qv[3]}
"""

        with open(f"workflow/report/{sample}.abstract.txt", "a") as output:
            output.write(abstract_text)

    if os.path.exists(
        f"results/Scaffolding/Scaffolding_stats/GFAstats/{sample}.yahs_scaffolds_hap1.fa.stats"
    ):
        data_hap1 = []
        with open(
            f"results/Scaffolding/Scaffolding_stats/GFAstats/{sample}.yahs_scaffolds_hap1.fa.stats",
            "r",
        ) as file:
            # Iterate through each line
            for line in file:
                # Use regular expression to search for lines containing patterns
                match = re.search(r"# scaffolds: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap1.append(number)
                match = re.search(r"Scaffold N50: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap1.append(number)
                match = re.search(r"Total scaffold length: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap1.append(number)
                match = re.search(r"# contigs: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap1.append(number)
                match = re.search(r"Largest scaffold: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap1.append(number)

        data_hap2 = []
        with open(
            f"results/Scaffolding/Scaffolding_stats/GFAstats/{sample}.yahs_scaffolds_hap2.fa.stats",
            "r",
        ) as file:
            # Iterate through each line
            for line in file:
                # Use regular expression to search for lines containing patterns
                match = re.search(r"# scaffolds: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap2.append(number)
                match = re.search(r"Scaffold N50: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap2.append(number)
                match = re.search(r"Total scaffold length: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap2.append(number)
                match = re.search(r"# contigs: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap2.append(number)
                match = re.search(r"Largest scaffold: (\d+)", line)
                if match:
                    number = match.group(1)
                    data_hap2.append(number)

        data_busco = []
        with open(
            f"results/Scaffolding/Scaffolding_stats/Compleasm/Hap1/{sample}.summary.txt",
            "r",
        ) as file:

            # Iterate through each line
            for line in file:
                # Use regular expression to search for lines containing "S:", "D:", "F:", "M:", "N:"
                match = re.search(r"(S|D|F|M|N):(\d+\.\d+)%", line)
                if match:
                    # Extract the percentage and the letter (S, D, F, M, or N)
                    percentage = float(match.group(2))
                    # Append the percentage to the percentages list
                    data_busco.append(percentage)

                if "N:" in line:
                    n_number = line.split("N:")[1].replace("\n", "")
                    data_busco.append(n_number)

        data_qv = []
        with open(
            f"results/Scaffolding/Scaffolding_stats/Merqury/{sample}.qv", "r"
        ) as file:
            for line in file:
                data_qv = line.split("\t")
                break

        abstract_text = f"""--- After YASH Scaffolding ---
N Scaffolds -> Hap1: {data_hap1[0]} Hap2: {data_hap2[0]}
Assembled Bases (bp) -> Hap1: {data_hap1[1]} Hap2: {data_hap2[1]}
N50 -> Hap1: {data_hap1[2]} Hap2: {data_hap2[2]}
N Contigs -> Hap1: {data_hap1[3]} Hap2: {data_hap2[3]}
Largest scaffold (bp) -> Hap1: {data_hap1[4]} Hap2: {data_hap2[4]}
Busco (Hap1. only) -> C(%): {data_busco[0]} D(%): {data_busco[1]} F(%): {data_busco[2]} M(%): {data_busco[3]} Genes:  {data_busco[4]}
QV (Hap1. only) -> {data_qv[3]}
"""

        with open(f"workflow/report/{sample}.abstract.txt", "a") as output:
            output.write(abstract_text)

def parse_gfastats(sample):
    df_completo = pd.DataFrame()
    # Define a regular expression pattern to extract key-value pairs
    pattern = re.compile(r"([^:]+):\s*([^:]+)")

    columns_to_keep = [
        "asm",
        "# scaffolds",
        "Total scaffold length",
        "Scaffold N50",
        "Scaffold L50",
        "Largest scaffold",
        "Smallest scaffold",
        "# contigs",
        "Total contig length",
        "Contig N50",
        "Contig L50",
        "Largest contig",
        "Smallest contig",
        "# gaps in scaffolds",
        "Total gap length in scaffolds",
        "Gap N50 in scaffolds",
        "Gap L50 in scaffolds",
        "Largest gap in scaffolds",
        "Smallest gap in scaffolds",
        "# segments",
        "Total segment length",
    ]

    if os.path.exists(f"results/Assembly/Genome_Stats/GFAstats/{sample}.bp.hap1.p_ctg.fa.stats") and os.path.exists(f"results/Assembly/Genome_Stats/GFAstats/{sample}.a_ctg.fa.stats"):

        # Read the text file
        with open(
            f"results/Assembly/Genome_Stats/GFAstats/{sample}.bp.hap1.p_ctg.fa.stats",
            "r",
        ) as file:
            lines = file.readlines()

        # Define a regular expression pattern to extract key-value pairs
        pattern = re.compile(r"([^:]+):\s*([^:]+)")

        # Initialize a dictionary to store extracted data
        data = {}

        # Iterate through each line and extract key-value pairs
        data["asm"] = "Solo_Hap1"
        for line in lines:
            match = pattern.match(line)
            if match:
                key = match.group(1).strip()
                value = match.group(2).strip()
                data[key] = value

        # Create a DataFrame from the extracted data
        df_hifi_hap1 = pd.DataFrame([data])

        columns_to_keep = [
            "asm",
            "# scaffolds",
            "Total scaffold length",
            "Scaffold N50",
            "Scaffold L50",
            "Largest scaffold",
            "Smallest scaffold",
            "# contigs",
            "Total contig length",
            "Contig N50",
            "Contig L50",
            "Largest contig",
            "Smallest contig",
            "# gaps in scaffolds",
            "Total gap length in scaffolds",
            "Gap N50 in scaffolds",
            "Gap L50 in scaffolds",
            "Largest gap in scaffolds",
            "Smallest gap in scaffolds",
            "# segments",
            "Total segment length",
        ]

        df_hifi_hap1 = df_hifi_hap1[columns_to_keep]

        # Read the text file
        with open(
            f"results/Assembly/Genome_Stats/GFAstats/{sample}.a_ctg.fa.stats",
            "r",
        ) as file:
            lines = file.readlines()

        # Initialize a dictionary to store extracted data
        data = {}

        # Iterate through each line and extract key-value pairs
        data["asm"] = "Solo_Hap1"
        for line in lines:
            match = pattern.match(line)
            if match:
                key = match.group(1).strip()
                value = match.group(2).strip()
                data[key] = value

        # Create a DataFrame from the extracted data
        df_hifi_hap2 = pd.DataFrame([data])

        df_hifi_hap2 = df_hifi_hap2[columns_to_keep]

        # Concat both dataframes
        df_completo = pd.concat([df_hifi_hap1, df_hifi_hap2])

    if os.path.exists(
        f"results/Assembly/Genome_Stats/GFAstats/{sample}.hic.hap1.p_ctg.fa.stats"
    ):
        # Read the text file
        with open(
            f"results/Assembly/Genome_Stats/GFAstats/{sample}.hic.hap1.p_ctg.fa.stats",
            "r",
        ) as file:
            lines = file.readlines()

        # Initialize a dictionary to store extracted data
        data = {}

        # Iterate through each line and extract key-value pairs
        data["asm"] = "Phased_Hap1"
        for line in lines:
            match = pattern.match(line)
            if match:
                key = match.group(1).strip()
                value = match.group(2).strip()
                data[key] = value

        # Create a DataFrame from the extracted data
        df_hic_hap1 = pd.DataFrame([data])

        df_hic_hap1 = df_hic_hap1[columns_to_keep]

        df_completo = pd.concat([df_completo, df_hic_hap1])

        with open(
            f"results/Assembly/Genome_Stats/GFAstats/{sample}.hic.hap2.p_ctg.fa.stats",
            "r",
        ) as file:
            lines = file.readlines()

        # Initialize a dictionary to store extracted data
        data = {}

        # Iterate through each line and extract key-value pairs
        data["asm"] = "Phased_Hap2"
        for line in lines:
            match = pattern.match(line)
            if match:
                key = match.group(1).strip()
                value = match.group(2).strip()
                data[key] = value

        # Create a DataFrame from the extracted data
        df_hic_hap2 = pd.DataFrame([data])

        df_hic_hap2 = df_hic_hap2[columns_to_keep]

        df_completo = pd.concat([df_completo, df_hic_hap2])

    if os.path.exists(
        f"results/Scaffolding/Scaffolding_stats/GFAstats/{sample}.yahs_scaffolds_hap1.fa.stats"
    ):
        # Read the text file
        with open(
            f"results/Scaffolding/Scaffolding_stats/GFAstats/{sample}.yahs_scaffolds_hap1.fa.stats",
            "r",
        ) as file:
            lines = file.readlines()

        # Initialize a dictionary to store extracted data
        data = {}

        # Iterate through each line and extract key-value pairs
        data["asm"] = "Scaffolding_Hap1"
        for line in lines:
            match = pattern.match(line)
            if match:
                key = match.group(1).strip()
                value = match.group(2).strip()
                data[key] = value

        # Create a DataFrame from the extracted data
        df_scaf_hap1 = pd.DataFrame([data])

        df_scaf_hap1 = df_scaf_hap1[columns_to_keep]

        df_completo = pd.concat([df_completo, df_scaf_hap1])

        with open(
            f"results/Scaffolding/Scaffolding_stats/GFAstats/{sample}.yahs_scaffolds_hap2.fa.stats",
            "r",
        ) as file:
            lines = file.readlines()

        # Initialize a dictionary to store extracted data
        data = {}

        # Iterate through each line and extract key-value pairs
        data["asm"] = "Scaffolding_Hap2"
        for line in lines:
            match = pattern.match(line)
            if match:
                key = match.group(1).strip()
                value = match.group(2).strip()
                data[key] = value

        # Create a DataFrame from the extracted data
        df_scaf_hap2 = pd.DataFrame([data])

        df_scaf_hap2 = df_scaf_hap2[columns_to_keep]

        df_completo = pd.concat([df_completo, df_scaf_hap2])

    df_completo.to_csv(f"workflow/report/{sample}.gfastats.csv", index=False)

def parse_fcsadaptors(sample):
    files = [
        "results/Decontamination/FCS-Adaptor/Solo_Asm_Primary/fcs_adaptor_report.txt",
        "results/Decontamination/FCS-Adaptor/Solo_Asm_Alt/fcs_adaptor_report.txt",
        "results/Decontamination/FCS-Adaptor/Phased_Asm_Hap1/fcs_adaptor_report.txt",
        "results/Decontamination/FCS-Adaptor/Phased_Asm_Hap2/fcs_adaptor_report.txt",
    ]

    df_completo = pd.DataFrame(
        columns=["#accession", "length", "action", "range", "name", "asm"]
    )

    for file in files:
        if os.path.exists(file):
            df = pd.read_csv(file, sep="\t")
            df.fillna("None", inplace=True)
            df["asm"] = str(file).split("/")[2].split("-")[1]
            df_completo = pd.concat([df_completo, df])

    # Fill remaining NaN values with "None"
    df_completo.fillna("None", inplace=True)

    # Reorder columns with 'asm' as the first column
    cols = list(df_completo.columns)
    cols = ["asm"] + [col for col in cols if col != "asm"]
    df_completo = df_completo[cols]

    df_completo.to_csv(f"workflow/report/{sample}.fcsadaptor.csv", index=False)


def parse_fcsgx(sample, taxid):
    files = [
        f"results/Decontamination/Contaminants/Solo_Asm_Primary/{sample}.p_ctg.{taxid}.fcs_gx_report.txt",
        f"results/Decontamination/Contaminants/Solo_Asm_Alt/{sample}.a_ctg.{taxid}.fcs_gx_report.txt",
        f"results/Decontamination/Contaminants/Phased_Asm_Hap1/{sample}.hic.hap1.p_ctg.{taxid}.fcs_gx_report.txt",
        f"results/Decontamination/Contaminants/Phased_Asm_Hap2/{sample}.hic.hap2.p_ctg.{taxid}.fcs_gx_report.txt",
    ]

    df_completo = pd.DataFrame(
        columns=[
            "asm",
            "#seq_id",
            "start_pos",
            "end_pos",
            "seq_len",
            "action",
            "div",
            "agg_cont_cov",
            "top_tax_name",
        ]
    )

    for file in files:
        if os.path.exists(file):
            df = pd.read_csv(file, sep="\t", skiprows=1)
            assembly = str(file).split("/")[3]
            df.loc[0, 'asm'] = assembly
            df.fillna("None", inplace=True)
            df_completo = pd.concat([df_completo, df])

    # Fill remaining NaN values with "None"
    df_completo.fillna("None", inplace=True)

    # Reorder columns with 'asm' as the first column
    cols = list(df_completo.columns)
    cols = ["asm"] + [col for col in cols if col != "asm"]
    df_completo = df_completo[cols]

    df_completo.to_csv(f"workflow/report/{sample}.fcsgx.csv", index=False)

def parse_flagstats(sample):
    files = [f"results/Scaffolding/Initial_Contacts/Hap1/{sample}.R1.bam.flagstats",
             f"results/Scaffolding/Initial_Contacts/Hap1/{sample}.R2.bam.flagstats",
             f"results/Scaffolding/Initial_Contacts/Hap1/{sample}.merged.bam.flagstats",
             f"results/Scaffolding/Initial_Contacts/Hap2/{sample}.R1.bam.flagstats",
             f"results/Scaffolding/Initial_Contacts/Hap2/{sample}.R2.bam.flagstats",
             f"results/Scaffolding/Initial_Contacts/Hap2/{sample}.merged.bam.flagstats",
             f"results/Scaffolding/Final_Contacts/Hap1/{sample}.R1.bam.flagstats",
             f"results/Scaffolding/Final_Contacts/Hap1/{sample}.R2.bam.flagstats",
             f"results/Scaffolding/Final_Contacts/Hap1/{sample}.merged.bam.flagstats",
             f"results/Scaffolding/Final_Contacts/Hap2/{sample}.R1.bam.flagstats",
             f"results/Scaffolding/Final_Contacts/Hap2/{sample}.R2.bam.flagstats",
             f"results/Scaffolding/Final_Contacts/Hap2/{sample}.merged.bam.flagstats"]

    names = ["Initial_Contacts_Hap1_R1",
             "Initial_Contacts_Hap1_R2",
             "Initial_Contacts_Hap1_Merged",
             "Initial_Contacts_Hap2_R1",
             "Initial_Contacts_Hap2_R2",
             "Initial_Contacts_Hap2_Merged",
             "Final_Contacts_Hap1_R1",
             "Final_Contacts_Hap1_R2",
             "Final_Contacts_Hap1_Merged",
             "Final_Contacts_Hap2_R1",
             "Final_Contacts_Hap2_R2",
             "Final_Contacts_Hap2_Merged"]

    df_completo = pd.DataFrame()

    i = 0

    for file in files:
        if os.path.exists(file):
            # Read the tab-separated file into a pandas DataFrame
            df = pd.read_csv(file, sep="\t", header=None)
            df.columns = ["Value1", "Value2", "Description"]

            # Concatenate 'Value1' and 'Value2' into a new column
            df["Reads"] = df["Value1"].astype(str) + " | " + df["Value2"].astype(str)
            df.drop(["Value1", "Value2"], axis=1, inplace=True)

            # Keep the original order of 'Description' column
            df["Description"] = pd.Categorical(
                df["Description"], categories=df["Description"], ordered=True
            )

            df = df.pivot_table(
                index=None, values=["Reads"], columns="Description", aggfunc="first"
            )

            df.insert(0, "file", names[i])

            i = i + 1

            df_completo = pd.concat([df_completo, df])

    df_completo.to_csv(f"workflow/report/{sample}.flagstats.csv", index=False)

copy_files(snakemake.config["sample"])
parse_fastqc(snakemake.config["sample"])
parse_nanoplot(snakemake.config["sample"])
parse_genomescope2(snakemake.config["sample"])
parse_merqury(snakemake.config["sample"])
parse_compleasm(snakemake.config["sample"])
abstract(snakemake.config["sample"])
parse_gfastats(snakemake.config["sample"])
parse_fcsadaptors(snakemake.config["sample"])
parse_fcsgx(snakemake.config["sample"], snakemake.config["taxid"])
parse_flagstats(snakemake.config["sample"])
