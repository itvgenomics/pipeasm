#!/bin/bash
set -e

cat << 'EOF'

### Pipeasm: a tool for automated large chromosome-scale genome assembly and evaluation
#Version: 1.1.0
#Authors: Bruno Marques Silva, Fernanda de Jesus Trindade, Lucas Eduardo Costa Canesin, Giordano Souza, Alexandre Aleixo, Gisele Nunes, Renato Renison Moreira-Oliveira
#Bioinformatics Advances, Volume 6, Issue 1, 2026, vbaf326, https://doi.org/10.1093/bioadv/vbaf326

###    Copyright (C) 2024  Renato Oliveira
###
###    This program is free software: you can redistribute it and/or modify
###    it under the terms of the GNU General Public License as published by
###    the Free Software Foundation, either version 3 of the License, or
###    any later version.
###
###    This program is distributed in the hope that it will be useful,
###    but WITHOUT ANY WARRANTY; without even the implied warranty of
###    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###    GNU General Public License for more details.
###
###    You should have received a copy of the GNU General Public License
###    along with this program.  If not, see <http://www.gnu.org/licenses/>.

EOF

FULL=true
TRIMMING_QC=false
KMER_EVAL=false
ASSEMBLY=false
SCAFFOLDING=false

SETUNLOCK=""
SETUNP=""
SETSLURM="false"
PARTITION=""

USED_T=false
USED_J=false

while [ "$1" != "" ]; do
    case $1 in
    -c)
        shift
        CONFIGFILE=$1
        ;;
    -d)
        shift
        WORKDIR=$1
        ;;
    -t)
        shift
        THREADS=$1
        USED_T=true
        ;;
    -j)
        shift
        THREADS=$1
        USED_J=true
        ;;
    -s)
        shift
        SNAKEFILE=$1
        ;;
    -trimming_qc)
        TRIMMING_QC=true
        FULL=false
        ;;
    -kmer_eval)
        KMER_EVAL=true
        FULL=false
        ;;
    -assembly)
        ASSEMBLY=true
        FULL=false
        ;;
    -scaffolding)
        SCAFFOLDING=true
        FULL=false
        ;;
    -slurm)
        SETSLURM=true
        ;;
    -partition)
        shift
        PARTITION=$1
        ;;
    -np)
        SETNP="-np"
        ;;
    -unlock)
        SETUNLOCK="--unlock"
        ;;
    *)
        echo "ERROR: Invalid argument definition. Ensure there are no spaces within your strings and that all arguments are correctly specified."
        exit 1
        ;;
    esac
    shift
done

if [ "$SETSLURM" = true ]; then
    if [ "$USED_J" = false ]; then
        echo "ERROR: When using -slurm, you must use -j (not -t) to set the number of jobs submitted to the queue."
        exit 1
    fi
    if [ "$USED_T" = true ]; then
        echo "ERROR: -t is not allowed with -slurm. Use -j instead."
        exit 1
    fi
fi

# Ensure THREADS and WORKDIR have values
if [ -z "$WORKDIR" ] || [ -z "$THREADS" ]; then
    cat << 'EOF'
ERROR: Missing one of the required arguments: -d (Work Directory), -t (# Threads)

## Usage: ./Pipeasm.sh -c <config.yaml> -d </path/to/work/dir> -s </path/to/snakefile> -t <# threads>
#-d </path/to/work/dir> = Path to your working directory where all the workflow file are
#-c </path/to/config.yaml> = Overwrite the default configuration file with all nedded parameters (config/config.yaml).
#-s </path/to/snakefile> = Overwrite the default snakefile path (workflow/Snakefile)
#-t <int> = Number of threads to use
#-j <int> = Number of jobs submitted to the queue (only for -slurm)

# You can choose a Pipeasm step with:
    -trimming_qc (for Trimming and Quality Control);
    -kmer_eval (for k-mer Evaluation stats/plots);
    -assembly (for all Assembly and Decontamination steps);
    -scaffolding (for all Scaffolding and Hi-C Map steps);

# Only the -d and -t flags are required if you want to use the Snakemake default parameters
# You can perform a dry-run (build only the DAG, no rule will be run) with -np and unlock the directory, if needed, with --unlock

EOF

    exit 1
fi

if [ "$SETSLURM" = true ] && [ -z "$PARTITION" ]; then
    echo "ERROR: -partition flag is required when using -slurm."
    exit 1
fi

WORKDIR=$(realpath "$WORKDIR")

# Set default value for SNAKEFILE if not provided
if [ -z "$SNAKEFILE" ]; then
    SNAKEFILE="$WORKDIR/workflow/Snakefile"
fi

# Set default value for CONFIGFILE if not provided
if [ -z "$CONFIGFILE" ]; then

    CONFIGFILE="$WORKDIR/config/config.yaml"
fi

CONFIGFILE=$(realpath "$CONFIGFILE")
SNAKEFILE=$(realpath "$SNAKEFILE")

lines=(
    'run_trimming_qc: "No"'
    'run_kmer_evaluation: "No"'
    'run_asm: "No"'
    'run_auto_scaffolding: "No"'
)

for line in "${lines[@]}"; do
    if ! grep -qF "$line" "$CONFIGFILE"; then
        echo "$line" >> "$CONFIGFILE"
    fi
done

if [ "$FULL" = true ]; then
    sed -i 's/run_trimming_qc: "No"/run_trimming_qc: "Yes"/g' "$CONFIGFILE"
    sed -i 's/run_kmer_evaluation: "No"/run_kmer_evaluation: "Yes"/g' "$CONFIGFILE"
    sed -i 's/run_asm: "No"/run_asm: "Yes"/g' "$CONFIGFILE"
    sed -i 's/run_auto_scaffolding: "No"/run_auto_scaffolding: "Yes"/g' "$CONFIGFILE"
fi

if [ "$TRIMMING_QC" = true ]; then
    sed -i 's/run_trimming_qc: "No"/run_trimming_qc: "Yes"/g' "$CONFIGFILE"
fi

if [ "$KMER_EVAL" = true ]; then
    sed -i 's/run_kmer_evaluation: "No"/run_kmer_evaluation: "Yes"/g' "$CONFIGFILE"
fi

if [ "$ASSEMBLY" = true ]; then
    sed -i 's/run_asm: "No"/run_asm: "Yes"/g' "$CONFIGFILE"
fi

if [ "$SCAFFOLDING" = true ]; then
    sed -i 's/run_auto_scaffolding: "No"/run_auto_scaffolding: "Yes"/g' "$CONFIGFILE"
fi

echo "INFO: Arguments:"
echo "      WORKDIR $WORKDIR"
echo "      CONFIGFILE $CONFIGFILE"
echo "      SNAKEFILE $SNAKEFILE"
echo "      THREADS $THREADS"

echo "INFO: Creating temp dirs."

# Create the singularity cache and temporary files directories
mkdir -p $WORKDIR/singularity $WORKDIR/tmp

# Run script to fetch the Singularity images
python $WORKDIR/workflow/scripts/singularity.py --config $CONFIGFILE

echo "INFO: Running Pipeasm."

if [ "$SETSLURM" = true ]; then
    # Run the Pipeline
    export SINGULARITY_CACHEDIR=$WORKDIR/singularity && \
    export SINGULARITY_TMPDIR=$WORKDIR/tmp && \
    export TMPDIR=$WORKDIR/tmp && \
    sed "s|{WORKDIR}|$WORKDIR|g; s|{THREADS}|$THREADS|g; s|{PARTITION}|$PARTITION|g" \
        $WORKDIR/config/slurm_params.yaml > $WORKDIR/profiles/slurm/config.yaml && \
    snakemake -d $WORKDIR --snakefile $SNAKEFILE \
        --configfile $CONFIGFILE \
        --profile $WORKDIR/profiles/slurm/ \
        $SETUNLOCK $SETNP
else
    # Run the Pipeline
    export SINGULARITY_CACHEDIR=$WORKDIR/singularity && \
    export SINGULARITY_TMPDIR=$WORKDIR/tmp && \
    export TMPDIR=$WORKDIR/tmp && \
    sed "s|{WORKDIR}|$WORKDIR|g; s|{THREADS}|$THREADS|g" \
        $WORKDIR/config/local_params.yaml > $WORKDIR/profiles/local/config.yaml && \
    snakemake -d $WORKDIR --snakefile $SNAKEFILE \
        --configfile $CONFIGFILE \
        --profile $WORKDIR/profiles/local/ \
        --cores $THREADS $SETUNLOCK $SETNP
fi
