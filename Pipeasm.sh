#!/bin/bash
cat << 'EOF'

### Pipeasm - a tool for automated large genome assembly and analysis 
#Authors: Trindade F., Silva B. M., Canesin L., Souza Junior R. O., Oliveira R. 2024

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
        ;;
    -s)
        shift
        SNAKEFILE=$1
        ;;
    --trimming_qc)
        TRIMMING_QC=true
        FULL=false
        ;;
    --kmer_eval)
        KMER_EVAL=true
        FULL=false
        ;;
    --assembly)
        ASSEMBLY=true
        FULL=false
        ;;
    --scaffolding)
        SCAFFOLDING=true
        FULL=false
        ;;
    *)
        echo "ERROR: Invalid argument definition. Ensure there are no spaces within your strings and that all arguments are correctly specified."
        exit 1
        ;;
    esac
    shift
done

# Ensure THREADS and WORKDIR have values
if [ -z "$WORKDIR" ] || [ -z "$THREADS" ]; then
    cat << 'EOF'
ERROR: Missing one of the required arguments: -d (Work Directory), -t (# Threads)

## Usage: ./Pipeasm.sh -c <config.yaml> -d </path/to/work/dir> -s </path/to/snakefile> -t <# threads>
#-d </path/to/work/dir> = Path to your working directory where all the workflow file are
#-c </path/to/config.yaml> = Overwrite the default configuration file with all nedded parameters (config/config.yaml).
#-s </path/to/snakefile> = Overwrite the default snakefile path (workflow/Snakefile)
#-t <int> = Number of threads to use

# You can choose a Pipeasm step with: 
    --trimming_qc (for Trimming and Quality Control);
    --kmer_eval (for k-mer Evaluation stats/plots);
    --assembly (for all Assembly and Decontamination steps);
    --scaffolding (for all Scaffolding and Hi-C Map steps);

# Only the -d and -t flags are required if you want to use the Snakemake default parameters
    
EOF

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

echo "INFO: Running Pipeasm."

# Run the Pipeline
export SINGULARITY_CACHEDIR=$WORKDIR/singularity && \
export TMPDIR=$WORKDIR/tmp && \
snakemake -d $WORKDIR --snakefile $SNAKEFILE --configfile $CONFIGFILE --use-singularity --singularity-args "-B $WORKDIR:/data --pwd /data" --cores $THREADS
