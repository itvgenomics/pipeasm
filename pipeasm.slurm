#!/bin/bash
#SBATCH --job-name=pipeasm        # Job name
#SBATCH --output=output.txt           # Standard output file
#SBATCH --error=error.txt             # Standard error file
#SBATCH --partition=op1_slurm    # Partition or queue name
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --mem=2G
#SBATCH --cpus-per-task=4
#SBATCH --time=96:00:00                # Time limit hrs:min:sec

printf "Begin: `date`\n";
TBEGIN=`echo "print time();" | perl`

printf "\n"
# Load Snakemake and Singularity (e.g., via environment modules or conda)

module load snakemake

# Move to workdir
cd $SLURM_SUBMIT_DIR

# Create temp dirs for proper execution
mkdir -p $SLURM_SUBMIT_DIR/singularity $SLURM_SUBMIT_DIR/tmp

### RECOMMENDED ###
# Perform a dry run to see what jobs will be executed without actually running them.
# This will check all input files and build the Singularity images. After this command you can run the following, without -np
# bash Pipeasm.sh -d . -c config/config.yaml -s workflow/Snakefile -j <n_jobs> -slurm -partition op1_slurm -np

# Example Singularity image directory to set in config.yaml (sif_dir): "/scratch/sif_files/pipeasm"
# Example FCS-GX database path to set in config.yaml (gxdb): "/mnt/gen/opt/utils/lucas/S4-Databases/fcs-gx_testDB/gxdb"
# Example command to run Pipeasm with slurm:
bash Pipeasm.sh -d . -c config/config.yaml -s workflow/Snakefile -j <n_jobs> -slurm -partition op1_slurm


TEND=`echo "print time();" | perl`

printf "\n"
printf "End: `date`\n";
printf "Run time (s): `expr $TEND - $TBEGIN`\n";
printf "Run time (min): `expr $(( ($TEND - $TBEGIN)/60 ))`\n";
