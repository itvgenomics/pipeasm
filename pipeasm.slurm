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
# Activate conda envs

module load snakemake

# Move to workdir
cd $SLURM_SUBMIT_DIR

# Create temp dirs for proper execution
mkdir -p $SLURM_SUBMIT_DIR/singularity $SLURM_SUBMIT_DIR/tmp

### RECOMMENDED ###
# Perform a dry run to see what jobs will be executed without actually running them.
# This is will check all input files, build singularity images and download the seeds. After this command you can run the following, without -np
# bash Pipeasm.sh -d . -c config/config.yaml -s workflow/Snakefile -t <threads> -np

# Caminho dos arquivos do singularity para colocar no config.yaml (sif_dir): "/scratch/sif_files/pipeasm"
# Caminho para colocar no gxdb: "/mnt/gen/opt/utils/lucas/S4-Databases/fcs-gx_testDB/gxdb"
# Example command to run Pipeasm with slurm:
bash Pipeasm.sh -d . -c config/config.yaml -s workflow/Snakefile -j <n_jobs> -slurm -partition op1_slurm


TEND=`echo "print time();" | perl`

printf "\n"
printf "End: `date`\n";
printf "Run time (s): `expr $TEND - $TBEGIN`\n";
printf "Run time (min): `expr $(( ($TEND - $TBEGIN)/60 ))`\n";
