#!/bin/bash
#SBATCH --job-name=pipeasm        # Job name
#SBATCH --output=output.txt           # Standard output file
#SBATCH --error=error.txt             # Standard error file
#SBATCH --partition=op1_slurm    # Partition or queue name
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks-per-node=1           # Number of tasks per node
#SBATCH --mem=2G
#SBATCH --cpus-per-task=4

printf "Begin: `date`\n";
TBEGIN=`echo "print time();" | perl`

printf "\n"
# Activate conda envs

source /mnt/gen/opt/softwares/miniconda3/etc/profile.d/conda.sh
conda activate /mnt/gen/opt/softwares/miniforge3/envs/organpipe

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
bash Pipeasm.sh -d . -c config/config.yaml -s workflow/Snakefile -t <threads> -slurm


TEND=`echo "print time();" | perl`

printf "\n"
printf "End: `date`\n";
printf "Run time (s): `expr $TEND - $TBEGIN`\n";
printf "Run time (min): `expr $(( ($TEND - $TBEGIN)/60 ))`\n";
