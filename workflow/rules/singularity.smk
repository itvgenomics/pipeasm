import subprocess

def build_sifs():
    command = "singularity build resources/fcs-gx.sif docker://ncbi/fcs-gx:latest"
    subprocess.run(command, shell=True, check=True)
    command = "singularity build resources/fcs-adaptor.sif docker://ncbi/fcs-adaptor:latest"
    subprocess.run(command, shell=True, check=True)


rule build_sifs:
    output:
        "resources/fcs-gx.sif",
        "resources/fcs-adaptor.sif"
    run:
        build_sifs()