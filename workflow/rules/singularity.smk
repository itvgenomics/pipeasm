import subprocess

def build_fcsgx():
    command = "singularity build resources/fcs-gx.sif docker://ncbi/fcs-gx:latest"
    subprocess.run(command, shell=True, check=True)


rule build_fcsgx_sif:
    output:
        "resources/fcs-gx.sif"
    threads: 1
    run:
        build_fcsgx()
