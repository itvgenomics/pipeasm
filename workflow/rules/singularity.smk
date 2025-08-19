import subprocess

def build_fcsgx():
    command = "singularity build resources/fcs-gx.sif docker://ncbi/fcs-gx:latest"
    subprocess.run(command, shell=True, check=True)

def build_fcsadaptor():
    command = "singularity build resources/fcs-adaptor.sif docker://ncbi/fcs-adaptor:latest"
    subprocess.run(command, shell=True, check=True)

rule build_fcsgx_sif:
    output:
        "resources/fcs-gx.sif"
    threads: 1
    run:
        build_fcsgx()

rule build_fcsadaptor_sif:
    output:
        "resources/fcs-adaptor.sif"
    threads: 1
    run:
        build_fcsadaptor()
