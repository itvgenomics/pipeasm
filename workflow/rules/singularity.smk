import subprocess

def build_mitohifi():
    command = "singularity build resources/mitohifi.sif docker://ghcr.io/marcelauliano/mitohifi:master"
    subprocess.run(command, shell=True, check=True)

def build_fcsgx():
    command = "singularity build resources/fcs-gx.sif docker://ncbi/fcs-gx:latest"
    subprocess.run(command, shell=True, check=True)


rule build_mitohifi_sif:
    output:
        "resources/mitohifi.sif"
    threads: 1
    run:
        build_mitohifi()

rule build_fcsgx_sif:
    output:
        "resources/fcs-gx.sif"
    threads: 1
    run:
        build_fcsgx()
