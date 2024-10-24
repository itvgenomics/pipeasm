rule download_buscodb:
    output:
        directory("resources/{sample}_buscodb")
    threads: 1
    params:
        buscodb= config['buscodb']
    log:
        "logs/{sample}.download_buscodb.log"
    singularity:
        "docker://huangnengcsu/compleasm:v0.2.5"
    shell:
        """
        bash -c 'mkdir {output} && \
        compleasm download {params} --library_path {output} 2> {log}'
        """

# rule remove_buscodb:
#     input:
#         directory("resources/{sample}_buscodb"),
#         files_remove_buscodb()
#     output:
#         "logs/{sample}_buscodb.check"
#     threads: 1
#     shell:
#         "rm -r {input[0]} && touch {output}"
