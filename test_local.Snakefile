import os

rule all:
    input: "/tmp/pipe_sf_test.txt"

rule test_path:
    output: "/tmp/pipe_sf_test.txt"
    params:
        wsf = workflow.snakefile,
        pdir = os.path.dirname(os.path.abspath(workflow.snakefile))
    shell:
        """
        echo "wsf={params.wsf}" > {output}
        echo "pdir={params.pdir}" >> {output}
        """
