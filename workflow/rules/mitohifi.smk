__author__ = "Magdalena Z"
__copyright__ = "Copyright 2024, Uppsala Universitet"
__email__ = "magdalena.z@scilifelab.uu.se"
__license__ = "GPL-3"


'''
# Rule for running Parabricks DeepVariant
rule run_deepvariant:
    input:
        bam=config['parabricks']['input_bam'],
        ref=config['parabricks']['reference']
    output:
        vcf=config['parabricks']['output_vcf']
    params:
        chunk_size=config['parabricks']['chunk_size'],
        num_gpus=config['parabricks']['num_gpus'],
        num_workers=config['parabricks']['num_workers']
    shell:
        """
        parabricks deepvariant \
          --input-bam {input.bam} \
          --output-vcf {output.vcf} \
          --reference {input.ref} \
          --chunk-size {params.chunk_size} \
          --num-gpus {params.num_gpus} \
          --num-workers {params.num_workers}
        """
'''


# Rule for running MitoHiFi with Docker
rule run_mitohifi:
    input:
        fasta=config.get("mitohifi", {}).get("fasta", ""),  # Replace with actual input
    output:
        outfile="snv_indels/mitohifi/{sample}_{type}.mitohifi.res",  # This file doesn't actually exist, I just put it for the workflow to work, the output needed by the program is the outdir
    params:
        genome=config.get("mitohifi", {}).get("genome", ""),
        extra=config.get("mitohifi", {}).get("extra", ""),
        outdir=config["mitohifi"]["outdir"],  # Replace with actual output directory
        #outfolder=config.get("mitohifi", {}).get("outfolder", ""),
    log:
        "snv_indels/mitohifi/{sample}_{type}.mitohifi.log",
    benchmark:
        repeat(
            "snv_indels/mitohifi/{sample}_{type}.out.benchmark.tsv",
            config.get("mitohifi", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("mitohifi", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("mitohifi", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mitohifi", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("mitohifi", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("mitohifi", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("mitohifi", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("mitohifi", {}).get("container", config["default_container"])
    message:
        "{rule}: Calls SNVs on {input.fasta} using mitohifi"
    shell:
        """
        mkdir -p {output.outdir}
        mitohifi \
          -i {input.fasta} \
          -o {output.outdir}  >& {log}
        """


# Use ruleorder if needed to prioritize certain rules
# ruleorder: run_deepvariant > run_mitohifi
