__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2022, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule bwa_mem:
    input:
        reads=[
            "mitochondrial/gatk_sam_to_fastq/{sample}_{type}_1.fastq",
            "mitochondrial/gatk_sam_to_fastq/{sample}_{type}_2.fastq",    
        ],
        idx=lambda wildcards: [
            config.get("bwa_mem", {}).get("amb", "").get(wildcards.mt_ref, ""),
            config.get("bwa_mem", {}).get("ann", "").get(wildcards.mt_ref, ""),
            config.get("bwa_mem", {}).get("bwt", "").get(wildcards.mt_ref, ""),
            config.get("bwa_mem", {}).get("pac", "").get(wildcards.mt_ref, ""),
            config.get("bwa_mem", {}).get("sa", "").get(wildcards.mt_ref, ""),
        ]

    output:
        bam=temp("mitochondrial/bwa_mem/{sample}_{type}_{mt_ref}.bam"),
    params:
        extra=config.get("bwa_mem", {}).get("extra", ""),
        mt_ref=lambda wildcards: config["mt_reference"][wildcards.mt_ref],
    log:
        "mitochondrial/bwa_mem/{sample}_{type}_{mt_ref}.bam.log",
    benchmark:
        repeat(
            "mitochondrial/bwa_mem/{sample}_{type}_{mt_ref}.output.benchmark.tsv",
            config.get("bwa_mem", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("bwa_mem", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bwa_mem", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bwa_mem", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bwa_mem", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bwa_mem", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bwa_mem", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bwa_mem", {}).get("container", config["default_container"])
    conda:
        "../envs/bwa.yaml"
    message:
        "{rule}: Align chrM reads from {wildcards.sample}_{wildcards.type} to the {wildcards.mt_ref} reference"
    wrapper:
        "v1.3.1/bio/bwa/mem"

