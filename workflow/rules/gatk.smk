__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2022, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"



## Subset to ChrM reads and Create unmapped Bam file

rule gatk_print_reads:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        bam="mito_snv_indels/gatk_print_reads/{sample}_{type}.bam",
    params:
        extra=config.get("gatk_print_reads", {}).get("extra", ""),
        interval=config.get("gatk_print_reads", {}).get("interval", ""),
    log:
        "mito_snv_indels/gatk_print_reads/{sample}_{type}.bam.log",
    benchmark:
        repeat(
            "mito_snv_indels/gatk_print_reads/{sample}_{type}.output.benchmark.tsv",
            config.get("gatk_print_reads", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("gatk_print_reads", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_print_reads", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_print_reads", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_print_reads", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_print_reads", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_print_reads", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_print_reads", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: Extract {params.interval} reads from {input.bam}"
    shell:
        "gatk PrintReads "
        "-I {input.bam} "
        "-R {input.ref} "
        "-L {params.interval} "
        "-O {output.bam} "
        "--read-filter MateOnSameContigOrNoMappedMateReadFilter "
        "--read-filter MateUnmappedAndUnmappedReadFilter "
        "--read-index {input.bai} "
        

rule gatk_revert_sam:
    input:
        bam="mito_snv_indels/gatk_print_reads/{sample}_{type}.bam",
    output:
        bam="mito_snv_indels/gatk_revert_sam/{sample}_{type}.bam",
    params:
        extra=config.get("gatk_revert_sam", {}).get("extra", ""),
    log:
        "mito_snv_indels/gatk_revert_sam/{sample}_{type}.bam.log",
    benchmark:
        repeat(
            "mito_snv_indels/gatk_revert_sam/{sample}_{type}.output.benchmark.tsv",
            config.get("gatk_revert_sam", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("gatk_revert_sam", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_revert_sam", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_revert_sam", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_revert_sam", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_revert_sam", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_revert_sam", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_revert_sam", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: Revert {input.bam} to unmapped BAM"
    shell:
        "gatk RevertSam "
        "--INPUT {input.bam} "
        "--OUTPUT_BY_READGROUP false "
        "--OUTPUT {output.bam} "
        "--VALIDATION_STRINGENCY LENIENT "
        "--ATTRIBUTE_TO_CLEAR FT "
        "--ATTRIBUTE_TO_CLEAR CO "
        "--SORT_ORDER queryname "
        "--RESTORE_ORIGINAL_QUALITIES false "






