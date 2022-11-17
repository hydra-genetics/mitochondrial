__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2022, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule haplocheck:
    input:
        vcf="mitochondrial/gatk_select_variants/{sample}_{type}.vcf",
    output:
        haplocheck_report="mitochondrial/haplocheck/{sample}_{type}.contamination.txt",
        haplocheck_raw_report="mitochondrial/haplocheck/{sample}_{type}.contamination.raw.txt",
        haplocheck_html_report="mitochondrial/haplocheck/{sample}_{type}.contamination.html",
    params:
        extra=config.get("haplocheck", {}).get("extra", ""),
    log:
        "mitochondrial/haplocheck/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "mitochondrial/haplocheck/{sample}_{type}.output.benchmark.tsv",
            config.get("haplocheck", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("haplocheck", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("haplocheck", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("haplocheck", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("haplocheck", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("haplocheck", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("haplocheck", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("haplocheck", {}).get("container", config["default_container"])
    conda:
        "../envs/haplocheck.yaml"
    message:
        "{rule}: Estimate contamination from {input.vcf} using haplocheck"
    shell:
        "haplocheck --raw --out {output.haplocheck_report} {input.vcf} &> {log}"



