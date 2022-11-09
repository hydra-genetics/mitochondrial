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
        bam="mitochondrial/gatk_print_reads/{sample}_{type}.bam",
    params:
        extra=config.get("gatk_print_reads", {}).get("extra", ""),
        interval=config.get("gatk_print_reads", {}).get("interval", ""),
    log:
        "mitochondrial/gatk_print_reads/{sample}_{type}.bam.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_print_reads/{sample}_{type}.bam.benchmark.tsv",
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
        "{rule}: Extract {params.interval} reads from {input.bam} using gatk PrintReads"
    shell:
        "(gatk --java-options '-Xmx3g' PrintReads "
        "-I {input.bam} "
        "-R {input.ref} "
        "-L {params.interval} "
        "-O {output.bam} "
        "--read-filter MateOnSameContigOrNoMappedMateReadFilter "
        "--read-filter MateUnmappedAndUnmappedReadFilter "
        "--read-index {input.bai}) &>{log}"
        

rule gatk_revert_sam:
    input:
        bam="mitochondrial/gatk_print_reads/{sample}_{type}.bam",
    output:
        bam="mitochondrial/gatk_revert_sam/{sample}_{type}.bam",
    params:
        extra=config.get("gatk_revert_sam", {}).get("extra", ""),
    log:
        "mitochondrial/gatk_revert_sam/{sample}_{type}.bam.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_revert_sam/{sample}_{type}.bam.benchmark.tsv",
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
        "{rule}: Revert {input.bam} to unmapped BAM using using the gatk4 wrapper of Picard's RevertSam "
    shell:
        "(gatk --java-options '-Xmx3g' RevertSam "
        "--INPUT {input.bam} "
        "--OUTPUT_BY_READGROUP false "
        "--OUTPUT {output.bam} "
        "--VALIDATION_STRINGENCY LENIENT "
        "--ATTRIBUTE_TO_CLEAR FT "
        "--ATTRIBUTE_TO_CLEAR CO "
        "--SORT_ORDER queryname "
        "--RESTORE_ORIGINAL_QUALITIES false) &>{log} "


rule gatk_sam_to_fastq:
    input:
        bam="mitochondrial/gatk_revert_sam/{sample}_{type}.bam",
    output:
        fq1="mitochondrial/gatk_sam_to_fastq/{sample}_{type}_1.fastq",
        fq2="mitochondrial/gatk_sam_to_fastq/{sample}_{type}_2.fastq",
    params:
        extra=config.get("gatk_sam_to_fastq", {}).get("extra", ""),
    log:
        "mitochondrial/gatk_sam_to_fastq/{sample}_{type}.fastq.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_sam_to_fastq/{sample}_{type}.fastq.benchmark.tsv",
            config.get("gatk_sam_to_fastq", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("gatk_sam_to_fastq", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_sam_to_fastq", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_sam_to_fastq", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_sam_to_fastq", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_sam_to_fastq", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_sam_to_fastq", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_sam_to_fastq", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: Convert {input.bam} to fastq using the gatk4 wrapper of Picard's SamtoFastq"
    shell:
        "gatk --java-options '-Xmx5g' SamToFastq "
        "-INPUT {input.bam} "
        "-F {output.fq1} "
        "-F2 {output.fq2} "
        "-NON_PF true "


# Merge the unmapped BAM with to the aligned BAM (separate merging for BAM mapped to mt referenece and  mt_shifted )


rule gatk_merge_bam_alignment:
    input:
        bam="mitochondrial/bwa_mem/{sample}_{type}_{mt_ref}.bam",
        ubam="mitochondrial/gatk_revert_sam/{sample}_{type}.bam",
        ref=lambda wildcards: config.get("mt_reference", {}).get(wildcards.mt_ref, "")
    output:
        bam="mitochondrial/gatk_merge_bam_alignment/{sample}_{type}_{mt_ref}.bam",
    params:
        bwa_commandline=lambda wildcards, input: get_pg_info(input.bam)[1],
        bwa_version=lambda wildcards, input: get_pg_info(input.bam)[0],
        extra=config.get("gatk_merge_bam_alignment", {}).get("extra", ""),
    log:
        "mitochondrial/gatk_merge_bam_alignment/{sample}_{type}_{mt_ref}.bam.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_merge_bam_alignment/{sample}_{type}_{mt_ref}.bam.benchmark.tsv",
            config.get("gatk_merge_bam_alignment", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("gatk_merge_bam_alignment", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_merge_bam_alignment", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_merge_bam_alignment", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_merge_bam_alignment", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_merge_bam_alignment", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_merge_bam_alignment", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_merge_bam_alignment", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: Merge the ALigned {input.bam} and unmapped {input.ubam} using the gatk4 wrapper of Picard's MergeBamAlignment"
    shell:
        "(gatk --java-options '-Xmx4g' MergeBamAlignment  "
        "-VALIDATION_STRINGENCY SILENT "
        "-EXPECTED_ORIENTATIONS FR "
        "-ATTRIBUTES_TO_RETAIN X0 "
        "-ATTRIBUTES_TO_REMOVE NM "
        "-ATTRIBUTES_TO_REMOVE MD "
        "-ALIGNED_BAM {input.bam} "
        "-UNMAPPED_BAM {input.ubam} "
        "-OUTPUT {output.bam} "
        "-REFERENCE_SEQUENCE {input.ref} "
        "-PAIRED_RUN true "
        "-SORT_ORDER unsorted "
        "-IS_BISULFITE_SEQUENCE false "
        "-ALIGNED_READS_ONLY false "
        "-CLIP_ADAPTERS false "
        "-MAX_RECORDS_IN_RAM 2000000 "
        "-ADD_MATE_CIGAR true "
        "-MAX_INSERTIONS_OR_DELETIONS -1 "
        "-PRIMARY_ALIGNMENT_STRATEGY MostDistant "
        "-PROGRAM_RECORD_ID bwamem " 
        "-PROGRAM_GROUP_VERSION '{params.bwa_version}' " 
        "-PROGRAM_GROUP_COMMAND_LINE '{params.bwa_commandline}' " 
        "-PROGRAM_GROUP_NAME bwamem " 
        "-UNMAPPED_READ_STRATEGY COPY_TO_TAG "
        "-ALIGNER_PROPER_PAIR_FLAGS true "
        "-UNMAP_CONTAMINANT_READS true "
        "-ADD_PG_TAG_TO_READS false) &> {log}"


rule gatk_mark_duplicates:
    input:
        bam="mitochondrial/gatk_merge_bam_alignment/{sample}_{type}_{mt_ref}.bam",
    output:
        bam="mitochondrial/gatk_mark_duplicates/{sample}_{type}_{mt_ref}.bam",
        metrics="mitochondrial/gatk_mark_duplicates/{sample}_{type}_{mt_ref}.metrics.txt"
    params:
        extra=config.get("gatk_mark_duplicates", {}).get("extra", ""),
    log:
        "mitochondrial/gatk_mark_duplicates/{sample}_{type}_{mt_ref}.bam.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_mark_duplicates/{sample}_{type}_{mt_ref}.bam.benchmark.tsv",
            config.get("gatk_mark_duplicates", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("gatk_mark_duplicates", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_mark_duplicates", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_mark_duplicates", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_mark_duplicates", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_mark_duplicates", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_mark_duplicates", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_mark_duplicates", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: Mark duplicates in {input} using the gatk4 wrapper of Picard's MarkDuplicates"
    shell:
        "(gatk --java-options '-Xmx4g' MarkDuplicates "
        "-INPUT {input.bam} "
        "-OUTPUT {output.bam} "
        "-METRICS_FILE {output.metrics} "
        "-VALIDATION_STRINGENCY SILENT "
        "-OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 "
        "-ASSUME_SORT_ORDER queryname " 
        "-CLEAR_DT false "
        "-ADD_PG_TAG_TO_READS false) &> {log}"


rule gatk_sort_sam:
    input:
        bam="mitochondrial/gatk_mark_duplicates/{sample}_{type}_{mt_ref}.bam",
    output:
        bam="mitochondrial/gatk_sort_sam/{sample}_{type}_{mt_ref}.bam",
    params:
        extra=config.get("gatk_sort_sam", {}).get("extra", ""),
    log:
        "mitochondrial/gatk_sort_sam/{sample}_{type}_{mt_ref}.bam.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_sort_sam/{sample}_{type}_{mt_ref}.bam.benchmark.tsv",
            config.get("gatk_sort_sam", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("gatk_sort_sam", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_sort_sam", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_sort_sam", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_sort_sam", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_sort_sam", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_sort_sam", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_sort_sam", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: Sort {input.bam} by coordinates using the gatk4 wrapper of Picard's SortSam"
    shell:
        "(gatk --java-options '-Xmx4g' SortSam  "
        "-INPUT {input.bam} "
        "-OUTPUT {output.bam} "
        "-SORT_ORDER coordinate " 
        "-CREATE_INDEX true "
        "-MAX_RECORDS_IN_RAM 300000) &> {log}"



rule gatk_collect_wgs_metrics:
    input:
        bam="mitochondrial/gatk_sort_sam/{sample}_{type}_{mt_ref}.bam",
        ref=lambda wildcards: config.get("mt_reference", {}).get(wildcards.mt_ref, "")
    output:
        metrics="mitochondrial/gatk_collect_wgs_metrics/{sample}_{type}_{mt_ref}.metrics.txt",
        t_sensitivity="mitochondrial/gatk_collect_wgs_metrics/{sample}_{type}_{mt_ref}.theoretical_sensitivity.txt"
    params:
        coverage_cap=config.get("gatk_collect_wgs_metrics").get("coverage_cap", ""),
        extra=config.get("gatk_collect_wgs_metrics", {}).get("extra", ""),
        read_length=config.get("gatk_collect_wgs_metrics").get("read_length"),
    log:
        "mitochondrial/gatk_collect_wgs_metrics/{sample}_{type}_{mt_ref}.metrics.txt.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_collect_wgs_metrics/{sample}_{type}_{mt_ref}.metrics.txt.benchmark.tsv",
            config.get("gatk_collect_wgs_metrics", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("gatk_collect_wgs_metrics", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_collect_wgs_metrics", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_collect_wgs_metrics", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_collect_wgs_metrics", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_collect_wgs_metrics", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_collect_wgs_metrics", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_collect_wgs_metrics", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: Collect coverage and performance metrics for {input.bam} using the gatk4 wrapper of Picard's CollectWgsMetrics"
    shell:
        "(gatk --java-options '-Xmx2g' CollectWgsMetrics "
        "-INPUT {input.bam} "
        "-VALIDATION_STRINGENCY SILENT "
        "-REFERENCE_SEQUENCE {input.ref} "
        "-OUTPUT {output.metrics} "
        "-USE_FAST_ALGORITHM true "
        "-READ_LENGTH {params.read_length} "
        "-COVERAGE_CAP  {params.coverage_cap} "
        "-INCLUDE_BQ_HISTOGRAM true "
        "-THEORETICAL_SENSITIVITY_OUTPUT {output.t_sensitivity}) &> {log}"
        


rule gatk_extract_average_coverage:
    input:
        metrics="mitochondrial/gatk_collect_wgs_metrics/{sample}_{type}_{mt_ref}.metrics.txt",
    output:
        mean="mitochondrial/gatk_collect_wgs_metrics/{sample}_{type}_{mt_ref}.mean_coverage.txt",
        median="mitochondrial/gatk_collect_wgs_metrics/{sample}_{type}_{mt_ref}.median_coverage.txt"
    params:
        extra=config.get("gatk_collect_wgs_metrics", {}).get("extra", ""),
    log:
        "mitochondrial/gatk_collect_wgs_metrics/{sample}_{type}_{mt_ref}.coverage.txt.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_collect_wgs_metrics/{sample}_{type}_{mt_ref}.output.benchmark.tsv",
            config.get("gatk_collect_wgs_metrics", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("gatk_collect_wgs_metrics", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_collect_wgs_metrics", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_collect_wgs_metrics", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_collect_wgs_metrics", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_collect_wgs_metrics", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_collect_wgs_metrics", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_extract_average_coverage", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: Extract the mean and median coverage from {input.metrics}"
    script:
        "../scripts/gatk_extract_average_coverage.py"
