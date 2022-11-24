__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2022, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"



## Subset to ChrM reads and Create unmapped Bam file

rule gatk_print_reads:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
    output:
        bam=temp("mitochondrial/gatk_print_reads/{sample}_{type}.bam"),
        bai=temp("mitochondrial/gatk_print_reads/{sample}_{type}.bai")
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
        "-L {params.interval} "
        "-O {output.bam} "
        "--read-filter MateOnSameContigOrNoMappedMateReadFilter "
        "--read-filter MateUnmappedAndUnmappedReadFilter "
        "--read-index {input.bai}) &>{log}"
        

rule gatk_revert_sam:
    input:
        bam="mitochondrial/gatk_print_reads/{sample}_{type}.bam",
    output:
        bam=temp("mitochondrial/gatk_revert_sam/{sample}_{type}.bam"),
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
        fq1=temp("mitochondrial/gatk_sam_to_fastq/{sample}_{type}_1.fastq"),
        fq2=temp("mitochondrial/gatk_sam_to_fastq/{sample}_{type}_2.fastq"),
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
        "(gatk --java-options '-Xmx5g' SamToFastq "
        "-INPUT {input.bam} "
        "-F {output.fq1} "
        "-F2 {output.fq2} "
        "-NON_PF true) &> {log} "


# Merge the unmapped BAM with to the aligned BAM (separate merging for BAM mapped to mt referenece and  mt_shifted )


rule gatk_merge_bam_alignment:
    input:
        bam="mitochondrial/bwa_mem_mito/{sample}_{type}_{mt_ref}.bam",
        ubam="mitochondrial/gatk_revert_sam/{sample}_{type}.bam",
        ref=lambda wildcards: config.get("mt_reference", {}).get(wildcards.mt_ref, "")
    output:
        bam=temp("mitochondrial/gatk_merge_bam_alignment/{sample}_{type}_{mt_ref}.bam"),
    params:
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
        """
        version=`(gatk PrintReadsHeader -I {input.bam} --verbosity 'ERROR'  \
        -O /dev/stdout | grep 'PN:bwa' | cut -f4 | cut -f2 -d':')`

        command=`gatk PrintReadsHeader -I {input.bam}  --verbosity 'ERROR' \
        -O /dev/stdout | grep 'PN:bwa' | cut -f5 | cut -f2- -d':'`
        
        (gatk --java-options '-Xmx4g' MergeBamAlignment \
        -VALIDATION_STRINGENCY SILENT \
        -EXPECTED_ORIENTATIONS FR \
        -ATTRIBUTES_TO_RETAIN X0 \
        -ATTRIBUTES_TO_REMOVE NM \
        -ATTRIBUTES_TO_REMOVE MD \
        -ALIGNED_BAM {input.bam} \
        -UNMAPPED_BAM {input.ubam} \
        -OUTPUT {output.bam} \
        -REFERENCE_SEQUENCE {input.ref} \
        -PAIRED_RUN true \
        -SORT_ORDER unsorted \
        -IS_BISULFITE_SEQUENCE false \
        -ALIGNED_READS_ONLY false \
        -CLIP_ADAPTERS false \
        -MAX_RECORDS_IN_RAM 2000000 \
        -ADD_MATE_CIGAR true \
        -MAX_INSERTIONS_OR_DELETIONS -1 \
        -PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        -PROGRAM_RECORD_ID "bwamem" \
        -PROGRAM_GROUP_VERSION "$version" \
        -PROGRAM_GROUP_COMMAND_LINE "$command" \
        -PROGRAM_GROUP_NAME "bwamem" \
        -UNMAPPED_READ_STRATEGY COPY_TO_TAG \
        -ALIGNER_PROPER_PAIR_FLAGS true \
        -UNMAP_CONTAMINANT_READS true \
        -ADD_PG_TAG_TO_READS false) &> {log}

        """


rule gatk_mark_duplicates:
    input:
        bam="mitochondrial/gatk_merge_bam_alignment/{sample}_{type}_{mt_ref}.bam",
    output:
        bam=temp("mitochondrial/gatk_mark_duplicates/{sample}_{type}_{mt_ref}.bam"),
        metrics=temp("mitochondrial/gatk_mark_duplicates/{sample}_{type}_{mt_ref}.metrics.txt")
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
        bam=temp("mitochondrial/gatk_sort_sam/{sample}_{type}_{mt_ref}.bam"),
        bai=temp("mitochondrial/gatk_sort_sam/{sample}_{type}_{mt_ref}.bai"),
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

## Collect coverage metrics for MT 

rule gatk_collect_wgs_metrics:
    input:
        bam="mitochondrial/gatk_sort_sam/{sample}_{type}_{mt_ref}.bam",
        ref=lambda wildcards: config.get("mt_reference", {}).get(wildcards.mt_ref, "")
    output:
        metrics=temp("mitochondrial/gatk_collect_wgs_metrics/{sample}_{type}_{mt_ref}.metrics.txt"),
        t_sensitivity=temp("mitochondrial/gatk_collect_wgs_metrics/{sample}_{type}_{mt_ref}.theoretical_sensitivity.txt")
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


# Call MT variants with mutect2

rule gatk_mutect2:
    input:
        bam="mitochondrial/gatk_sort_sam/{sample}_{type}_{mt_ref}.bam",
        bai="mitochondrial/gatk_sort_sam/{sample}_{type}_{mt_ref}.bai",
        ref=lambda wildcards: config.get("mt_reference", {}).get(wildcards.mt_ref, ""),
    output:
        stats = temp("mitochondrial/gatk_mutect2/{sample}_{type}_{mt_ref}.vcf.stats"),
        vcf=temp("mitochondrial/gatk_mutect2/{sample}_{type}_{mt_ref}.vcf"),
        idx=temp("mitochondrial/gatk_mutect2/{sample}_{type}_{mt_ref}.vcf.idx"),
    params:
        extra=config.get("gatk_mutect2", {}).get("extra", ""),
        interval=lambda wildcards: config.get("gatk_mutect2", {}).get("interval", {}).get(wildcards.mt_ref, ""),
        max_reads_per_alignment_start=config.get("gatk_mutect2", {}).get("max_reads_per_alignment_start", 75),
    log:
        "mitochondrial/gatk_mutect2/{sample}_{type}_{mt_ref}.vcf.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_mutect2/{sample}_{type}_{mt_ref}.vcf.benchmark.tsv",
            config.get("gatk_mutect2", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("gatk_mutect2", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_mutect2", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_mutect2", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_mutect2", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_mutect2", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_mutect2", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_mutect2", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: Call mitochondrial variants using Mutect2 on {input.bam}"
    shell:
        "(gatk --java-options '-Xmx3g' Mutect2 "
        "-R {input.ref} "
        "-I {input.bam} "
        "--read-filter MateOnSameContigOrNoMappedMateReadFilter "
        "--read-filter MateUnmappedAndUnmappedReadFilter "
        "-O {output.vcf} "
        "-L {params.interval} "
        "--annotation StrandBiasBySample "
        "--mitochondria-mode "
        "--max-reads-per-alignment-start {params.max_reads_per_alignment_start} "
        "--max-mnp-distance 0) &> {log}"

# Liftover the variant from the control region in the shifted VCF to the mt reference genome coordinaged 
# and merge with the non-control region chrM VCF

rule gatk_lift_over_vcf:
    input:
        shifted_vcf="mitochondrial/gatk_mutect2/{sample}_{type}_mt_shifted.vcf",
        ref=config.get("mt_reference", {}).get("mt", ""),
        shift_back_chain=config.get("gatk_lift_over_vcf", {}).get("shift_back_chain", ""),
    output:
        shifted_back_vcf=temp("mitochondrial/gatk_lift_over_vcf/{sample}_{type}.shifted_back.vcf"),
        shifted_back_idx=temp("mitochondrial/gatk_lift_over_vcf/{sample}_{type}.shifted_back.vcf.idx"),
        reject=temp("mitochondrial/gatk_lift_over_vcf/{sample}_{type}.reject.vcf"),
    params:
        extra=config.get("gatk_lift_over_vcf", {}).get("extra", ""),
    log:
        "mitochondrial/gatk_lift_over_vcf/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_lift_over_vcf/{sample}_{type}.vcf.benchmark.tsv",
            config.get("gatk_lift_over_vcf", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("gatk_lift_over_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_lift_over_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_lift_over_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_lift_over_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_lift_over_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_lift_over_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_lift_over_vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: Lift over shifted {input.shifted_vcf} using chain file to the original unshifted fasta reference {input.ref}"
    shell:
        "(gatk --java-options '-Xmx3g' LiftoverVcf "
        "-I {input.shifted_vcf} "
        "-O {output.shifted_back_vcf} "
        "-R {input.ref} "
        "--CHAIN {input.shift_back_chain} "
        "--REJECT {output.reject}) &> {log} "


rule gatk_merge_vcfs:
    input:
        shifted_vcf="mitochondrial/gatk_lift_over_vcf/{sample}_{type}.shifted_back.vcf",
        vcf="mitochondrial/gatk_mutect2/{sample}_{type}_mt.vcf",
    output:
        vcf=temp("mitochondrial/gatk_merge_vcfs/{sample}_{type}.vcf"),
        idx=temp("mitochondrial/gatk_merge_vcfs/{sample}_{type}.vcf.idx"),
    params:
        extra=config.get("gatk_merge_vcfs", {}).get("extra", ""),
    log:
        "mitochondrial/gatk_merge_vcfs/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_merge_vcfs/{sample}_{type}.vcf.benchmark.tsv",
            config.get("gatk_merge_vcfs", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("gatk_merge_vcfs", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_merge_vcfs", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_merge_vcfs", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_merge_vcfs", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_merge_vcfs", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_merge_vcfs", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_merge_vcfs", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        """"{rule}: Combine the VCF for the shifted back control region {input.shifted_vcf}
        with the VCF for the rest of chrM {input.vcf}"""
    shell:
        "(gatk --java-options '-Xmx3g' MergeVcfs "
        "-I {input.shifted_vcf} "
        "-I {input.vcf} "
        "-O {output.vcf}) &> {log}"


rule gatk_merge_stats:
    input:
        stats="mitochondrial/gatk_mutect2/{sample}_{type}_mt.vcf.stats",
        shifted_stats="mitochondrial/gatk_mutect2/{sample}_{type}_mt_shifted.vcf.stats",
    output:
        merged_stats=temp("mitochondrial/gatk_merge_stats/{sample}_{type}.vcf.stats"),
    params:
        extra=config.get("gatk_merge_stats", {}).get("extra", ""),
    log:
        "mitochondrial/gatk_merge_stats/{sample}_{type}.vcf.stats.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_merge_stats/{sample}_{type}.vcf.stats.benchmark.tsv",
            config.get("gatk_merge_stats", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("gatk_merge_stats", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_merge_stats", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_merge_stats", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_merge_stats", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_merge_stats", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_merge_stats", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_merge_stats", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: Merge thg mutect2 vcf stats files {input.shifted_stats} and {input.stats}"
    shell:
        "(gatk --java-options '-Xmx3g' MergeMutectStats " 
        "--stats {input.shifted_stats} "
        "--stats {input.stats} "
        "-O {output.merged_stats}) &> {log}"


rule gatk_filter_mutect_calls_mt:
    input:
        vcf="mitochondrial/gatk_merge_vcfs/{sample}_{type}.vcf",
        ref=config.get("mt_reference", {}).get("mt", ""),
        mutect_stats="mitochondrial/gatk_merge_stats/{sample}_{type}.vcf.stats",
    output:
        vcf=temp("mitochondrial/gatk_filter_mutect_calls_mt/{sample}_{type}.vcf"),
        idx=temp("mitochondrial/gatk_filter_mutect_calls_mt/{sample}_{type}.vcf.idx"),
        stats=temp("mitochondrial/gatk_filter_mutect_calls_mt/{sample}_{type}.vcf.filteringStats.tsv"),
    params:
        extra=config.get("gatk_filter_mutect_calls_mt", {}).get("extra", ""),
        max_alt_allele_count=config.get("gatk_filter_mutect_calls_mt", {}).get("max_alt_allele_count", 4),
        vaf_filter_threshold=config.get("gatk_filter_mutect_calls_mt", {}).get("vaf_filter_threshold", 0),
        f_score_beta=config.get("gatk_filter_mutect_calls_mt", {}).get("f_score_beta", 1.0),
        contamination=0.0,
    log:
        "mitochondrial/gatk_filter_mutect_calls_mt/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_filter_mutect_calls_mt/{sample}_{type}.vcf.benchmark.tsv",
            config.get("gatk_filter_mutect_calls_mt", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("gatk_filter_mutect_calls_mt", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_filter_mutect_calls_mt", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_filter_mutect_calls_mt", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_filter_mutect_calls_mt", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_filter_mutect_calls_mt", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_filter_mutect_calls_mt", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_filter_mutect_calls_mt", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: Filter {input.vcf} using FilterMutectCalls in mitochondria mode"
    shell:
        "(gatk --java-options '-Xmx3g' FilterMutectCalls  "
        "-V {input.vcf}  "
        "-R {input.ref}  "
        "-O {output.vcf}  "
        "--stats {input.mutect_stats}  "
        "{params.extra}  "
        "--mitochondria-mode  "
        "--max-alt-allele-count {params.max_alt_allele_count}  "
        "--min-allele-fraction {params.vaf_filter_threshold}  "
        "--f-score-beta {params.f_score_beta} "
        "--contamination-estimate {params.contamination}) &> {log}"


rule gatk_variant_filtration:
    input:
        vcf="mitochondrial/gatk_filter_mutect_calls_mt/{sample}_{type}.vcf",
        blacklisted_sites=config.get("mt_reference", {}).get("blacklist", "")
    output:
        filtered_vcf=temp("mitochondrial/gatk_variant_filtration/{sample}_{type}.vcf"),
        idx=temp("mitochondrial/gatk_variant_filtration/{sample}_{type}.vcf.idx"),
    params:
        extra=config.get("gatk_variant_filtration", {}).get("extra", ""),
    log:
        "mitochondrial/gatk_variant_filtration/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_variant_filtration/{sample}_{type}.vcf.benchmark.tsv",
            config.get("gatk_variant_filtration", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("gatk_variant_filtration", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_variant_filtration", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_variant_filtration", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_variant_filtration", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_variant_filtration", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_variant_filtration", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_variant_filtration", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: Mask blacklisted ChrM sites in {input.vcf} using GATK's VariantFiltration"
    shell:
        "(gatk --java-options '-Xmx3g' VariantFiltration  "
        "-V {input.vcf}  "
        "-O {output.filtered_vcf}  "
        "--apply-allele-specific-filters  "
        "--mask {input.blacklisted_sites}  "
        "--mask-name 'blacklisted_site') &> {log}"


rule gatk_left_align_and_trim_variants:
    input:
        vcf="mitochondrial/gatk_variant_filtration/{sample}_{type}.vcf",
        ref=config.get("mt_reference", {}).get("mt", ""),
    output:
        vcf=temp("mitochondrial/gatk_left_align_and_trim_variants/{sample}_{type}.vcf"),
        idx=temp("mitochondrial/gatk_left_align_and_trim_variants/{sample}_{type}.vcf.idx"),
    params:
        extra=config.get("gatk_left_align_and_trim_variants", {}).get("extra", ""),
    log:
        "mitochondrial/gatk_left_align_and_trim_variants/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_left_align_and_trim_variants/{sample}_{type}.vcf.benchmark.tsv",
            config.get("gatk_left_align_and_trim_variants", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("gatk_left_align_and_trim_variants", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_left_align_and_trim_variants", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_left_align_and_trim_variants", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_left_align_and_trim_variants", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_left_align_and_trim_variants", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_left_align_and_trim_variants", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_left_align_and_trim_variants", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: Left align alleles and split multiallelic sites in {input.vcf}"
    shell:
        "(gatk --java-options '-Xmx3g' LeftAlignAndTrimVariants "
        "-R {input.ref} "
        "-V {input.vcf} "
        "-O {output.vcf} "
        "--split-multi-allelics "
        "--dont-trim-alleles "
        "--keep-original-ac) &> {log}"


rule gatk_select_variants:
    input:
        vcf="mitochondrial/gatk_left_align_and_trim_variants/{sample}_{type}.vcf",
    output:
        vcf=temp("mitochondrial/gatk_select_variants/{sample}_{type}.vcf"),
        idx=temp("mitochondrial/gatk_select_variants/{sample}_{type}.vcf.idx"),
    params:
        extra=config.get("gatk_select_variants", {}).get("extra", ""),
    log:
        "mitochondrial/gatk_select_variants/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_select_variants/{sample}_{type}.vcf.benchmark.tsv",
            config.get("gatk_select_variants", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("gatk_select_variants", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_select_variants", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_select_variants", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_select_variants", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_select_variants", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_select_variants", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_select_variants", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: Exclude filtered sites in {input.vcf}"
    shell:
        "(gatk --java-options '-Xmx3g' SelectVariants "
        "-V {input.vcf} "
        "-O {output.vcf} "
        "--exclude-filtered) &> {log}"


# Filter comtamination using the gatk_filter_mutect_calls_mt rule again

use rule gatk_filter_mutect_calls_mt as gatk_filter_contamination with:
    input:
        vcf="mitochondrial/gatk_select_variants/{sample}_{type}.vcf",
        ref=config.get("mt_reference", {}).get("mt", ""),
        mutect_stats="mitochondrial/gatk_merge_stats/{sample}_{type}.vcf.stats", 
        contamination="mitochondrial/haplocheck/{sample}_{type}.contamination.txt",
    output:
        vcf=temp("mitochondrial/gatk_filter_contamination/{sample}_{type}.vcf"),
        idx=temp("mitochondrial/gatk_filter_contamination/{sample}_{type}.vcf.idx"),
        stats=temp("mitochondrial/gatk_filter_contamination/{sample}_{type}.vcf.filteringStats.tsv")
    params:
        extra=config.get("gatk_filter_mutect_calls_mt", {}).get("extra", ""),
        max_alt_allele_count=config.get("gatk_filter_mutect_calls_mt", {}).get("max_alt_allele_count", 4),
        vaf_filter_threshold=config.get("gatk_filter_mutect_calls_mt", {}).get("vaf_filter_threshold", 0),
        f_score_beta=config.get("gatk_filter_mutect_calls_mt", {}).get("f_score_beta", 1.0),
        contamination=lambda wildcards, input: get_contamination_estimate(wildcards, input.contamination),
    log:
        "mitochondrial/gatk_filter_contamination/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_filter_contamination/{sample}_{type}.vcf.benchmark.tsv",
            config.get("gatk_filter_mutect_calls_mt", {}).get("benchmark_repeats", 1)
        )
    message:
        "{rule}: Filter {input.vcf} using FilterMutectCalls in mitochondria mode with a contamination estimate"


use rule gatk_left_align_and_trim_variants as gatk_split_multi_allelic_sites with:
    input:
        vcf="mitochondrial/gatk_filter_contamination/{sample}_{type}.vcf",
        ref=config.get("mt_reference", {}).get("mt", ""),
    output:
        vcf=temp("mitochondrial/gatk_split_multi_allelic_sites/{sample}_{type}.vcf"),
        idx=temp("mitochondrial/gatk_split_multi_allelic_sites/{sample}_{type}.vcf.idx"),
    log:
        "mitochondrial/gatk_split_multi_allelic_sites/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "mitochondrial/gatk_split_multi_allelic_sites/{sample}_{type}.vcf.benchmark.tsv",
            config.get("gatk_left_align_and_trim_variants", {}).get("benchmark_repeats", 1)
        )
    message:
        "{rule}: Filter {input.vcf} using FilterMutectCalls in mitochondria mode with a contamination estimate"