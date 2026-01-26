"""
Rules for calculating RPKM values from featureCounts outputs.
Generates RPKM files for both all reads and uniquely mapped reads.
"""
import re


def is_valid_paired_run_tag(run_tag):
    if run_tag in SAMPLES_DF.index and run_tag not in RUN_TAGS:
        return is_paired_end(run_tag)
    if run_tag in RUN_TAG_SAMPLES:
        samples = RUN_TAG_SAMPLES[run_tag]
        if samples:
            return is_paired_end(samples[0])
    return False


def is_valid_single_run_tag(run_tag):
    if run_tag in SAMPLES_DF.index and run_tag not in RUN_TAGS:
        return is_single_end(run_tag) or is_nanopore(run_tag)
    if run_tag in RUN_TAG_SAMPLES:
        samples = RUN_TAG_SAMPLES[run_tag]
        if samples:
            return is_single_end(samples[0]) or is_nanopore(samples[0])
    return False


PAIRED_RUN_TAGS = [rt for rt in EFFECTIVE_SAMPLES_TO_PROCESS if is_valid_paired_run_tag(rt)]
SINGLE_RUN_TAGS = [rt for rt in EFFECTIVE_SAMPLES_TO_PROCESS if is_valid_single_run_tag(rt)]


rule rpkm_paired_all:
    wildcard_constraints:
        run_tag = "|".join([re.escape(rt) for rt in PAIRED_RUN_TAGS]) if PAIRED_RUN_TAGS else "^$"
    input:
        counts = get_processing_path("{run_tag}/{run_tag}_counts_paired_all.txt")
    output:
        rpkm = get_processing_path("{run_tag}/{run_tag}_rpkm_paired_all.txt")
    log:
        get_processing_path("{run_tag}/logs/rpkm_paired_all.log")
    conda:
        "../../envs/rpkm.yaml"
    script:
        "rules/scripts/calc_rpkm.py"


rule rpkm_paired_unique:
    wildcard_constraints:
        run_tag = "|".join([re.escape(rt) for rt in PAIRED_RUN_TAGS]) if PAIRED_RUN_TAGS else "^$"
    input:
        counts = get_processing_path("{run_tag}/{run_tag}_counts_paired_unique.txt")
    output:
        rpkm = get_processing_path("{run_tag}/{run_tag}_rpkm_paired_unique.txt")
    log:
        get_processing_path("{run_tag}/logs/rpkm_paired_unique.log")
    conda:
        "../../envs/rpkm.yaml"
    script:
        "rules/scripts/calc_rpkm.py"


rule rpkm_paired_complete:
    wildcard_constraints:
        run_tag = "|".join([re.escape(rt) for rt in PAIRED_RUN_TAGS]) if PAIRED_RUN_TAGS else "^$"
    input:
        rpkm_all = get_processing_path("{run_tag}/{run_tag}_rpkm_paired_all.txt"),
        rpkm_unique = get_processing_path("{run_tag}/{run_tag}_rpkm_paired_unique.txt")
    output:
        flag = get_processing_path("{run_tag}/rpkm_paired_complete.flag")
    log:
        get_processing_path("{run_tag}/logs/rpkm_paired_complete.log")
    shell:
        """
        echo "RPKM calculation completed for paired-end sample {wildcards.run_tag}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "All reads RPKM: {input.rpkm_all}" >> {output.flag}
        echo "Unique reads RPKM: {input.rpkm_unique}" >> {output.flag}
        """


rule rpkm_single_all:
    wildcard_constraints:
        run_tag = "|".join([re.escape(rt) for rt in SINGLE_RUN_TAGS]) if SINGLE_RUN_TAGS else "^$"
    input:
        counts = get_processing_path("{run_tag}/{run_tag}_counts_single_all.txt")
    output:
        rpkm = get_processing_path("{run_tag}/{run_tag}_rpkm_single_all.txt")
    log:
        get_processing_path("{run_tag}/logs/rpkm_single_all.log")
    conda:
        "../../envs/rpkm.yaml"
    script:
        "rules/scripts/calc_rpkm.py"


rule rpkm_single_unique:
    wildcard_constraints:
        run_tag = "|".join([re.escape(rt) for rt in SINGLE_RUN_TAGS]) if SINGLE_RUN_TAGS else "^$"
    input:
        counts = get_processing_path("{run_tag}/{run_tag}_counts_single_unique.txt")
    output:
        rpkm = get_processing_path("{run_tag}/{run_tag}_rpkm_single_unique.txt")
    log:
        get_processing_path("{run_tag}/logs/rpkm_single_unique.log")
    conda:
        "../../envs/rpkm.yaml"
    script:
        "rules/scripts/calc_rpkm.py"


rule rpkm_single_complete:
    wildcard_constraints:
        run_tag = "|".join([re.escape(rt) for rt in SINGLE_RUN_TAGS]) if SINGLE_RUN_TAGS else "^$"
    input:
        rpkm_all = get_processing_path("{run_tag}/{run_tag}_rpkm_single_all.txt"),
        rpkm_unique = get_processing_path("{run_tag}/{run_tag}_rpkm_single_unique.txt")
    output:
        flag = get_processing_path("{run_tag}/rpkm_single_complete.flag")
    log:
        get_processing_path("{run_tag}/logs/rpkm_single_complete.log")
    shell:
        """
        echo "RPKM calculation completed for single-end sample {wildcards.run_tag}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "All reads RPKM: {input.rpkm_all}" >> {output.flag}
        echo "Unique reads RPKM: {input.rpkm_unique}" >> {output.flag}
        """
