"""
Prepare the annotation for downstream rules.

* If the user supplies a GFF/GFF3 annotation, convert it to GTF via gffread.
* Validate that the configured feature and attribute exist in the final GTF.
"""

import os
from pathlib import Path

from sample_utils import validate_feature_attribute

ANNOTATION_NEEDS_CONVERSION = config["processing_annotation_source"] != config["processing_gtf_file"]


def _annotation_source(wildcards=None):
    return config["processing_annotation_source"]


def _annotation_gtf(wildcards=None):
    return config["processing_gtf_file"]


def _validation_flag(wildcards=None):
    return os.path.join(config["processing_dir"], "reference", ".annotation_validated")


def _validate(gtf_path, feature_type=None, attribute_type=None):
    validate_feature_attribute(
        gtf_path,
        feature_type or config.get("feature_type"),
        attribute_type or config.get("attribute_type"),
    )


if ANNOTATION_NEEDS_CONVERSION:
    rule prepare_annotation:
        input:
            source=_annotation_source,
        output:
            gtf=_annotation_gtf,
            flag=_validation_flag,
        params:
            feature_type=lambda wildcards: config.get("feature_type"),
            attribute_type=lambda wildcards: config.get("attribute_type"),
        threads:
            1
        conda:
            "../../envs/gffread.yaml"
        singularity:
            config.get("singularity_image", "")
        run:
            shell("mkdir -p $(dirname {output.gtf})")
            shell("gffread {input.source} -T -o {output.gtf}")
            _validate(output.gtf, params.feature_type, params.attribute_type)
            Path(output.flag).touch()
else:
    rule prepare_annotation:
        input:
            gtf=_annotation_gtf,
        output:
            gtf=_annotation_gtf,
            flag=_validation_flag,
        params:
            feature_type=lambda wildcards: config.get("feature_type"),
            attribute_type=lambda wildcards: config.get("attribute_type"),
        threads:
            1
        run:
            _validate(input.gtf, params.feature_type, params.attribute_type)
            Path(output.flag).touch()
