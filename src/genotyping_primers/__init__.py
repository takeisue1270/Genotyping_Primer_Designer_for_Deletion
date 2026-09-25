"""Genotyping primers for a defined deletion.

Give it a sequence and the block to remove; get back the PCR that separates the
edited allele from the unedited one by product size, with both primers kept
clear of any repair template's homology arms.

    from genotyping_primers import DesignConfig, design, prepare, target_from_file

    target, notes = target_from_file("locus.fa")
    target, deletion = prepare(target, 4799, 5300)   # 0-based half-open
    result = design(target, deletion, DesignConfig(homology_arm=43))
    print(result.best.wt_product, result.best.ko_product)
"""

from .config import DesignConfig
from .design import (
    DesignError,
    design,
    edited_sequence,
    prepare,
    resolve_span,
    target_from_file,
    target_from_genome,
    with_insert,
)
from .models import Deletion, DesignResult, Primer, PrimerPair, SpecificityHits, Target
from .report import text_report, write_maps, write_tsv

__version__ = "0.1.0"

__all__ = [
    "Deletion",
    "DesignConfig",
    "DesignError",
    "DesignResult",
    "Primer",
    "PrimerPair",
    "SpecificityHits",
    "Target",
    "__version__",
    "design",
    "edited_sequence",
    "prepare",
    "resolve_span",
    "target_from_file",
    "target_from_genome",
    "text_report",
    "with_insert",
    "write_maps",
    "write_tsv",
]
