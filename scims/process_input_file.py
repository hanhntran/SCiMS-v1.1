"""
process_input_file.py: Process idxstats files and classify chromosomal sex.

Uses the analytical multinomial likelihood classifier.
"""

import os
import logging

import numpy as np
import pandas as pd

from .analytical import classify_sample as classify_analytical
from .utils import extract_sample_id

logger = logging.getLogger("scims")


def process_input_file(
    input_file: str,
    scaffold_ids: list,
    args,
) -> dict:
    """
    Process a single input file (.idxstats) and return classification results.

    Parameters
    ----------
    input_file : str
        Path to the .idxstats file.
    scaffold_ids : list
        List of scaffold names to use.
    args : argparse.Namespace
        Must contain: homogametic_scaffold, heterogametic_scaffold, ZW, threshold.
        May contain: mismapping_rate.

    Returns
    -------
    dict with classification results, or a dict with 'Status' key on failure.
    """
    sample_id = extract_sample_id(os.path.basename(input_file))

    try:
        idxstats = pd.read_table(input_file, header=None, index_col=0)

        # Subset to scaffolds of interest (keep only rows that exist)
        available = [s for s in scaffold_ids if s in idxstats.index]
        if not available:
            raise ValueError(f"None of the specified scaffolds found in {input_file}")
        idxstats = idxstats.loc[available]

        is_zw = getattr(args, "ZW", False)
        threshold = getattr(args, "threshold", 0.95)
        mismapping_rate = getattr(args, "mismapping_rate", 1e-4)

        # Analytical multinomial likelihood classification
        classification_info = classify_analytical(
            idxstats=idxstats,
            scaffolds=scaffold_ids,
            homogametic_id=args.homogametic_scaffold,
            heterogametic_id=args.heterogametic_scaffold,
            is_zw=is_zw,
            threshold=threshold,
            mismapping_rate=mismapping_rate,
        )
        return _format_result(sample_id, classification_info, is_zw)

    except Exception as exc:
        logger.error(f"Error processing {input_file}: {exc}")
        return {
            "SCiMS_ID": sample_id or "Unknown",
            "Status": f"Failed: {exc}",
        }


def _format_result(sample_id: str, info: dict, is_zw: bool) -> dict:
    """Format analytical classifier output into standard SCiMS result dict."""
    if is_zw:
        return {
            "SCiMS_ID": sample_id,
            "SCiMS_sex": info["SCiMS predicted sex"],
            "SCiMS_reads_mapped": info["Total reads mapped"],
            "SCiMS_reads_mapped_to_Z": info["Reads mapped to Z"],
            "SCiMS_reads_mapped_to_W": info["Reads mapped to W"],
            "SCiMS_male_post_prob": np.round(info["Posterior probability of being male"], 6),
            "SCiMS_female_post_prob": np.round(info["Posterior probability of being female"], 6),
            "SCiMS_log_likelihood_ratio": info.get("Log-likelihood ratio", np.nan),
        }
    else:
        return {
            "SCiMS_ID": sample_id,
            "SCiMS_sex": info["SCiMS predicted sex"],
            "SCiMS_reads_mapped": info["Total reads mapped"],
            "SCiMS_reads_mapped_to_X": info["Reads mapped to X"],
            "SCiMS_reads_mapped_to_Y": info["Reads mapped to Y"],
            "SCiMS_male_post_prob": np.round(info["Posterior probability of being male"], 6),
            "SCiMS_female_post_prob": np.round(info["Posterior probability of being female"], 6),
            "SCiMS_log_likelihood_ratio": info.get("Log-likelihood ratio", np.nan),
        }