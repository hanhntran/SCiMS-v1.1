"""
analytical.py: Training-free sex classification using multinomial likelihood.

Classifies host chromosomal sex by computing the likelihood of observed
per-chromosome read counts under male and female generative models.
No training data required  the expected read distributions are derived
directly from chromosome sizes and ploidy.

Mathematical framework:
    For each chromosome c with length L_c and copy number n_c under sex S:
        p_c^S = (n_c * L_c) / sum_j(n_j * L_j)

    Log-likelihood ratio:
        Lambda = sum_c N_c * [log(p_c^male) - log(p_c^female)]

    Posterior (with uniform prior):
        P(male | data) = sigmoid(Lambda) = 1 / (1 + exp(-Lambda))
"""

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger("scims")


def compute_chromosome_probs(
    chrom_lengths: dict,
    homogametic_id: str,
    heterogametic_id: str,
    sex: str,
    is_zw: bool = False,
    mismapping_rate: float = 1e-4,
) -> dict:
    """
    Compute the expected read probability for each chromosome under a sex model.

    In a diploid organism, the probability of a read landing on chromosome c
    is proportional to (copy_number  length). Sex chromosomes have different
    copy numbers depending on sex:

        XY system: Female = XX (homo2, hetero0), Male = XY (homo1, hetero1)
        ZW system: Male = ZZ (homo2, hetero0), Female = ZW (homo1, hetero1)

    Parameters
    ----------
    chrom_lengths : dict
        {chrom_name: length_in_bp} for all scaffolds.
    homogametic_id : str
        Scaffold ID for the homogametic sex chromosome (X or Z).
    heterogametic_id : str
        Scaffold ID for the heterogametic sex chromosome (Y or W).
    sex : str
        "male" or "female".
    is_zw : bool
        If True, use ZW system (male=ZZ, female=ZW). Default: XY.
    mismapping_rate : float
        Probability floor for zero-ploidy chromosomes, accounting for
        microbial mismapping and alignment errors.

    Returns
    -------
    dict of {chrom_name: probability}
    """
    # Determine ploidy for each chromosome under this sex model
    ploidy = {}
    for chrom in chrom_lengths:
        if chrom == homogametic_id:
            if is_zw:
                # ZW: male is ZZ (homogametic), female is ZW
                ploidy[chrom] = 2 if sex == "male" else 1
            else:
                # XY: female is XX (homogametic), male is XY
                ploidy[chrom] = 2 if sex == "female" else 1
        elif chrom == heterogametic_id:
            if is_zw:
                ploidy[chrom] = 0 if sex == "male" else 1
            else:
                ploidy[chrom] = 0 if sex == "female" else 1
        else:
            # Autosome: diploid in both sexes
            ploidy[chrom] = 2

    # Compute effective size (ploidy  length) for each chromosome
    effective = {c: ploidy[c] * chrom_lengths[c] for c in chrom_lengths}
    total = sum(effective.values())

    # Convert to probabilities, applying mismapping floor
    probs = {}
    for c in chrom_lengths:
        if effective[c] == 0:
            probs[c] = mismapping_rate
        else:
            probs[c] = effective[c] / total

    # Renormalize so probabilities sum to 1
    total_prob = sum(probs.values())
    probs = {c: probs[c] / total_prob for c in probs}

    return probs


def classify_sample(
    idxstats: pd.DataFrame,
    scaffolds: list,
    homogametic_id: str,
    heterogametic_id: str,
    is_zw: bool = False,
    threshold: float = 0.95,
    mismapping_rate: float = 1e-4,
) -> dict:
    """
    Classify chromosomal sex using the analytical multinomial model.

    Parameters
    ----------
    idxstats : pd.DataFrame
        Indexed by chromosome name. Column 0 = length, column 1 = mapped reads.
    scaffolds : list
        Scaffold names to use (autosomes + sex chromosomes).
    homogametic_id : str
        Scaffold ID for the homogametic sex chromosome (X or Z).
    heterogametic_id : str
        Scaffold ID for the heterogametic sex chromosome (Y or W).
    is_zw : bool
        If True, use ZW sex determination system.
    threshold : float
        Posterior probability threshold for confident calls (default: 0.95).
    mismapping_rate : float
        Probability floor for zero-ploidy chromosomes (default: 1e-4).

    Returns
    -------
    dict with keys:
        SCiMS predicted sex : str ("male", "female", or "uncertain")
        Posterior probability of being male : float
        Posterior probability of being female : float
        Total reads mapped : int
        Reads mapped to homogametic : int
        Reads mapped to heterogametic : int
        Log-likelihood ratio : float
        Rx : float (for compatibility with existing output)
        Ry : float (for compatibility with existing output)
    """
    # Subset to scaffolds of interest
    df = idxstats.loc[idxstats.index.isin(scaffolds)].copy()

    # Extract chromosome lengths and read counts
    chrom_lengths = {}
    read_counts = {}
    for chrom in scaffolds:
        if chrom in df.index:
            chrom_lengths[chrom] = int(df.loc[chrom].iloc[0])
            read_counts[chrom] = int(df.loc[chrom].iloc[1])
        else:
            chrom_lengths[chrom] = 0
            read_counts[chrom] = 0

    total_reads = sum(read_counts.values())
    homo_reads = read_counts.get(homogametic_id, 0)
    hetero_reads = read_counts.get(heterogametic_id, 0)

    # Handle edge case: no reads at all
    if total_reads == 0:
        return _empty_result(homogametic_id, heterogametic_id, is_zw)

    # Compute expected probabilities under each sex model
    p_male = compute_chromosome_probs(
        chrom_lengths, homogametic_id, heterogametic_id,
        sex="male", is_zw=is_zw, mismapping_rate=mismapping_rate,
    )
    p_female = compute_chromosome_probs(
        chrom_lengths, homogametic_id, heterogametic_id,
        sex="female", is_zw=is_zw, mismapping_rate=mismapping_rate,
    )

    # Compute log-likelihood ratio: Lambda = sum N_c * log(p_male_c / p_female_c)
    log_ratio = 0.0
    for chrom in scaffolds:
        n = read_counts.get(chrom, 0)
        if n > 0 and p_male[chrom] > 0 and p_female[chrom] > 0:
            log_ratio += n * (np.log(p_male[chrom]) - np.log(p_female[chrom]))

    # Posterior via logistic function (numerically stable)
    if log_ratio > 500:
        P_male_post = 1.0
    elif log_ratio < -500:
        P_male_post = 0.0
    else:
        P_male_post = 1.0 / (1.0 + np.exp(-log_ratio))

    P_female_post = 1.0 - P_male_post

    # Classification
    if P_male_post >= threshold:
        predicted_sex = "male"
    elif P_female_post >= threshold:
        predicted_sex = "female"
    else:
        predicted_sex = "uncertain"


    # Build result dict matching existing SCiMS output format
    if is_zw:
        return {
            "Total reads mapped": total_reads,
            "Reads mapped to Z": homo_reads,
            "Reads mapped to W": hetero_reads,
            "Posterior probability of being male": np.round(P_male_post, 6),
            "Posterior probability of being female": np.round(P_female_post, 6),
            "Log-likelihood ratio": np.round(log_ratio, 3),
            "SCiMS predicted sex": predicted_sex,
        }
    else:
        return {
            "Total reads mapped": total_reads,
            "Reads mapped to X": homo_reads,
            "Reads mapped to Y": hetero_reads,
            "Posterior probability of being male": np.round(P_male_post, 6),
            "Posterior probability of being female": np.round(P_female_post, 6),
            "Log-likelihood ratio": np.round(log_ratio, 3),
            "SCiMS predicted sex": predicted_sex,
        }


def _empty_result(homogametic_id: str, heterogametic_id: str, is_zw: bool) -> dict:
    """Return an uncertain result when no reads are available."""
    if is_zw:
        return {
            "Total reads mapped": 0,
            "Reads mapped to Z": 0, "Reads mapped to W": 0,
            "Posterior probability of being male": 0.5,
            "Posterior probability of being female": 0.5,
            "Log-likelihood ratio": 0.0,
            "SCiMS predicted sex": "uncertain",
        }
    else:
        return {
            "Total reads mapped": 0,
            "Reads mapped to X": 0, "Reads mapped to Y": 0,
            "Posterior probability of being male": 0.5,
            "Posterior probability of being female": 0.5,
            "Log-likelihood ratio": 0.0,
            "SCiMS predicted sex": "uncertain",
        }

