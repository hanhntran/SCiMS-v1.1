# classification.py : core logic for classification process

from scipy.stats import gaussian_kde
import logging
import numpy as np
import pandas as pd
from .helpers import (
    compute_coverage_ratio_rx,
    compute_joint_posterior,
    determine_sex_with_joint_posteriors,
    load_training_data
) 

def process_sample_xy(
        idxstats: pd.DataFrame,
        x_id: str,
        y_id: str,
        male_kde: gaussian_kde,
        female_kde: gaussian_kde,
        threshold: float
) -> dict:
    """
    Given an idxstats DataFrame, XY chromosome IDs (x_id, y_id), and KDE models for male/female,
    computes coverage ratios (Rx, Ry) and determines sex using joint posteriors.
    """
    # Indices of X/Y
    x_index = idxstats.index.get_loc(x_id) if x_id in idxstats.index else None
    y_index = idxstats.index.get_loc(y_id) if y_id in idxstats.index else None
    
    # Idenfity which rows are autosomes (excluding X and Y)
    autosome_indices = [
        i for i in range(len(idxstats))
        if i != x_index and (y_index is None or i != y_index)]
    
    # Coverage ratio for X vs autosomes (Rx)
    Rx = compute_coverage_ratio_rx(idxstats, autosome_indices, x_index)

    # Coverage ratio for Y vs X+Y (Ry)
    x_count = idxstats.loc[x_id].iloc[1] if x_id in idxstats.index else 0    # column 1 = mapped reads column 0 = chromosome length
    y_count = idxstats.loc[y_id].iloc[1] if y_id in idxstats.index else 0
    x_length = idxstats.loc[x_id].iloc[0]
    y_length = idxstats.loc[y_id].iloc[0]
    total_count = idxstats.iloc[:, 1].sum()
    total_xy = x_count + y_count

    if total_xy == 0:
        Ry = np.nan
        logging.warning(f"No reads mapped to X or Y in {x_id} or {y_id}. Skipping Ry ratio calculation.")
    else:
        factor = x_length / y_length # This factor is used to take into account that X and Y are different lengths, so we need to normalize the coverage ratio by the length of the chromosomes so that the coverage ratio is comparable between the two chromosomes and not just because of artifactual differences in length
        Ry = (y_count / total_xy) * factor

    # Compute posteriors
    P_male = 0.5
    P_female = 0.5

    P_male_posterior, P_female_posterior = compute_joint_posterior(Rx, Ry, male_kde, female_kde, P_male, P_female)
    
    inferred_sex = determine_sex_with_joint_posteriors(P_male_posterior, P_female_posterior, threshold)

    return {
        'Rx': Rx,
        'Ry': Ry,
        'Total reads mapped': total_count,
        'Reads mapped to X': x_count,
        'Reads mapped to Y': y_count,
        'Posterior probability of being male': np.round(P_male_posterior, 3),
        'Posterior probability of being female': np.round(P_female_posterior, 3),
        'SCiMS predicted sex': inferred_sex
    }


def process_sample_zw(
        idxstats: pd.DataFrame,
        z_id: str,
        w_id: str,
        male_kde: gaussian_kde,
        female_kde: gaussian_kde,
        threshold: float
) -> dict:
    """
    Given an idxstats DataFrame, ZW chromosome IDs (z_id, w_id), and KDE models for male/female,
    computes coverage ratios (Rz, Rw) and determines sex using joint posteriors.
    """
    z_index = idxstats.index.get_loc(z_id) if z_id in idxstats.index else None
    w_index = idxstats.index.get_loc(w_id) if w_id in idxstats.index else None

    # Identify which rows are autosomes (excluding Z and W)
    autosome_indices = [
        i for i in range(len(idxstats))
        if i != z_index and (w_index is None or i != w_index)
    ]

    # Coverage ratio for Z vs autosomes (Rz)
    Rz = compute_coverage_ratio_rx(idxstats, autosome_indices, z_index)

    # Coverage ratio for W vs Z+W (Rw)
    z_count = idxstats.loc[z_id].iloc[1] if z_id in idxstats.index else 0
    z_length = idxstats.loc[z_id].iloc[0]
    #print(z_length)
    w_count = idxstats.loc[w_id].iloc[1] if w_id in idxstats.index else 0
    w_length = idxstats.loc[w_id].iloc[0]
    #print(w_length)
    total_count = idxstats.iloc[:, 1].sum()
    total_zw = z_count + w_count

    if total_zw == 0:
        Rw = np.nan
        logging.warning(f"No reads mapped to Z or W in {z_id} or {w_id}. Skipping Rw ratio calculation.")
    else:
        #Rw = (w_count / total_zw) * (z_length / w_length)
        #print(Rw)
        factor = z_length / w_length # This factor is used to take into account that Z and W are different lengths, so we need to normalize the coverage ratio by the length of the chromosomes so that the coverage ratio is comparable between the two chromosomes and not just because of artifactual differences in length
        #print(factor)
        Rw = (w_count / total_zw) * factor
        #print(Rw)

    # Compute posteriors
    P_male = 0.5
    P_female = 0.5

    P_male_posterior, P_female_posterior = compute_joint_posterior(Rz, Rw, male_kde, female_kde, P_male, P_female)

    inferred_sex = determine_sex_with_joint_posteriors(P_male_posterior, P_female_posterior, threshold)

    return {
        'Rz': Rz,
        'Rw': Rw,
        'Total reads mapped': total_count,
        'Reads mapped to Z': z_count,
        'Reads mapped to W': w_count,
        'Posterior probability of being male': np.round(P_male_posterior, 3),
        'Posterior probability of being female': np.round(P_female_posterior, 3),
        'SCiMS predicted sex': inferred_sex
    }
