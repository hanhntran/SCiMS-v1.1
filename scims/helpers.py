# helpers.py: functions to process coverage ratios (Rx and Ry)and compute joint posteriors

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import os
import pkg_resources

def compute_joint_posterior(
        Rx: float,
        Ry: float,
        male_kde: gaussian_kde,
        female_kde: gaussian_kde,
        P_male: float,
        P_female: float,
        eps: float = 1e-10
) -> tuple:
    """
    Given coverage ratios for a sample (Rx, Ry) and KDE models for male/female
    plus prior probabilities for each, compute the joint posterior probabilities.
    """
    point = np.array([Rx, Ry])
    male_likelihood = male_kde(point)
    female_likelihood = female_kde(point)

    P_male_joint = male_likelihood * P_male
    P_female_joint = female_likelihood * P_female

    total_joint = P_male_joint + P_female_joint + eps

    P_male_posterior = P_male_joint / total_joint
    #P_female_posterior = P_female_joint / total_joint
    P_female_posterior = 1 - P_male_posterior

    return P_male_posterior[0], P_female_posterior[0]


def determine_sex_with_joint_posteriors(
        P_male_posterior: float,
        P_female_posterior: float,
        threshold: float
) -> str:
    """
    Given the joint posterior probabilities for a sample, determine the sex.
    """
    if P_male_posterior >= threshold:
        return 'male'
    elif P_female_posterior >= threshold:
        return 'female'
    else:
        return 'uncertain'


def compute_coverage_ratio_rx(idxstats: pd.DataFrame, autosome_indices: list, x_index: int) -> float:
    """
    Computes the coverage ratio for the sepecified chromosome index 
    relative to the mean coverage of provided autosome indices.
    """
    coverage_X = idxstats.iloc[x_index, 1] / idxstats.iloc[x_index, 0]
    #coverage_Y = idxstats.iloc[y_index, 1] / idxstats.iloc[y_index, 0]
    #coverage_chrom = idxstats.iloc[autosome_indices, 1] / idxstats.iloc[autosome_indices, 0] # .iloc[:, 0] = chromosome length, .iloc[:, 1] = number of mapped reads
    
    autosome_covs = [
        idxstats.iloc[i, 1] / idxstats.iloc[i, 0]
        for i in autosome_indices
    ]
    Rx = coverage_X / np.mean(autosome_covs) if autosome_covs else np.nan
    return Rx


def load_training_data(training_path: str = None) -> pd.DataFrame:
    """
    Reads training data from a tab-delimited file.
    Expects columns like 'actual_sex_xy' and 'actual_sex_zw', among others.

    If no training_path is provided or if the provided path matches the default
    filename and is not found on disk, the default training data bundled with the
    package is used.
    """
    default_filename = "training_data_hmp_1000x_normalizedXY.txt"
    
    # Use default if no training path is provided
    if training_path is None or training_path == default_filename:
        training_path = pkg_resources.resource_filename('scims', f"data/{default_filename}")
    else:
        # If a training_path was provided but the file doesn't exist locally,
        # try to look for it in the package data.
        if not os.path.exists(training_path):
            pkg_data_path = pkg_resources.resource_filename('scims', f"data/{training_path}")
            if os.path.exists(pkg_data_path):
                training_path = pkg_data_path
            else:
                raise FileNotFoundError(f"Training data file not found: {training_path}")
    
    training_data = pd.read_csv(training_path, sep='\t')
    return training_data


