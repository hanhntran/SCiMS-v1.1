import os
import logging
import pandas as pd
import numpy as np

from .classification import (
    process_sample_xy,
    process_sample_zw
)

from .utils import (
    extract_sample_id
)

def process_idxstats_file(idxstats_file, scaffold_ids, args, kde_male_joint, kde_female_joint):
    """
    Process a single idxstats file and return a dictionary with classification results.
    """
    sample_id = extract_sample_id(os.path.basename(idxstats_file))
    try:
        idxstats = pd.read_table(idxstats_file, header=None, index_col=0)
        # Subset to scaffolds of interest
        idxstats = idxstats.loc[scaffold_ids]

        if not args.ZW:  # Default to XY system
            classification_info = process_sample_xy(
                idxstats,
                x_id=args.homogametic_scaffold,
                y_id=args.heterogametic_scaffold,
                male_kde=kde_male_joint,
                female_kde=kde_female_joint,
                threshold=args.threshold
            )
            result = {
                'SCiMS_ID': sample_id,
                'SCiMS_sex': classification_info['SCiMS predicted sex'],
                'SCiMS_reads_mapped': classification_info['Total reads mapped'],
                'SCiMS_reads_mapped_to_X': classification_info['Reads mapped to X'],
                'SCiMS_reads_mapped_to_Y': classification_info['Reads mapped to Y'],
                'SCiMS_male_post_prob': np.round(classification_info['Posterior probability of being male'], 3),
                'SCiMS_female_post_prob': np.round(classification_info['Posterior probability of being female'], 3)
            }
        else:  # ZW system
            classification_info = process_sample_zw(
                idxstats,
                z_id=args.homogametic_scaffold,
                w_id=args.heterogametic_scaffold,
                male_kde=kde_male_joint,
                female_kde=kde_female_joint,
                threshold=args.threshold
            )
            result = {
                'SCiMS_ID': sample_id,
                'SCiMS_sex': classification_info['SCiMS predicted sex'],
                'SCiMS_reads_mapped': classification_info['Total reads mapped'],
                'SCiMS_reads_mapped_to_Z': classification_info['Reads mapped to Z'],
                'SCiMS_reads_mapped_to_W': classification_info['Reads mapped to W'],
                'SCiMS_male_post_prob': np.round(classification_info['Posterior probability of being male'], 3),
                'SCiMS_female_post_prob': np.round(classification_info['Posterior probability of being female'], 3)
            }
    except Exception as exc:
        logging.error(f"Error processing {idxstats_file}: {exc}")
        result = {
            'SCiMS_ID': sample_id or 'Unknown',
            'Status': f'Failed: {exc}'
        }
    return result
