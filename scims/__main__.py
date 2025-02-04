"""
SCiMS: Sex Calling in Metagenomic Sequencing

This script classifies host sex using metagenomic sequencing data alone. Metagenomic samples are classified into male, female, or uncertain 
based on coverage ratios of putative sex chromosomes (X/Y or Z/W). It uses a kernel density estimation (KDE) approach, comparing coverage ratios 
against training data. In the XY system, 'male' is heterogametic (XY); in the ZW system, 'female' is heterogametic (ZW).

Author: Hanh Tran

Version: 1.1.0
"""

import argparse
import logging
import os
import re
import numpy as np
import pandas as pd
import sys

from scipy.stats import gaussian_kde

from .utils import (
    read_metadata,
    find_sample_id_column,
    read_master_file,
    extract_sample_id
)

from .helpers import (
    load_training_data
)

from .classification import (
    process_sample_xy,
    process_sample_zw
)


# -------------------------------------------------------------------------------------------------
# CONFIGURE LOGGING
# -------------------------------------------------------------------------------------------------
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


# -------------------------------------------------------------------------------------------------
# MAIN FUNCTION
# -------------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="SCiMS: Sex Calling in Metagenomic Sequencing")
    parser.add_argument('--scaffolds', dest="scaffold_ids_file", required=True, help='Path to the text file containing scaffold IDs')
    parser.add_argument('--metadata', required=True, help='Path to the metadata file')
    parser.add_argument('--master_file', required=True, help='Path to the master file')
    parser.add_argument('--homogametic_id', dest="x_id", required=True, help='ID of the homogametic chromosome in reference genome (e.g. X or Z)')
    parser.add_argument('--heterogametic_id', dest="y_id", required=True, help='ID of the heterogametic chromosome in reference genome (e.g. Y or W)')
    parser.add_argument('--output', dest="output_file", required=True, help='Path to the output file')
    parser.add_argument('--system', dest="system", choices=['XY', 'ZW'], required=True, help='Sex determination system')
    parser.add_argument('--threshold', dest="threshold", type=float, default=0.95, help='Probability threshold for determining sex. Default is 0.95.')
    parser.add_argument('--id_column', dest="id_column", help='User-specified sample ID column name in metadata. If not provided, SCiMS will attempt to locate possible column names from the metadata file.')
    parser.add_argument('--training_data', dest="training_data", help='Path to the training data file if available. SCiMS will use the default training data if not provided.', default="training_data_hmp_1000x_normalizedXY.txt")
    args = parser.parse_args()

    print("=================================================")
    print("""
    _|_|_|   _|_|_|  _|_|_|  _|      _|   _|_|_|  
    _|      _|         _|    _|_|  _|_|   _|        
    _|_|_|  _|         _|    _|  _|  _|   _|_|_|    
        _|  _|         _|    _|      _|       _|  
    _|_|_|   _|_|_|  _|_|_|  _|      _|   _|_|_|    """)
    print("=================================================")
    print("SCiMS: Sex Calling in Metagenomic Sequencing")

    # 1. Load metadata
    try:
        metadata = read_metadata(args.metadata)
    except Exception as e:
        logging.error(f"Failed to read metadata: {e}")
        sys.exit(1)

    # Use the improved function
    try:
        sample_id_col = find_sample_id_column(metadata, args.id_column)
    except ValueError as e:
        logging.error(str(e))
        sys.exit(1)

    logging.info(f"Using '{sample_id_col}' as the sample ID column.")
    
    # 2. Load master file (containing list of idxstats paths)
    idxstats_files = read_master_file(args.master_file)

    # 3. Load training data
    training_data = load_training_data(args.training_data)

    # Extract XY/ZW coverage ratio columns for 'male' and 'female'
    if args.system == 'XY':
        # Use the XY labeling
        male_rows = training_data[training_data['actual_sex'] == 'male']
        female_rows = training_data[training_data['actual_sex'] == 'female']
        # Then extract coverage arrays for KDE
        male_data = male_rows[['Rx', 'Ry']].dropna().values.T
        female_data = female_rows[['Rx', 'Ry']].dropna().values.T
    elif args.system == 'ZW':
        # Use the ZW labeling
        male_rows = training_data[training_data['actual_sex_zw'] == 'male']
        female_rows = training_data[training_data['actual_sex_zw'] == 'female']
        # Then extract coverage arrays for KDE
        male_data = male_rows[['Rz', 'Rw']].dropna().values.T
        female_data = female_rows[['Rz', 'Rw']].dropna().values.T


    # Build KDE models for male and female
    kde_male_joint = gaussian_kde(male_data)
    kde_female_joint = gaussian_kde(female_data)

    # 4. Load scaffold IDs
    with open(args.scaffold_ids_file, 'r') as sf:
        scaffold_ids = [line.strip() for line in sf if line.strip()]

    # 5. Classify each idxstats file
    results = []
    for idxstats_file in idxstats_files:
        sample_id = None
        try:
            sample_id = extract_sample_id(os.path.basename(idxstats_file))
            idxstats = pd.read_table(idxstats_file, header=None, index_col=0)

            # Subser to scaffolds of interest only
            idxstats = idxstats.loc[scaffold_ids]

            # Classification logic
            if args.system == 'XY':
                classification_info = process_sample_xy(
                    idxstats,
                    x_id=args.x_id,
                    y_id=args.y_id,
                    male_kde=kde_male_joint,
                    female_kde=kde_female_joint,
                    threshold=args.threshold
                )
                results.append({
                    'SCiMS sample ID': sample_id,
                    'SCiMS predicted sex': classification_info['SCiMS predicted sex'],
                    'Total reads mapped': classification_info['Total reads mapped'],
                    'Reads mapped to X': classification_info['Reads mapped to X'],
                    'Reads mapped to Y': classification_info['Reads mapped to Y'],
                    'Posterior probability of being male': np.round(classification_info['Posterior probability of being male'], 3),
                    'Posterior probability of being female': np.round(classification_info['Posterior probability of being female'], 3),
                    'Status': 'Success'
                })

            # ZW system
            elif args.system == 'ZW':
                classification_info = process_sample_zw(
                    idxstats,
                    z_id=args.x_id,
                    w_id=args.y_id,
                    male_kde=kde_male_joint,
                    female_kde=kde_female_joint,
                    threshold=args.threshold
                )
                results.append({
                    'SCiMS sample ID': sample_id,
                    'SCiMS predicted sex': classification_info['SCiMS predicted sex'],
                    'Total reads mapped': classification_info['Total reads mapped'],
                    'Reads mapped to Z': classification_info['Reads mapped to Z'],
                    'Reads mapped to W': classification_info['Reads mapped to W'],
                    'Posterior probability of being male': np.round(classification_info['Posterior probability of being male'], 3),
                    'Posterior probability of being female': np.round(classification_info['Posterior probability of being female'], 3),
                    'Status': 'Success'
                    #'Rz': classification_info['Rz'],
                    #'Rw': classification_info['Rw']
                })

        except Exception as exc:
            logging.error(f"Error processing {idxstats_file}: {exc}")
            results.append({
                'SCiMS sample ID': sample_id or 'Unknown',
                'Status': f'Failed: {exc}'
            })
    
    # Merge results with metadata
    results_df = pd.DataFrame(results)

    merged_df = pd.merge(metadata, results_df, left_on=sample_id_col, right_on='SCiMS sample ID', how='left')
    
    # 6. Write results to output file
    merged_df.to_csv(args.output_file, sep='\t', index=False)
    logging.info(f"Results written to {args.output_file}")


if __name__ == "__main__":
    main()