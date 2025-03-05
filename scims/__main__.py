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

from .process_idxstats import (
    process_idxstats_file
)


# -------------------------------------------------------------------------------------------------
# CONFIGURE LOGGING
# -------------------------------------------------------------------------------------------------
#logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


# -------------------------------------------------------------------------------------------------
# MAIN FUNCTION
# -------------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="SCiMS: Sex Calling in Metagenomic Sequences")

    # Common arguments
    parser.add_argument('--scaffolds', dest="scaffold_ids_file", required=True, help='Path to the text file containing scaffold IDs')
    parser.add_argument('--homogametic_id', dest="x_id", required=True, help='ID of the homogametic chromosome (e.g. X or Z)')
    parser.add_argument('--heterogametic_id', dest="y_id", required=True, help='ID of the heterogametic chromosome (e.g. Y or W)')
    parser.add_argument('--system', dest="system", choices=['XY', 'ZW'], required=True, help='Sex determination system')
    parser.add_argument('--threshold', dest="threshold", type=float, default=0.95, help='Probability threshold for determining sex')
    parser.add_argument('--training_data', dest="training_data", help='Path to the training data file', default="training_data_hmp_1000x_normalizedXY.txt")
    
    # Mode options
    # Default is single-sample mode
    parser.add_argument('--idxstats_file', help='Path to a single idxstats file (default mode)')
    # Option to process multiple files
    parser.add_argument('--multiple', action='store_true', help="Flag to process multiple files from a folder")
    parser.add_argument('--idxstats_folder', help='Path to the folder containing idxstats files (required if --multiple is set)')
    
    # For metadata update (optional, used only in multiple mode)
    parser.add_argument('--metadata', help='Path to the metadata file (optional, used in multiple mode)')
    parser.add_argument('--id_column', dest="id_column", help='User-specified sample ID column name in metadata')
    
    parser.add_argument('--output', dest="output_file", required=True, help='Path to the output file')

    # Add log file argument
    parser.add_argument('--log', dest="log_file", help='Path to the log file')
    args = parser.parse_args()
    #Create logger and set its level
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # File handler if log file is specified
    if args.log_file:
        file_handler = logging.FileHandler(args.log_file, mode='w')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    # Print header
    logger.info("=================================================")
    logger.info(
        """
        _|_|_|   _|_|_|  _|_|_|  _|      _|   _|_|_|  
        _|      _|         _|    _|_|  _|_|   _|        
        _|_|_|  _|         _|    _|  _|  _|   _|_|_|    
            _|  _|         _|    _|      _|       _|  
        _|_|_|   _|_|_|  _|_|_|  _|      _|   _|_|_|    
        
        ================================================="""
    )
    logger.info("SCiMS: Sex Calling in Metagenomic Sequences")
    logger.info("Version: 1.1.0")
    logger.info("=================================================")


    # 1. Load scaffold IDs (common to both modes)
    try:
        with open(args.scaffold_ids_file, 'r') as sf:
            scaffold_ids = [line.strip() for line in sf if line.strip()]
    except Exception as e:
        logger.error(f"Failed to read scaffold IDs: {e}")
        sys.exit(1)

    # 2. Load training data and build KDE models (common to both modes)
    try:
        training_data = load_training_data(args.training_data)
    except Exception as e:
        logger.error(f"Failed to load training data: {e}")
        sys.exit(1)
    
    if args.system == 'XY':
        male_rows = training_data[training_data['actual_sex'] == 'male']
        female_rows = training_data[training_data['actual_sex'] == 'female']
        male_data = male_rows[['Rx', 'Ry']].dropna().values.T
        female_data = female_rows[['Rx', 'Ry']].dropna().values.T
    else:  # ZW system
        male_rows = training_data[training_data['actual_sex_zw'] == 'male']
        female_rows = training_data[training_data['actual_sex_zw'] == 'female']
        male_data = male_rows[['Rz', 'Rw']].dropna().values.T
        female_data = female_rows[['Rz', 'Rw']].dropna().values.T

    kde_male_joint = gaussian_kde(male_data)
    kde_female_joint = gaussian_kde(female_data)

    results = []

    # Check mode and required arguments
    if args.multiple:
        if not args.idxstats_folder:
            logger.error("For multiple-sample mode, --idxstats_folder is required.")
            sys.exit(1)
        # Get all idxstats files in the provided folder (filter by extension)
        folder_files = [os.path.join(args.idxstats_folder, f) for f in os.listdir(args.idxstats_folder) if f.endswith(".idxstats")]
        if not folder_files:
            logger.error("No idxstats files found in the folder.")
            sys.exit(1)
        for idxstats_file in folder_files:
            results.append(process_idxstats_file(idxstats_file, scaffold_ids, args, kde_male_joint, kde_female_joint))
    else:
        # Single sample mode (default)
        if not args.idxstats_file:
            logger.error("For single-sample mode, --idxstats_file is required.")
            sys.exit(1)
        results.append(process_idxstats_file(args.idxstats_file, scaffold_ids, args, kde_male_joint, kde_female_joint))
        logger.info(f"SCiMS predicting sex for {args.idxstats_file}: {results}")
    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    
    # If in multiple mode and metadata is provided, update the metadata file
    if args.multiple and args.metadata:
        try:
            metadata = read_metadata(args.metadata)
            sample_id_col = find_sample_id_column(metadata, args.id_column)
            merged_df = pd.merge(metadata, results_df, left_on=sample_id_col, right_on='SCiMS_ID', how='left')
            merged_df.drop(columns=['SCiMS_ID'], inplace=True)
            merged_df.to_csv(args.output_file, sep='\t', index=False)
            logger.info(f"Updated metadata with classification results written to {args.output_file}")
        except Exception as e:
            logger.error(f"Error updating metadata: {e}")
            sys.exit(1)
    else:
        # For single-sample mode or if no metadata is provided, simply output the results
        results_df.to_csv(args.output_file, sep='\t', index=False)
        logger.info(f"Results written to {args.output_file}")

if __name__ == "__main__":
    main()