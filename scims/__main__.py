"""
SCiMS: Sex Calling in Metagenomic Sequencing

This script classifies host sex using metagenomic sequencing data alone. 
Metagenomic samples are classified into male, female, or uncertain based on 
coverage ratios of putative sex chromosomes (X/Y or Z/W). It uses a kernel 
density estimation (KDE) approach, comparing coverage ratios against training data.
In the XY system, 'male' is heterogametic (XY); in the ZW system, 'female' is 
heterogametic (ZW).

Author: Hanh Tran
Version: {__version__}
"""

import argparse
import logging
import os
import pandas as pd
import sys
from scims import __version__ 
from scipy.stats import gaussian_kde

from .utils import (
    read_metadata,
    find_sample_id_column,
    extract_sample_id
)
from .helpers import load_training_data
from .process_idxstats import process_idxstats_file

def main():
    parser = argparse.ArgumentParser(description="SCiMS: Sex Calling in Metagenomic Sequences")
    
    # Mode selection: default is single-sample mode.
    parser.add_argument('--idxstats_file', dest="idxstats_file", help='Path to a single idxstats file (default mode)')
    parser.add_argument('--idxstats_folder', dest="idxstats_folder", help='Path to the folder containing idxstats files for multiple-sample mode')
    
    parser.add_argument('--scaffolds', dest="scaffold_ids_file", required=True, help='Path to the text file containing scaffold IDs')
    parser.add_argument('--homogametic_scaffold', dest="homogametic_scaffold", required=True, help='ID of the homogametic scaffold (e.g. X or Z)')
    parser.add_argument('--heterogametic_scaffold', dest="heterogametic_scaffold", required=True, help='ID of the heterogametic scaffold (e.g. Y or W)')
    parser.add_argument('--ZW', dest="ZW", action="store_true", help='Switch to ZW system (default is XY)')
    parser.add_argument('--threshold', dest="threshold", type=float, default=0.5, help='Probability threshold for determining sex')
    parser.add_argument('--output_dir', dest="output_dir", required=True, help='Path to the output directory')
    parser.add_argument('--training_data', dest="training_data", help='Path to the training data file', default="training_data_hmp_1000x_normalizedXY.txt")
    
    # Optional metadata update (only used in multiple-sample mode)
    parser.add_argument('--metadata', dest="metadata", help='Path to the metadata file (optional, used in multiple-sample mode)')
    parser.add_argument('--id_column', dest="id_column", help='User-specified sample ID column name in metadata')
    
    # New boolean flag for log output (default is False)
    parser.add_argument('--log', dest="log", action="store_true", help='If set, a log file is written to the output directory (scims.log)')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')    
    args = parser.parse_args()
    
    # Check metadata parameters early
    if args.metadata and not args.id_column:
        print("Error: When providing a metadata file, you must also specify the id column using --id_column.")
        sys.exit(1)
    
    # Setup logging after parsing arguments
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    # Clear any existing handlers
    if logger.hasHandlers():
        logger.handlers.clear()
    
    # Add console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Add file handler if the --log flag is set
    if args.log:
        log_file_path = os.path.join(args.output_dir, "scims.log")
        file_handler = logging.FileHandler(log_file_path, mode='w')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        logger.info(f"Log file created at: {log_file_path}")
    
    logger.info(" \n=================================================")
    logger.info("""
    _|_|_|   _|_|_|  _|_|_|  _|      _|   _|_|_|  
    _|      _|         _|    _|_|  _|_|   _|        
    _|_|_|  _|         _|    _|  _|  _|   _|_|_|    
        _|  _|         _|    _|      _|       _|  
    _|_|_|   _|_|_|  _|_|_|  _|      _|   _|_|_|    
    =================================================""")
    logger.info("SCiMS: Sex Calling in Metagenomic Sequencing")
    logger.info(f"Version: {__version__}")
    logger.info("=================================================")
    
    # Validate mode: either a single file or folder must be provided
    if args.idxstats_folder:
        mode = "multiple"
    elif args.idxstats_file:
        mode = "single"
    else:
        logger.error("You must specify either --idxstats_file for single-sample mode or --idxstats_folder for multiple-sample mode.")
        sys.exit(1)
    
    # Load scaffold IDs
    try:
        with open(args.scaffold_ids_file, 'r') as sf:
            scaffold_ids = [line.strip() for line in sf if line.strip()]
    except Exception as e:
        logger.error(f"Failed to read scaffold IDs: {e}")
        sys.exit(1)
    
    # Load training data and build KDE models
    try:
        training_data = load_training_data(args.training_data)
    except Exception as e:
        logger.error(f"Failed to load training data: {e}")
        sys.exit(1)
    
    if args.ZW:
        # Use ZW system columns from training data
        male_rows = training_data[training_data['actual_sex_zw'] == 'male']
        female_rows = training_data[training_data['actual_sex_zw'] == 'female']
        male_data = male_rows[['Rz', 'Rw']].dropna().values.T
        female_data = female_rows[['Rz', 'Rw']].dropna().values.T
    else:
        # Default to XY system
        male_rows = training_data[training_data['actual_sex'] == 'male']
        female_rows = training_data[training_data['actual_sex'] == 'female']
        male_data = male_rows[['Rx', 'Ry']].dropna().values.T
        female_data = female_rows[['Rx', 'Ry']].dropna().values.T
    
    kde_male_joint = gaussian_kde(male_data)
    kde_female_joint = gaussian_kde(female_data)
    
    if mode == "multiple":
        try:
            folder_files = [os.path.join(args.idxstats_folder, f) 
                            for f in os.listdir(args.idxstats_folder) if f.endswith(".idxstats")]
        except Exception as e:
            logger.error(f"Error reading idxstats folder: {e}")
            sys.exit(1)
            
        if not folder_files:
            logger.error("No idxstats files found in the provided folder.")
            sys.exit(1)
        
        all_results = []  # To collect results for optional metadata merging
        for idxstats_file in folder_files:
            result = process_idxstats_file(idxstats_file, scaffold_ids, args, kde_male_joint, kde_female_joint)
            all_results.append(result)
            # Build output dictionary using consistent keys
            out_dict = {
                "SCiMS_ID": result.get("SCiMS_ID"),
                "SCiMS_predicted_sex": result.get("SCiMS_sex"),
                "SCiMS_male_post_prob": result.get("SCiMS_male_post_prob"),
                "SCiMS_female_post_prob": result.get("SCiMS_female_post_prob")
            }
            base_name = os.path.basename(idxstats_file).split('.')[0]
            output_file = os.path.join(args.output_dir, f"{base_name}_results.txt")
            pd.DataFrame([out_dict]).to_csv(output_file, sep='\t', index=False)
            logger.info(f"Results written to {output_file}")
        
        # Merge metadata once after processing all files
        if args.metadata:
            try:
                results_df = pd.DataFrame(all_results)
                metadata = read_metadata(args.metadata)
                sample_id_col = find_sample_id_column(metadata, args.id_column)
                merged_df = pd.merge(metadata, results_df, left_on=sample_id_col, right_on='SCiMS_ID', how='left')
                merged_df.drop(columns=['SCiMS_ID'], inplace=True)
                metadata_basename = os.path.basename(args.metadata).split('.')[0]
                metadata_file = os.path.join(args.output_dir, f"{metadata_basename}_scims_updated.txt")
                merged_df.to_csv(metadata_file, sep='\t', index=False)
                logger.info(f"Updated metadata with classification results written to {metadata_file}")
            except Exception as e:
                logger.error(f"Error updating metadata: {e}")
                sys.exit(1)
    else:
        # Single-sample mode: process one file and write the filtered result
        result = process_idxstats_file(args.idxstats_file, scaffold_ids, args, kde_male_joint, kde_female_joint)
        out_dict = {
            "SCiMS_ID": result.get("SCiMS_ID"),
            "SCiMS_predicted_sex": result.get("SCiMS_sex"),
            "SCiMS_male_post_prob": result.get("SCiMS_male_post_prob"),
            "SCiMS_female_post_prob": result.get("SCiMS_female_post_prob")
        }
        results_df = pd.DataFrame([out_dict])
        base_name = os.path.basename(args.idxstats_file).split('.')[0]
        output_file = os.path.join(args.output_dir, f"{base_name}_results.txt")
        results_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Results written to {output_file}")
    
if __name__ == "__main__":
    main()
