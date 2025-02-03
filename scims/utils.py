# utils.py: functions to process metadata and master file for SCiMS

import os
import logging
import pandas as pd

def read_metadata(metadata_path: str) -> pd.DataFrame:
    """
    Read metadata from a file.
    Assume a tab-delimitedd file.
    """
    return pd.read_csv(metadata_path, sep='\t')

def normalize_colname(colname: str) -> str:
    """Helper to normalize column names by stripping non-alphanumeric and converting to lowercase."""
    return "".join(ch.lower() for ch in colname if ch.isalnum())

def find_sample_id_column(
    metadata: pd.DataFrame,
    user_specified_col: str = None
) -> str:
    """
    Attempts to locate the sample ID column in a robust way:
      1. If user_specified_col is provided, verify it exists.
      2. Otherwise, try known column names (case-insensitive, ignoring underscores).
      3. If still not found, raise a ValueError with a helpful message.
    """
    if user_specified_col:
        # User told us explicitly which column
        if user_specified_col in metadata.columns:
            return user_specified_col
        raise ValueError(
            f"User-specified column '{user_specified_col}' not found in metadata. "
            f"Available columns: {list(metadata.columns)}"
        )

    # Attempt auto-detection
    known_patterns = {"sampleid", "id", "run", "sra", "featureid"}

    # Build a map from normalized -> original column name
    col_map = {normalize_colname(c): c for c in metadata.columns}
    intersection = known_patterns.intersection(col_map.keys())

    if intersection:
        # If there's exactly one match, return it
        if len(intersection) == 1:
            matched_pattern = next(iter(intersection))
            chosen = col_map[matched_pattern]
            logging.info(f"Auto-detected sample ID column: '{chosen}'")
            return chosen
        else:
            # More than one possible match
            # We'll pick the first arbitrarily and warn the user
            matched_pattern = next(iter(intersection))
            chosen = col_map[matched_pattern]
            logging.warning(
                f"Multiple columns could match sample ID patterns: {intersection}. "
                f"Defaulting to '{chosen}'. Please specify --id-column if this is incorrect."
            )
            return chosen

    # No match found
    raise ValueError(
        "No valid sample ID column found. "
        "Please use the '--id-column' parameter or rename one of your columns. "
        f"Metadata columns are: {list(metadata.columns)}"
    )


def read_master_file(master_file_path: str) -> list:
    """
    Reads the master file containing paths to the idxstats output,
    one file path per line. Returns a list of file paths.
    """
    with open(master_file_path, 'r') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]
    return lines


def extract_sample_id(file_path: str) -> str:
    """
    Returns the base filename (minus extension) as the sample ID.
    Example: 'my_sample_001.idxstats' -> 'my_sample_001'
    """
    base_name = os.path.basename(file_path)      # e.g. 'my_sample_001.idxstats'
    root, ext = os.path.splitext(base_name)      # ('my_sample_001', '.idxstats')
    return root
