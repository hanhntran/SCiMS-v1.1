import pytest
import pandas as pd
import os
import logging
from scims.utils import (
    read_metadata,
    normalize_colname,
    find_sample_id_column,
    read_master_file,
    extract_sample_id
)

def test_normalize_colname():
    """Test that column names are correctly stripped of special chars and lowercased."""
    assert normalize_colname("Sample_ID") == "sampleid"
    assert normalize_colname("Run #1") == "run1"
    assert normalize_colname("ID") == "id"
    assert normalize_colname("already_clean") == "alreadyclean"

def test_extract_sample_id():
    """Test extracting the base filename without extension."""
    # Standard case
    assert extract_sample_id("/path/to/sample1.idxstats") == "sample1"
    # Case with just filename
    assert extract_sample_id("sample2.txt") == "sample2"
    # Case with multiple dots
    assert extract_sample_id("sample3.tar.gz") == "sample3.tar"


def test_read_metadata(tmp_path):
    """Test reading a TSV metadata file."""
    # Create a dummy metadata file
    meta_file = tmp_path / "metadata.tsv"
    meta_file.write_text("id\tsex\nsample1\tmale\nsample2\tfemale", encoding="utf-8")

    # Run function
    df = read_metadata(str(meta_file))

    # Assertions
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert list(df.columns) == ["id", "sex"]
    assert len(df) == 2

def test_read_master_file(tmp_path):
    """Test reading a master file line by line."""
    master_file = tmp_path / "master.txt"
    content = "file1.idxstats\nfile2.idxstats\n  \n" # includes blank line to test stripping
    master_file.write_text(content, encoding="utf-8")

    lines = read_master_file(str(master_file))

    assert len(lines) == 2
    assert lines[0] == "file1.idxstats"
    assert "file2.idxstats" in lines


def test_find_sample_id_explicit_success():
    """Test when user manually specifies a valid column."""
    df = pd.DataFrame({'custom_name': ['s1', 's2'], 'other': [1, 2]})
    
    col = find_sample_id_column(df, user_specified_col='custom_name')
    assert col == 'custom_name'

def test_find_sample_id_explicit_missing():
    """Test error when user specifies a column that doesn't exist."""
    df = pd.DataFrame({'other': [1, 2]})
    
    with pytest.raises(ValueError, match="not found in metadata"):
        find_sample_id_column(df, user_specified_col='missing_col')

def test_find_sample_id_auto_detect():
    """Test auto-detection of common names like 'Sample_ID' or 'Run'."""
    # Case 1: 'Sample_ID' -> normalizes to 'sampleid' (match)
    df1 = pd.DataFrame({'Sample_ID': ['s1'], 'Data': [1]})
    assert find_sample_id_column(df1) == 'Sample_ID'

    # Case 2: 'run' -> matches directly
    df2 = pd.DataFrame({'run': ['s1'], 'Data': [1]})
    assert find_sample_id_column(df2) == 'run'

def test_find_sample_id_auto_fail():
    """Test error when no recognizable ID column exists."""
    df = pd.DataFrame({'Age': [20], 'Height': [180]})
    
    with pytest.raises(ValueError, match="No valid sample ID column found"):
        find_sample_id_column(df)

def test_find_sample_id_ambiguous(caplog):
    """
    Test warning when multiple valid ID columns exist (e.g., 'id' and 'run').
    """
    df = pd.DataFrame({'id': ['s1'], 'run': ['s1']})
    
    with caplog.at_level(logging.WARNING):
        col = find_sample_id_column(df)
        
    # It should succeed (pick one arbitrarily) but log a warning
    assert col in ['id', 'run']
    assert "Multiple columns could match" in caplog.text