from platform import python_revision
import pytest
import pandas as pd
import numpy as np
from unittest.mock import MagicMock, patch
from types import SimpleNamespace

# import functions to test
from scims.process_idxstats import process_idxstats_file

@pytest.fixture
def mock_args():
    """
    create mock args namespace
    """
    return SimpleNamespace(
        ZW=False,
        homogametic_scaffold='chrX',
        heterogametic_scaffold='chrY',
        threshold=0.95
    )

@pytest.fixture
def mock_kdes():
    """Returns two dummy objects to represent the KDE models."""
    return MagicMock(), MagicMock()


@patch('scims.process_idxstats.process_sample_xy')
@patch('scims.process_idxstats.extract_sample_id')
@patch('scims.process_idxstats.pd.read_table')

## test for success case
def test_process_idxstats_file_success(mock_read_table, mock_extract, mock_process_xy, mock_args, mock_kdes):
    # 1. Setup Mocks
    mock_extract.return_value = "Sample_001"    
    # Create a dummy dataframe to simulate reading an idxstats file
    # Index must contain the scaffolds we request
    dummy_df = pd.DataFrame(
        {'length': [1000, 500], 'reads': [100, 50]}, 
        index=['chrX', 'chrY']
    )
    mock_read_table.return_value = dummy_df
    
    # Mock the return value of the classification logic
    mock_process_xy.return_value = {
        'SCiMS predicted sex': 'Female',
        'Total reads mapped': 1000,
        'Reads mapped to X': 500,
        'Reads mapped to Y': 10,
        'Posterior probability of being male': 0.01,
        'Posterior probability of being female': 0.99
    }

    # 2. Execution
    result = process_idxstats_file(
        idxstats_file="fake/path/Sample_001.idxstats",
        scaffold_ids=["chrX", 'chrY'],
        args=mock_args,
        kde_male_joint=mock_kdes[0],
        kde_female_joint=mock_kdes[1]
    )

    # 3. Assertions
    # Did it extract the ID correctly?
    assert result['SCiMS_ID'] == "Sample_001"
    # Did it map the classification output to the final result dictionary keys?
    assert result['SCiMS_sex'] == "Female"
    assert result['SCiMS_reads_mapped_to_X'] == 500
    assert result['SCiMS_male_post_prob'] == 0.01
    
    # Ensure correct function was called
    mock_process_xy.assert_called_once()



@patch('scims.process_idxstats.extract_sample_id')
@patch('scims.process_idxstats.pd.read_table')
def test_process_file_error_handling(mock_read_table, mock_extract, mock_args, mock_kdes):
    """
    Test that the function catches exceptions (like missing files or bad formats)
    and returns a status dictionary instead of crashing.
    """
    mock_extract.return_value = "Corrupt_Sample"
    
    # Simulate a file reading error
    mock_read_table.side_effect = Exception("File format not supported")

    result = process_idxstats_file(
        idxstats_file="bad_file.txt",
        scaffold_ids=["chrX"],
        args=mock_args,
        kde_male_joint=mock_kdes[0],
        kde_female_joint=mock_kdes[1]
    )

    # Should return a dict with 'Status', not crash
    assert result['SCiMS_ID'] == "Corrupt_Sample"
    assert "Failed" in result['Status']
    assert "File format not supported" in result['Status']