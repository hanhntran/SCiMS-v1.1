import pytest
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
from unittest.mock import MagicMock, patch

# import functions to test
from scims.classification import process_sample_xy, process_sample_zw

@pytest.fixture
def mock_kde():
    """
    create mock Kernel Density Estimation model 
    """
    return MagicMock(spec=gaussian_kde)

@pytest.fixture
def sample_idxstats():
    """
    create mock idxstats DataFrame
    """
    data = {
        'length': [1000, 1000, 500, 250],
        'mapped': [200, 100, 50, 10]
    }
    df = pd.DataFrame(data, index=['chr1', 'chr2', 'chrX', 'chrY'])
    return df
########################################################
# test for XY
########################################################
@patch('scims.classification.determine_sex_with_joint_posteriors')
@patch('scims.classification.compute_joint_posterior')
@patch('scims.classification.compute_coverage_ratio_rx')
def test_process_sample_xy(mock_rx, mock_posterior, mock_sex, sample_idxstats, mock_kde):
    """
    This test checks if the process_sample_xy function correctly processes a sample with XY chromosomes.
    """
    # mock return values
    mock_rx.return_value = 0.5
    mock_posterior.return_value = (0.1, 0.9)
    mock_sex.return_value = 'Female'

    # call the function
    result = process_sample_xy(
        idxstats=sample_idxstats,
        x_id='chrX',
        y_id='chrY',
        male_kde=mock_kde,
        female_kde=mock_kde,
        threshold=0.95
    )

    # Check if helper functions were called correctly
    mock_rx.assert_called_once()

    # Check manual calculations (Logic inside classification.py)
        # X count = 50, Y count = 10, Total XY = 60
        # X len = 500, Y len = 250
        # X_coverage = x_count / x_len = 50 / 500 = 0.1
        # Y_coverage = y_count / y_len = 10 / 250 = 0.04
        # chr1_coverage = chr1_count / chr1_len = 200 / 1000 = 0.2
        # chr2_coverage = chr2_count / chr2_len = 100 / 1000 = 0.1
        # Rx Formula: (X_coverage/chr1_coverage) + (X_coverage/chr2_coverage) = 0.1/0.2 + 0.1/0.1 = 0.5 + 1 = 1.5
        # Ry Formula: y_count / total_xy = (10 / 60) * (500 / 250) = 0.1666... * 2 = 0.333...
    expected_rx =  0.1/0.2 + 0.1/0.1
    expected_ry = (10 / 60) * (500 / 250)

    # Assert the results
    assert result['Rx'] == 0.5
    assert result['Ry'] == pytest.approx(expected_ry)
    assert result['SCiMS predicted sex'] == 'Female'
    assert result['Reads mapped to X'] == 50
    assert result['Reads mapped to Y'] == 10

@patch('scims.classification.determine_sex_with_joint_posteriors')
@patch('scims.classification.compute_joint_posterior')
@patch('scims.classification.compute_coverage_ratio_rx')
def test_process_sample_xy_zero_reads(mock_rx, mock_posterior, mock_sex, sample_idxstats, mock_kde):
    """
    Test edge case: No reads mapped to X or Y. 
    Should handle ZeroDivisionError and set Ry to NaN.
    """
    # Modify data to have 0 reads for X and Y
    sample_idxstats.loc['chrX', 'mapped'] = 0
    sample_idxstats.loc['chrY', 'mapped'] = 0
    
    # Setup Mocks
    mock_rx.return_value = 0.0
    mock_posterior.return_value = (0.5, 0.5) 
    mock_sex.return_value = 'Undetermined'

    result = process_sample_xy(
        idxstats=sample_idxstats,
        x_id='chrX',
        y_id='chrY',
        male_kde=mock_kde,
        female_kde=mock_kde,
        threshold=0.95
    )

    # Check that Ry is NaN (not crashing)
    assert np.isnan(result['Ry'])
    assert result['Reads mapped to X'] == 0
    assert result['Reads mapped to Y'] == 0


########################################################
# test for ZW
########################################################
@patch('scims.classification.determine_sex_with_joint_posteriors')
@patch('scims.classification.compute_joint_posterior')
@patch('scims.classification.compute_coverage_ratio_rx')
def test_process_sample_zw(mock_rx, mock_posterior, mock_sex, sample_idxstats, mock_kde):
    """
    Test ZW logic. 
    Note: We reuse the same dataframe but treat 'chrX' as Z and 'chrY' as W for simplicity.
    """
    mock_rx.return_value = 1.2
    mock_posterior.return_value = (0.8, 0.2)
    mock_sex.return_value = 'Female'

    result = process_sample_zw(
        idxstats=sample_idxstats,
        z_id='chrX', # Treating chrX as Z (reads=50, len=500)
        w_id='chrY', # Treating chrY as W (reads=10, len=250)
        male_kde=mock_kde,
        female_kde=mock_kde,
        threshold=0.95
    )

    # Manual Calc Check:
    # Z_coverage = z_count / z_len = 50 / 500 = 0.1
    # W_coverage = w_count / w_len = 10 / 250 = 0.04
    # chr1_coverage = chr1_count / chr1_len = 200 / 1000 = 0.2
    # chr2_coverage = chr2_count / chr2_len = 100 / 1000 = 0.1
    # Rz Formula: (Z_coverage/chr1_coverage) + (Z_coverage/chr2_coverage) = 0.1/0.2 + 0.1/0.1 = 0.5 + 1 = 1.5
    # Rw Formula:(W_reads / total_zw) = (10 / 60) * (500 / 250) = 0.1666... * 2 = 0.333...
    expected_rz =  0.1/0.2 + 0.1/0.1
    expected_rw =(10 / 60) * (500 / 250)

    assert result['Rz'] == 1.2
    assert result['Rw'] == pytest.approx(expected_rw)
    assert result['SCiMS predicted sex'] == 'Female'
    assert result['Reads mapped to Z'] == 50
    assert result['Reads mapped to W'] == 10
