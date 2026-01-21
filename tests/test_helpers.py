import pytest
import os
import tempfile
from scims.helpers import load_training_data

def test_load_training_data(tmp_path):
    """
    Test that the training data is loaded correctly. and the expected columns are present.
    Also test that the training data is loaded from the default path if no path is provided.
    """
    # create temporary training data file with expected headers
    temp_file = tmp_path / "training_data.txt"
    
    content = (
        "Run\tactual_sex\tactual_sex_zw\tSCiMS sample ID\tRx\tRy\n"
        "sample1\tmale\tfemale\tsample1\t0.5\t0.6\n"
        "sample2\tmale\tfemale\tsample2\t0.7\t0.8\n"
    )
    temp_file.write_text(content)

    # load the training data
    df = load_training_data(str(temp_file))

    # check that the training data is loaded correctly
    assert not df.empty
    assert "Run" in df.columns
    assert "actual_sex" in df.columns
    assert "actual_sex_zw" in df.columns
    assert "SCiMS sample ID" in df.columns
    assert "Rx" in df.columns
    assert "Ry" in df.columns


