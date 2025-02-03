import unittest
import os
import tempfile
from scims.helpers import load_training_data

class TestHelpers(unittest.TestCase):
    def test_load_training_data(self):
        # Create a temporary training data file with expected headers.
        content = (
            "Run\tactual_sex\tactual_sex_zw\tSCiMS sample ID\tRx\tRy\n"
            "sample1\tmale\tfemale\tsample1\t0.5\t0.6\n"
            "sample2\tmale\tfemale\tsample2\t0.7\t0.8\n"
        )
        with tempfile.NamedTemporaryFile('w+', delete=False) as tmpfile:
            tmpfile.write(content)
            tmpfile.flush()
            tmp_filename = tmpfile.name
        
        try:
            df = load_training_data(tmp_filename)
            self.assertFalse(df.empty)
            # Now check for one or more of the expected columns
            for col in ["Run", "actual_sex", "actual_sex_zw", "SCiMS sample ID", "Rx", "Ry"]:
                self.assertIn(col, df.columns)
        finally:
            os.remove(tmp_filename)

if __name__ == '__main__':
    unittest.main()
