import unittest
import pandas as pd
import os
import tempfile
from scims.utils import read_metadata, normalize_colname, find_sample_id_column

class TestUtils(unittest.TestCase):
    def test_normalize_colname(self):
        self.assertEqual(normalize_colname("Sample_ID"), "sampleid")
        self.assertEqual(normalize_colname("Gene-Name!"), "genename")

    def test_read_metadata(self):
        # Create a temporary metadata file
        content = "sample-id\tvalue\nsample1\t100\nsample2\t200\n"
        with tempfile.NamedTemporaryFile('w+', delete=False) as tmpfile:
            tmpfile.write(content)
            tmpfile.flush()
            tmp_filename = tmpfile.name
        try:
            df = read_metadata(tmp_filename)
            self.assertEqual(df.shape, (2, 2))
            self.assertIn("sample-id", df.columns)
        finally:
            os.remove(tmp_filename)

    def test_find_sample_id_column(self):
        # Create a simple DataFrame
        df = pd.DataFrame({
            "Sample-ID": ["sample1", "sample2"],
            "value": [10, 20]
        })
        sample_col = find_sample_id_column(df)
        self.assertEqual(sample_col.lower(), "sample-id".lower())

if __name__ == '__main__':
    unittest.main()
