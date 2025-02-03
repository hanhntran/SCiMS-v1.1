import unittest
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
from scims.classification import process_sample_xy

class TestClassification(unittest.TestCase):
    def test_process_sample_xy(self):
        # Create a mock idxstats DataFrame
        data = {
            0: [1000, 1000, 1000],
            1: [500, 100, 50]
        }
        idxstats = pd.DataFrame(data, index=["chrX", "chrY", "chr1"])
        
        # Create a simple KDE model for testing using 3 data points.
        # This avoids a singular covariance matrix error.
        sample_data = np.array([
            [0.5, 1.5, 2.5],
            [0.3, 1.2, 0.8]
        ])
        kde = gaussian_kde(sample_data)
        
        # Run the classification function
        result = process_sample_xy(
            idxstats=idxstats,
            x_id="chrX",
            y_id="chrY",
            male_kde=kde,
            female_kde=kde,
            threshold=0.95
        )
        
        # Validate that the result contains expected keys
        expected_keys = [
            'Rx', 'Ry', 'Total reads mapped', 'Reads mapped to X',
            'Reads mapped to Y', 'Posterior probability of being male',
            'Posterior probability of being female', 'SCiMS predicted sex'
        ]
        for key in expected_keys:
            self.assertIn(key, result)

if __name__ == '__main__':
    unittest.main()
