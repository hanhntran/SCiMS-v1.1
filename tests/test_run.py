import unittest
import os
import tempfile
import subprocess

class TestRun(unittest.TestCase):
    def test_full_run(self):
        # Create temporary files for testing in a temporary directory.
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create metadata file
            metadata_path = os.path.join(tmpdir, "metadata.txt")
            with open(metadata_path, "w") as f:
                f.write("sample-id\tvalue\nsample1\t100\n")
    
            # Create scaffolds file
            scaffolds_path = os.path.join(tmpdir, "scaffolds.txt")
            with open(scaffolds_path, "w") as f:
                f.write("chrX\nchrY\nchr1\n")
            
            # Create master file and idxstats file
            master_file_path = os.path.join(tmpdir, "master.txt")
            idxstats_path = os.path.join(tmpdir, "sample1.idxstats")
            with open(idxstats_path, "w") as f:
                f.write("chrX\t1000\t500\nchrY\t1000\t100\nchr1\t1000\t50\n")
            with open(master_file_path, "w") as f:
                f.write(idxstats_path + "\n")
            
            # Create training data file with expected columns
            training_data_path = os.path.join(tmpdir, "training_data.txt")
            with open(training_data_path, "w") as f:
                f.write(
                    "Run\tactual_sex\tactual_sex_zw\tSCiMS sample ID\tRx\tRy\n"
                    "sample1\tmale\tfemale\tsample1\t0.5\t0.6\n"
                    "sample2\tmale\tfemale\tsample2\t0.7\t0.8\n"
                )
            
            output_file = os.path.join(tmpdir, "output.txt")
            
            # Run the command-line tool (assuming it's installed as "scims")
            result = subprocess.run([
                "scims",
                "--metadata", metadata_path,
                "--scaffolds", scaffolds_path,
                "--master_file", master_file_path,
                "--training_data", training_data_path,
                "--system", "XY",
                "--homogametic_id", "chrX",
                "--heterogametic_id", "chrY",
                "--output", output_file
            ], capture_output=True, text=True)
            
            # Check that the command ran without errors
            self.assertEqual(result.returncode, 0, msg=f"stderr: {result.stderr}")
            self.assertTrue(os.path.exists(output_file))

if __name__ == '__main__':
    unittest.main()
