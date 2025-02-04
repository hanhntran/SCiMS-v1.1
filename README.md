# `SCiMS` : Sex Calling in Metagenomic Sequencing
The goal of this software package is to provide an intuitive and accessible tool for identifying the sex of a host organism based on the alignment of metagenomic sequences.  

## Overview

Metagenomic sequencing data often contains a mix of host and non-host sequences. SCiMS salvages the reads mapping statistics that align to the host genome and uses them to identify the sex of the host organism. SCiMS leverages robust statistical methods to accurately determine the sex of the host, providing host sex information for downstream analyses.

## Requirements

- Python 3.9+
- numpy, pandas, scipy, setuptools
- (Optional) samtools for generating `.idxstats` files 

## Installation instructions

The simpliest installation works through the [conda](https://docs.conda.io/en/latest/miniconda.html) installer that can maintain different versions of Python on the same machine. 

```
# Create a new conda environment with Python 3.9
conda create -n scims python=3.9

# Activate the environment
conda activate scims

# Install SCiMS
pip install git+https://github.com/hanhntran/SCiMS-v1.1
```

To confirm that the instillation was successful, run:
```
scims -h
```

## Usage
`SCiMS` can be used on any alignment data, regardless of the platform used for sequencing or the aligner that generated the alignment file. 

```
scims --master_file <master_file.txt> --scaffolds <scaffolds.txt> --metadata <metadata_file.txt> --system XY --homogametic_id chrom_id --heterogametic_id chrom_id --id_column sample-id --output <output_file.txt>
```
| Option             | Description                                                                          |
|--------------------|----------------------------------------------------------------------------------|
| -h, --help         | Show this help message and exit                                                      |
| --master_file      | Path to the master file containing paths to the .idxstats files for each sample         |
| --scaffolds        | Path to the scaffolds.txt file containing the scaffolds of interest                     |
| --metadata         | Path to the metadata file containing the sample IDs                                   |
| --heterogametic_id | The ID of the heterogametic sex chromosome                                             |
| --homogametic_id   | The ID of the homogametic sex chromosome                                               |
| --system           | The sex determination system (XY or ZW)                                                |
| --output           | Path to the output file                                                                |
| --threshold        | The threshold for the sex calling algorithm (default: 0.95)                             |
| --id_column        | The column name of the sample ID in the metadata file                                  |
| --training_data    | If you have a training dataset, you can specify the path to the training data here        |

## Required input files

### `scaffolds.txt`
Since most assemblies include scaffolds representing other DNA than simply genomic (ex. mitochondrial), it is necessary to define what scaffolds we are interested in using for our analysis. This can be specified with a ```scaffolds.txt``` file. This is a single-column text file where each row is a scaffold ID. Here is an example, 
```
NC_000001.11
NC_000002.12
NC_000003.12
NC_000004.12
NC_000005.10
NC_000006.12
NC_000007.14
NC_000008.11
...
``` 

### `.idxstats files`
A .idxstats file can easily be created with samtools. If you have a .bam file of interest, fun the following commands to generate the .idxstats file:

```shell
samtools index <bam_file>
```

```shell
samtools idxstats <bam_file> > <prefix>.idxstats
```
### `master_file.txt`
To run SCiMS, you will need to create a master file. This file should contain the paths to the .idxstats files for each sample. 

Example of a master_file.txt:
```
path/to/idxstats/file/sample1.idxstats
path/to/idxstats/file/sample2.idxstats
path/to/idxstats/file/sample3.idxstats
path/to/idxstats/file/sample4.idxstats
```

### `metadata_file.txt`
A metadata file is required to run SCiMS. This file should contain at least one columns, `sample-id`. The `sample-id` column should contain the sample IDs that are present in the .idxstats file. 

Example:
```
sample-id	feature
sample1		A
sample2		B
sample3		C
sample4		D

```

## Example run
Example files can be found in the ```test_data``` folder.

Change path to the ```test_data``` folder and run the following command:
```
scims --scaffolds GRCh38_scaffolds.txt --master_file test_master_file.txt \
      --metadata metadata_file.txt \
      --system XY \
      --homogametic_id NC_000023.11 \
      --heterogametic_id NC_000024.10 \
      --id_column sample-id \
      --output test_output.txt
```

Output:
```
=================================================

    _|_|_|   _|_|_|  _|_|_|  _|      _|   _|_|_|  
    _|      _|         _|    _|_|  _|_|   _|        
    _|_|_|  _|         _|    _|  _|  _|   _|_|_|    
        _|  _|         _|    _|      _|       _|  
    _|_|_|   _|_|_|  _|_|_|  _|      _|   _|_|_|    
    
=================================================
SCiMS: Sex Calling in Metagenomic Sequencing
2025-02-04 12:30:14,018 - INFO - Using 'Run' as the sample ID column.
2025-02-04 12:30:55,747 - INFO - Results written to test_output.txt
```


