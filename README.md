<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/othneildrew/Best-README-Template">
    <img src="static/scims_logo.png" alt="Logo" width="450" height="200">
  </a>

  <h1 align="center">Sex Calling in Metagenomic Sequences</h1>

  <p align="center">
    An tool for identifying the sex of a host organism based on the alignment of metagenomic sequences.
    <br />
    <br />
    <a href="https://github.com/hanhntran/SCiMS-v1.1/issues/new?labels=bug&template=bug-report---.md">Report Bug</a>
    &middot;
    <a href="https://github.com/hanhntran/SCiMS-v1.1/issues/new?labels=enhancement&template=feature-request---.md">Request Feature</a>
  </p>
</div>



<!-- ABOUT THE PROJECT -->
## About The Project

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
pip install scims
```

To confirm that the instillation was successful, run:
```
scims -h
```

## Usage
`SCiMS` can be used on any alignment data, regardless of the platform used for sequencing or the aligner that generated the alignment file. 

```
scims --idxstats_file <sample.idxstats> \
        --scaffolds <scaffolds.txt> \
        --metadata <metadata_file.txt> \
        --system <XY or ZW> \
        --homogametic_id <chrom_id> \
        --heterogametic_id <chrom_id> \
        --id_column <sample-id> \
        --output <output_file.txt>
```
| Option             | Description                                                                          |
|--------------------|----------------------------------------------------------------------------------|
| -h, --help         | Show this help message and exit                                                      |
| --idxstats_file    | Path to the .idxstats file for the sample                                             |
| --scaffolds        | Path to the scaffolds.txt file containing the scaffolds of interest                     |
| --heterogametic_id | The ID of the heterogametic sex chromosome                                             |
| --homogametic_id   | The ID of the homogametic sex chromosome                                               |
| --system           | The sex determination system (XY or ZW)                                                |
| --output           | Path to the output file                                                                |
| --threshold [OPTIONAL]        | The threshold for the sex calling algorithm (default: 0.95)                             |
| --training_data [OPTIONAL]    | If you have a training dataset, you can specify the path to the training data here        |
| --multiple [OPTIONAL]    | If you want to run SCiMS on multiple samples, you can specify this option [True or False, default: False]         |
| --metadata_file [OPTIONAL]    | If you have a metadata file and would like to add SCiMS predicted sex to the metadata file, you can specify the path to the metadata file here         |
| --id_column [OPTIONAL]        | The column name of the sample ID in the metadata file                                  |
| --log [OPTIONAL]    | Path to log file         |
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

### Running SCiMS on a single sample
Change path to the ```test_data``` folder and run the following command:
```
scims --idxstats_file ./idxstats_files/S79F300.idxstats \
      --scaffolds GRCh38_scaffolds.txt \
      --system XY \
      --homogametic_id NC_000023.11 \
      --heterogametic_id NC_000024.10 \
      --output test_output.txt
```

Output log:
```
2025-03-06 00:22:34,576 - INFO - Log file created at: out/scims.log
2025-03-06 00:22:34,576 - INFO -  
=================================================
2025-03-06 00:22:34,576 - INFO - 
    _|_|_|   _|_|_|  _|_|_|  _|      _|   _|_|_|  
    _|      _|         _|    _|_|  _|_|   _|        
    _|_|_|  _|         _|    _|  _|  _|   _|_|_|    
        _|  _|         _|    _|      _|       _|  
    _|_|_|   _|_|_|  _|_|_|  _|      _|   _|_|_|    
    =================================================
2025-03-06 00:22:34,576 - INFO - SCiMS: Sex Calling in Metagenomic Sequencing
2025-03-06 00:22:34,576 - INFO - Version: 1.0.0
2025-03-06 00:22:34,576 - INFO - =================================================
2025-03-06 00:22:34,591 - INFO - Results written to out/S79F300_results.txt
```
Output file:
```
$ cat out/S79F300_results.txt
```

### Running SCiMS on multiple samples

```
scims   --idxstats_folder idxstats_files/  \
        --scaffolds GRCh38_scaffolds.txt \
        --homogametic_id NC_000023.11 \
        --heterogametic_id NC_000024.10 \
        --output_dir out \
        --metadata metadata_file.txt \
        --id_column sample-id \
        --log log.txt
```

Output log:
```
2025-03-06 00:29:09,830 - INFO - Log file created at: out/scims.log
2025-03-06 00:29:09,830 - INFO -  
=================================================
2025-03-06 00:29:09,830 - INFO - 
    _|_|_|   _|_|_|  _|_|_|  _|      _|   _|_|_|  
    _|      _|         _|    _|_|  _|_|   _|        
    _|_|_|  _|         _|    _|  _|  _|   _|_|_|    
        _|  _|         _|    _|      _|       _|  
    _|_|_|   _|_|_|  _|_|_|  _|      _|   _|_|_|    
    =================================================
2025-03-06 00:29:09,830 - INFO - SCiMS: Sex Calling in Metagenomic Sequencing
2025-03-06 00:29:09,830 - INFO - Version: 1.0.0
2025-03-06 00:29:09,830 - INFO - =================================================
2025-03-06 00:29:09,845 - INFO - Results written to out/S28M1000000_results.txt
2025-03-06 00:29:09,846 - INFO - Updated metadata with classification results written to out/metadata_with_classification.txt
2025-03-06 00:29:09,848 - INFO - Results written to out/S56F150_results.txt
2025-03-06 00:29:09,849 - INFO - Updated metadata with classification results written to out/metadata_with_classification.txt
2025-03-06 00:29:09,851 - INFO - Results written to out/S79F300_results.txt
2025-03-06 00:29:09,852 - INFO - Updated metadata with classification results written to out/metadata_with_classification.txt
2025-03-06 00:29:09,854 - INFO - Results written to out/S90M250_results.txt
2025-03-06 00:29:09,855 - INFO - Updated metadata with classification results written to out/metadata_with_classification.txt
```
