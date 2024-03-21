
# aDNA Comparative Analysis

## Description

The "aDNA Comparative Analysis" project is aimed at comparing the ancient
DNA (aDNA) sequences between the AmtDB and the AADR. The primary objective is
to identify sequences present in the AADR Reich fasta file but missing in the
AmtDB, and subsequently retrieve these sequences. In the initial stage, we
managed to retrieve 404 sequences, among which 8 sequences had duplicates.
The duplicates differed only in the ending of the sequences, where additional
'N's were appended.

The second stage of the project involves creating a parsing function to take
metadata from AADR, MitoPatho, and parse the metadata into a format suitable for AmtDB.
From 920 missing sequences I was able to retrieve metadata for 561 sequences and 
325 sequence IDs had match with mitopatho data.

Subsequent stages of the project will focus on the usage of a pipeline that
includes a manual for mapping fastq files to a reference genome, utilizing
tools such as BWA and Samtools.

Additionally, the project will explore the use of EIGENSOFT for further analyses.

This comprehensive toolkit includes directories for data, source code, 
notebooks, and output results, designed to facilitate the comparative analysis of aDNA sequences across different datasets.

## Installation

To set up the project, ensure you have Python installed on your system. Then, clone the project and install the required Python packages:

```bash
git clone [<repository-url>](https://github.com/fulopjoz/aDNA_amtDB_AADR.git)
cd aDNA_Comparative_Analysis
./setup_project.sh
or 
mamba env create -f aDNA_env.yml

```

Activate the virtual environment if you are using one:

```bash
conda activate aDNA_env       # for Conda
source aDNA_env/bin/activate  # for Unix/Linux
aDNA_env\Scripts\activate     # for Windows
```

## Project Structure

- `data/`: Contains all the input data files, including FASTA sequences and metadata CSVs.
- `src/`: Source code for data loading, processing, and analysis.
- `notebooks/`: Jupyter notebooks for interactive analysis.
- `output/`: Output directory for results, including processed haplogroups and analysis reports.
- `pipeline/`: Bash scripts and other tools for automating analysis workflows.
- `tests/`: Unit tests for the Python code.

## Usage

### Setting Up Your Environment

Before running any scripts, ensure your working directory is set to the project's root:

```python
import os
os.chdir('/path/to/aDNA_Comparative_Analysis')
```

### Loading and Processing Data

```python
from src.data_loading import load_data, load_mt_dataset
from src.data_processing import DataProcessor

# Load your datasets
meta_amtDB, ids_seq_fasta = load_data('data/amtDB/amtdb_metadata.csv', 'data/amtDB/amtdb_1621-samples_7f_a0pkh.fasta')
```

### Analyzing Haplogroups

Example commands for analyzing haplogroups and integrating with MitoPatho:

```bash
haplogrep3 classify --tree phylotree-rcrs@17.0 --in output/missing_sequences_AmtDB.fasta --out output/analysis_result.hsd --extend-report
python3 mitopatho/mitopatho.py -i output/processed_haplogroups.hsd -o output/mitopatho_output.txt
```

### Further Analysis

Refer to the Jupyter notebooks within the `notebooks/` directory for interactive analysis examples and additional instructions on using the project for aDNA analysis.
