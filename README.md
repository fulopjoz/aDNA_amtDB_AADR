
# aDNA Comparative Analysis

## Description
The "aDNA Comparative Analysis" project is designed for the comparative analysis of ancient DNA (aDNA) sequences. It integrates various datasets, performs haplogroup classification, and utilizes the MitoPathoTool for mitochondrial DNA analysis. The project structure includes directories for data, source code, notebooks, output results, and more, providing a comprehensive toolkit for researchers in the field of genetic history and evolution.

## Installation

To set up the project, ensure you have Python installed on your system. Then, clone the project and install the required Python packages:

```bash
git clone <repository-url>
cd aDNA_Comparative_Analysis
pip install -r requirements.txt
```

Activate the virtual environment if you are using one:

```bash
source <env-name>/bin/activate  # for Unix/Linux
<env-name>\Scripts\activate     # for Windows
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

## Contributing

Contributions to the "aDNA Comparative Analysis" project are welcome. Please fork the repository, make your changes, and submit a pull request for review.

For major changes, please open an issue first to discuss what you would like to change.

## License

[MIT License](LICENSE.md) or specify another license here.
