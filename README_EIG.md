sd# EIGENSOFT 7.2.1 Setup and Execution Guide

## Introduction
This manual guides through setting up and running EIGENSOFT version 7.2.1, which can be downloaded from:
[https://github.com/DReichLab/EIG/archive/v7.2.1.tar.gz](https://github.com/DReichLab/EIG/archive/v7.2.1.tar.gz).

## Preparations
1. **Installation**: Ensure EIGENSOFT is installed on your system.
2. **Dependencies**: Install necessary Perl libraries and other dependencies as listed in the documentation.
To run `eigenstrat` using the EIGENSOFT package, you need to follow several steps, including installing dependencies, compiling the software, and preparing the input files. Below are detailed instructions:

### 1. Install Dependencies
Ensure you have the required dependencies installed. This includes `libgsl-dev`, `libopenblas-dev`, and `liblapack-dev`.

```sh
sudo apt-get update
sudo apt-get install libgsl-dev libopenblas-dev liblapack-dev liblapacke-dev
```

### 2. Clone the EIGENSOFT Repository

```sh
git clone https://github.com/DReichLab/EIG
or
wget https://github.com/DReichLab/EIG/archive/v7.2.1.tar.gz
tar -xvzf v7.2.1.tar.gz
cd EIG-7.2.1
or 
cd EIG
```

### 3. Compile the Software
Navigate to the `src` directory and compile the software.

change this line in Makefile

override LDLIBS += -lgsl -lopenblas -llapacke -llapack -lm -lpthread

```sh
cd src
make clean
make install
```
if problem try: 
export PATH=~/projects/aDNA_Comparative_Analysis/eigenstrat/bin:$PATH
source ~/.bashrc

also we need to change the command paths in files bin/s,artpca.perrl, bin/smarteigenstrat.perl
```bash
$command = "/home/dodo/projects/aDNA_Comparative_Analysis/EIG-7.2.1/bin/smartpca"; # Path to smartpca
$cmd = "/home/dodo/projects/aDNA_Comparative_Analysis/EIG-7.2.1/bin/smarteigenstrat -p $outfilename.par >$logfilename"; # Path to smarteigenstrat
```

### 4. Prepare Input Files
Ensure you have the necessary input files in the required formats:

- `example.geno` - Genotype file
- `example.snp` - SNP file
- `example.ind` - Individual file
- `example.pca.evec` - PCA output file (if needed)
- `example.chisq.par` - Parameter file for `smarteigenstrat`

### 5. Create Parameter File for `smartpca`
Create a parameter file for `smartpca` (`example.pca.par`):

```sh
cat > ../EIGENSTRAT/example.pca.par << EOL
genotypename: ../EIGENSTRAT/example.geno
snpname: ../EIGENSTRAT/example.snp
indivname: ../EIGENSTRAT/example.ind
evecoutname: ../EIGENSTRAT/example.pca.evec
evaloutname: ../EIGENSTRAT/example.eval
altnormstyle: NO
numoutevec: 5
numoutlieriter: 5
numoutlierevec: 10
outliersigmathresh: 6
qtmode: 0
EOL
```

### 6. Run `smartpca`
Run `smartpca` to perform principal component analysis.

```sh
./smartpca.perl -i ../EIGENSTRAT/example.geno -a ../EIGENSTRAT/example.snp -b ../EIGENSTRAT/example.ind -k 5 -o ../EIGENSTRAT/example.pca -p ../EIGENSTRAT/example.plot -e ../EIGENSTRAT/example.eval -l ../EIGENSTRAT/example_pca.log
```

### 7. Create Parameter File for `smarteigenstrat`
Create a parameter file for `smarteigenstrat` (`example.chisq.par`):

```sh
cat > ../EIGENSTRAT/example.chisq.par << EOL
genotypename: ../EIGENSTRAT/example.geno
snpname: ../EIGENSTRAT/example.snp
indivname: ../EIGENSTRAT/example.ind
pcaname: ../EIGENSTRAT/example.pca.evec
outputname: ../EIGENSTRAT/example.chisq
numpc: 5
qtmode: NO
EOL
```

### 8. Run `smarteigenstrat`
Run `smarteigenstrat` to perform the EIGENSTRAT stratification correction.

```sh
./smarteigenstrat.perl -i ../EIGENSTRAT/example.geno -a ../EIGENSTRAT/example.snp -b ../EIGENSTRAT/example.ind -q NO -p ../EIGENSTRAT/example.pca.evec -k 5 -o ../EIGENSTRAT/example.chisq -l ../EIGENSTRAT/example_eigenstrat.log
```



### Summary
1. **Install dependencies** using `apt-get`.
2. **Clone the repository**.
3. **Compile the software**.
4. **Prepare input files**.
5. **Create and configure parameter files**.
6. **Run `smartpca`**.
7. **Create parameter file for `smarteigenstrat`**.
8. **Run `smarteigenstrat`**.
9. **Debug if necessary**.

By following these steps, you should be able to run `eigenstrat` using the EIGENSOFT package successfully.sudo 



## Configuration
Update the Perl scripts to include the correct paths to the executable binaries:

```bash
$command = "/home/dodo/projects/aDNA_Comparative_Analysis/EIG-7.2.1/bin/smartpca"; # Path to smartpca
$cmd = "/home/dodo/projects/aDNA_Comparative_Analysis/EIG-7.2.1/bin/smarteigenstrat -p $outfilename.par >$logfilename"; # Path to smarteigenstrat
```
Directory Structure

Create and maintain a structured directory:

    bin/ for executable scripts.
    EIGENSTRAT/ for data files.

Input Data

Prepare your data files:

    Genotype data (.geno)
    SNP information (.snp)
    Sample information (.ind)

Running the Analysis
Principal Component Analysis (PCA)

Execute PCA using:

bash

./smartpca.perl -i ../EIGENSTRAT/example.geno -a ../EIGENSTRAT/example.snp -b ../EIGENSTRAT/example.ind -o ../EIGENSTRAT/example.pca -p ../EIGENSTRAT/example.plot -e ../EIGENSTRAT/example.eval -l ../EIGENSTRAT/example.log

Stratification Correction

Perform stratification correction using:

bash

./smarteigenstrat.perl -i ../EIGENSTRAT/example.geno -a ../EIGENSTRAT/example.snp -b ../EIGENSTRAT/example.ind -p ../EIGENSTRAT/example.pca -o ../EIGENSTRAT/example.chisq -l ../EIGENSTRAT/example.log




# EIGENSOFT: Differences between EIGENSTRAT, SMARTPCA, and POPGEN

## Overview

EIGENSOFT is a software package designed for analyzing population structure and stratification in genotype data. The key components of the package include programs such as EIGENSTRAT, SMARTPCA, and several tools under the POPGEN directory. Each of these components serves a distinct purpose, has its own method of operation, and comes with specific advantages and limitations.

## EIGENSTRAT

### Purpose
EIGENSTRAT is used for stratification correction in genetic association studies. It adjusts for population structure by computing principal components and using them to correct for stratification in association tests.

### How It Works
- **Input**: Genotype data, SNP file, individual file.
- **Output**: Principal components, corrected association statistics.
- **Process**:
  1. Compute principal components from genotype data.
  2. Use these principal components to correct for stratification when testing for associations between genotypes and phenotypes.

### Positives
- Corrects for population stratification, reducing false positives in association studies.
- Supports multiple file formats.

### Negatives
- Assumes that all individuals are unrelated.
- May not work well with datasets containing a small number of markers.

### Restrictions
- Primarily aimed at genome-scan datasets with at least 100,000 markers.

## SMARTPCA

### Purpose
SMARTPCA performs Principal Components Analysis (PCA) on genotype data to identify population structure.

### How It Works
- **Input**: Genotype data, SNP file, individual file.
- **Output**: Eigenvectors, eigenvalues, and various statistics.
- **Process**:
  1. Runs PCA on the input genotype data.
  2. Outputs eigenvectors and eigenvalues.
  3. Optionally removes outliers and computes additional statistics.

### Positives
- Identifies and visualizes population structure.
- Supports multi-threading for improved performance.
- Offers a fast mode for large datasets.

### Negatives
- Fast mode has limited functionality (e.g., no outlier removal, no significance testing).
- Assumes that samples are unrelated, though a small number of cryptically related individuals can be handled.

### Restrictions
- Issues with .ped files containing QTL data due to a limit on the number of population values.
- Multi-threading is not supported in fast mode.

## POPGEN

### Purpose
POPGEN includes various tools for analyzing population structure, such as `ploteig`, `twstats`, and `smartrel`.

### Key Tools
- **ploteig**: Constructs plots of the top principal components.
- **twstats**: Computes Tracy-Widom statistics to evaluate the significance of principal components.
- **smartrel**: Identifies related samples while accounting for population structure.

### How They Work
- **ploteig**: Uses gnuplot to visualize principal components.
- **twstats**: Uses statistical methods to determine the significance of each principal component.
- **smartrel**: Computes a correlatioVittles yellow net. This is the car. origins. I. Animal input in the snape. converted smart PC in. Snipping direct evolve data Savitri Snape and. individuals. I don't know. Humanity genes. Still 10. opened my python. file file. Subaru. that finally. mom. Totento. Cheating Balashan stop bottom Sultan like File the attention. smart PC. Subaru. file a store path file. Time was since a little treasure came in the weather. That's right. Hmm. Mm-hmm. Genotype genotypes. individual. smart. name Gen file. Hey, thank you. Dental email. The I am imagined just tell not if Tom Jenifiable threatened preventing it will Stimulus pity a remix because. it's not one. echo. What's a? good? good, good, good. to be I. That's a good job. Snip. type. And you don't type name, but I'm in. with named and then I sit and Female male. female, male. I thought I lost it in 10 poppies. Double double bottom, but we should get smart smart. PC 's. Then. I. Then start Stand there. There's a bone density or. a Okay. There's a steam against the transducers. chair. Yeah, Judas Peter Jacn matrix, removes top eigenvectors, and identifies related individuals based on correlation coefficients.

### Positives
- Provides comprehensive tools for detailed population structure analysis.
- `ploteig` and `twstats` enhance the interpretation of PCA results.
- `smartrel` helps in identifying and managing related individuals in structured populations.

### Negatives
- Some tools require additional software (e.g., gnuplot for `ploteig`).
- `twstats` is not suitable for datasets with ancestry-informative markers due to admixture-LD.

### Restrictions
- `twstats` should not be used on datasets with ancestry-informative markers.
- Visualization and plotting tools require external dependencies.

## Conclusion

Each component of the EIGENSOFT package serves a unique role in the analysis of population structure and stratification in genetic data. While EIGENSTRAT focuses on correcting for stratification in association studies, SMARTPCA is used to identify and visualize population structure through PCA. The POPGEN suite provides additional tools for plotting, statistical analysis, and relatedness detection. Understanding the strengths and limitations of each component helps in choosing the appropriate tool for specific genetic analysis tasks.

