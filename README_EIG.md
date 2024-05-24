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

