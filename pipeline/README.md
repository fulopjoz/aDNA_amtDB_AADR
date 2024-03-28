# Creating the corrected markdown file content

md_content = """

## Mapping FASTQ Files to Reference Fasta with BWA and Samtool

This manual describes the process of mapping FASTQ files (`.fastq.gz`) to a reference fasta file (`rcrs.fa`) using BWA and Samtools. It assumes the FASTQ files are either single-end sequenced or paired reads already merged with sequencing adapters removed. This is a common format for published sequencing data, but variations exist.

## Prerequisites

- BWA
- Samtools
- Python (for specific scripts mentioned)
- `FilterUniqueSAMCons.py` script from Mayer Kircher's 2010 paper

## Process Overview

1. **Mapping reads to reference**: Using BWA for initial alignment, then Samtools for converting SAM to BAM format, sorting, and filtering.
2. **Removing sequence duplicates**: Two approaches are provided, one using Samtools directly and another using a Python script from Mayer Kircher's 2010 paper.
3. **Optional filtering**: Based on read length and identity percentage. This step requires a specific unpublished Python script.
4. **Building consensus sequences**: Two methods are discussed, one using Samtools and another using ANGSD, with the need to post-process the header of the output file.

## Detailed Steps

### Mapping Reads to Reference

```bash
bwa aln -l 16500 -n 0.01 -o 2 -t 4 rcrs.fa fastq.gz | bwa samse rcrs.fa - fastq.gz| samtools view -h -Su -F 4 - | samtools sort -T temp_file name -o output_file_name.bam
```

This command performs the mapping of reads in `fastq.gz` to the reference genome `rcrs.fa` using BWA, with specific parameters for alignment. The output is piped through Samtools to convert to BAM format, filter, and sort.

### Removing Sequence Duplicates

Two approaches are suggested:

1. Using Samtools:

```bash
samtools rmdup -S output_file_name.bam output_file_name_unique.bam
```

1. Using the `FilterUniqueSAMCons.py` script:

```bash
samtools view -F 4 -h output_file_name.bam | FilterUniqueSAMCons.py | samtools view -h -Su - > output_file_name_unique.bam
```

### Optional Filtering (if script is available)

```bash
samtools calmd output_file_name_unique.bam rcrs.fa | percidentity_threshold.py 0.9 35 short.txt | samtools view -bS - > output_file_name_unique_90prc.bam
```

This step filters out reads shorter than 35bp or with less than 90% identity to the reference. 

### Building Consensus Sequences

Two methods are presented:

- Using Samtools:

```bash
samtools consensus -m simple -d 3 --min-MQ 30 --min-BQ 30 -c 0.5 output_file_name_unique.bam -o output_consensus_file_name.fa
```

- Using ANGSD:

```bash
angsd -doFasta 2 -doCounts 1 -i output_file_name_unique.bam -minQ 30 -minMapQ 30 -setMinDepth 3 -out output_consensus_file_name
gunzip output_consensus_file_name.fa.gz
```

### Post-Processing

Replace the header in the consensus fasta file to include the sample name:

```bash
sed -i "s/>/>$(basename output_consensus_file_name.fa .fa)/" output_consensus_file_name.fa
```

## Conclusion

This manual outlines the steps for mapping sequencing reads to a reference genome and generating a consensus sequence. The process involves using widely-used bioinformatics tools and can be adapted based on specific requirements or data formats.
"""

## Saving the corrected content to a markdown file

file_path = '/mnt/data/manual.md'
with open(file_path, 'w') as file:
    file.write(md_content.strip())

file_path
