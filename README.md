# Ld1S_RiboMethSeq
RiboMeth-seq analysis pipeline for *Leishmania donovani* Ld1S

This repository contains a bash-based workflow and associated scripts used to process RiboMeth-seq data, map reads to Ld1S rRNA, compute RiboMeth-related scores, and generate site-level outputs suitable for differential analysis between experimental conditions.

## Requirements

The pipeline was developed and tested using:

- smalt **v0.7.6**
- samtools **v1.9** installed via conda (bioconda)
- bedtools **v2.30.0**
- Python **v3.10.9** & dependencies: numpy, pandas, biopython
- Rscript **v4.5.2**
- bash (Linux environment)

Paths to executables and reference files are defined at the top of `LD_RMS_Analysis_Pipeline.sh` and can be edited as needed.

## Input data

The pipeline starts from raw paired-end FASTQ files (`.fastq.gz` or `.fq.gz`).

Supported naming patterns include:
- `*_R1.fastq.gz` / `*_R2.fastq.gz`
- `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz`
- equivalent `.fq.gz` variants

All FASTQ files must be placed in a directory named `FASTQ/` in the working directory.

## Running the pipeline

1. Clone the repository and move into it:
```bash
git clone https://github.com/mikaOlami/Ld1S_RiboMethSeq.git
cd Ld1S_RiboMethSeq
```
2. Create a working directory and place FASTQ files:
```mkdir FASTQ
# copy raw FASTQ files into FASTQ/
```
3. Run the pipeline from the working directory:
```
chmod +x scripts/LD_RMS_Analysis_Pipeline.sh # only once
./scripts/LD_RMS_Analysis_Pipeline.sh
```

The script operates relative to the current working directory.

## Pipeline overview

For each sample, the bash pipeline performs the following steps:

1. Mapping of paired-end reads to the Ld1S rRNA reference using SMALT
2. Filtering and BAM generation using samtools
3. Conversion of BAM to sorted BED format using bedtools
4. Extraction of:
- 5′ initiation counts
- 3′ end counts
5. Calculation of RiboMeth-related scores across samples

Multiple samples are processed in parallel, while all steps within each sample are executed sequentially.

## Output

During execution, the pipeline creates the following directories in the working directory:
* Bams/ - rRNA-aligned BAM files
* Beds/ - sorted BED files
* Logs/ - per-step and global log files

Additional per-sample outputs include:
* *.init - initiation (5′) count files
* *.3p - 3′ end count files
* CSV files containing site-level RiboMeth-related scores

## Differential analysis criteria

Sites are considered significantly regulated if they meet all of the following criteria:
* Mean control C score ≥ 0.2
* Fold change ≥ 1.2 (up) or ≤ 1/1.2 (-1.2) (down)
* Two-sided t-test p-value < 0.05

## Reference files

All reference files used by the pipeline are provided in the `DB/` directory, including the Ld1S rRNA FASTA sequence, genome size file, and pre-built SMALT index.

Details on reference construction and file contents are provided in `DB/README_DB.md`.
