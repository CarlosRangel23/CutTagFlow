# CutTagFlow

CutTagFlow is a reproducible and scalable Nextflow pipeline designed for the processing and analysis of CUT&Tag sequencing data, from raw reads to high-quality peak calling and downstream quality control.

This pipeline is optimized for histone modification profiling, supporting both:

- **Narrow peaks**: H3K4me3, H3K27ac  
- **Broad domains**: H3K27me3, H3K9me3, H3K36me3, H3K4me1, H3K4me2  

CutTagFlow is designed to work **without spike-in controls or IgG samples**, simplifying experimental design while maintaining robust analytical performance.

To ensure optimal data interpretation, the pipeline performs analyses **both with duplicate reads retained and with duplicates removed**, allowing users to choose the most appropriate result depending on signal characteristics.

## Execution environment

Currently, CutTagFlow is supported **only in HPC environments** with:

- A **job scheduler / queue manager** (e.g. SLURM, PBS, SGE)
- Containerized execution via **Singularity or Apptainer** (Docker images can be used through these)

⚠️ **Local execution is not yet supported**, including:
- Interactive local servers
- Personal workstations or laptops  

Support for local environments may be added in future versions.

---

## Installation

CutTagFlow does not require a traditional installation. You can simply clone the repository:

```bash
git clone https://github.com/your_username/CutTagFlow.git
cd CutTagFlow
```

### Requirements

* Nextflow (DSL2 compatible versions)
* Singularity or Apptainer (recommended for HPC environments)


## Usage
### Input data format

Input samples must be specified in a samplesheet following the format defined in (data/example_samplesheet.csv)[https://github.com/CarlosRangel23/CutTagFlow/blob/main/data/example_samplesheet.csv]

### Required columns

* The **`sample` column must follow the format**:
  ```
  sampleName_histoneMark
  ```
* The `histone_mark` column must match the corresponding modification.

The samplesheet must contain the following columns:

```
sample,fastq1,fastq2,histone_mark
```

### Example

```csv
sample,fastq1,fastq2,histone_mark
BPES2_H3K27Ac,/path/to/file_r1.fastq.gz,/path/to/file_r2.fastq.gz,H3K27Ac
BPES2_H3K4me3,/path/to/file_r1.fastq.gz,/path/to/file_r2.fastq.gz,H3K4me3
BPES4_H3K27Ac,/path/to/file_r1.fastq.gz,/path/to/file_r2.fastq.gz,H3K27Ac
BPES4_H3K4me3,/path/to/file_r1.fastq.gz,/path/to/file_r2.fastq.gz,H3K4me3
```

## Running
Before running, load Nextflow:

```bash
module load apps/binapps/nextflow/25.10.4
```

### Running the pipeline

```bash
nextflow run main.nf --samplesheet data/example_samplesheet.csv [options]
```


### Parameters

The pipeline is structured into two main stages:

* **First steps**: preprocessing and alignment
* **Second steps**: downstream analysis

### Available parameters

* **fastqc**
* **trim**
* **fastqc_trim**
* **alignment**
* **filtering**
* **peaks**
* **diffbind**
* **coverage**

***

## Suggested improvement

You can simplify execution by adding a global parameter:

```bash
nextflow run main.nf --samplesheet data/example_samplesheet.csv --all true
```


## Summary

* No spike-in or IgG required
* Supports narrow and broad histone marks
* Runs analyses with and without duplicates
* Designed for HPC environments with containers
* Modular execution with flexible parameter control