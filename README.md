# CutTagFlow
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](./LICENSE)
[![Follow on Twitter/X](https://img.shields.io/badge/X-@rangelpelaezc-black?logo=x)](https://x.com/rangelpelaezc)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-carlosrangelpelaez-0A66C2?logo=linkedin)](https://www.linkedin.com/in/carlosrangelpelaez/)
[![ORCID](https://img.shields.io/badge/ORCID-0000--0001--7697--1696-A6CE39?logo=orcid)](https://orcid.org/0000-0001-7697-1696)


## Introduction
CutTagFlow is a reproducible and scalable Nextflow pipeline designed for the processing and analysis of CUT&Tag sequencing data, from raw reads to high-quality peak calling and downstream quality control.

This pipeline is optimized for histone modification profiling, supporting both:

- **Narrow peaks**: H3K4me3, H3K27Ac  
- **Broad domains**: H3K27me3, H3K9me3, H3K36me3, H3K4me1, H3K4me2  

CutTagFlow is designed to work **without spike-in controls or IgG samples**, simplifying experimental design while maintaining robust analytical performance. To ensure optimal data interpretation, the pipeline performs analyses **both with duplicate reads retained and with duplicates removed**, allowing users to choose the most appropriate result depending on signal characteristics.

### Workflow
```mermaid
flowchart TD

  A[📁 FASTQ Files]
  B[📊 FastQC + MultiQC]
  C[✂️ Read Trimming]
  D[📊 FastQC + MultiQC]
  E[🧬 Bowtie2 Alignment]
  F[🔍 BAM Filtering]

  A --> B
  A --> C
  C --> D
  C --> E
  E --> F

  G[📦 Duplicated BAMs]
  H[📦 Deduplicated BAMs]

  F --> G
  F --> H
  G --> I[🎯 MACS3 Peak Calling]
  H --> I
  I --> J[📈 DiffBind Consensus Peaks]
  H --> K[⚖️ Scale Factor Calculation]
  K --> L[📉 DeepTools Coverage Tracks]
  H --> M[🧬 GATK Variant Calling]
  
  classDef qc fill:#D6EAF8,stroke:#2471A3,stroke-width:2px,color:#000;
  classDef align fill:#D5F5E3,stroke:#1E8449,stroke-width:2px,color:#000;
  classDef peak fill:#FADBD8,stroke:#C0392B,stroke-width:2px,color:#000;
  classDef downstream fill:#FCF3CF,stroke:#B7950B,stroke-width:2px,color:#000;
  class A,B,C,D qc;
  class E,F,G,H align;
  class I,J peak;
  class K,L,M downstream;
```

### Execution environment
CutTagFlow is highly portable and supports execution in both high-performance computing (HPC) clusters and personal local environments through pre-configured profiles.

#### 1. HPC Environments
For large-scale processing, the pipeline natively integrates with:
- A **job scheduler / queue manager** (e.g., SLURM, PBS, SGE).
- Containerized execution via **Singularity/Apptainer**.

You may change some specifications dependending on your HPC system. To launch the pipeline in an HPC environment using SLURM, run:

```bash
NXF_VER=25.10.4 nextflow run main.nf -profile slurm
```

#### 2. Local Environments (Workstations & Laptops)
To run the pipeline locally for testing, development, or smaller datasets, you only need to have the following dependencies installed:

- **Java** (version 17 or later)
- **Nextflow** (version 23 or later)
- **Apptainer/Singularity**

To launch the pipeline in your local environment, run:

```bash
NXF_VER=25.10.4 nextflow run main.nf -profile going_merry
```

**Smart Resource Management:** By default, the `going_merry` profile is designed to dynamically cap CPU usage at a maximum of **70% of your machine's total capacity**, leaving at least 30% completely free so your computer remains smooth and responsive. If you wish to change this threshold, you can pass your desired maximum percentage directly via the command line using the --max_cpu_percent parameter. For example, to cap it at 50% or push it to 90%:

```bash
# Run using only 50% of your local CPUs
NXF_VER=25.10.4 nextflow run main.nf -profile going_merry --max_cpu_percent 50

# Run using 90% of your local CPUs
NXF_VER=25.10.4 nextflow run main.nf -profile going_merry --max_cpu_percent 90
```

### Installation

CutTagFlow does not require a traditional installation. You can simply clone the repository:

```bash
git clone https://github.com/CarlosRangel23/CutTagFlow.git
cd CutTagFlow
```

## Usage
### Input data format

Input samples must be specified in a samplesheet following the format defined in [data/example_samplesheet.csv](https://github.com/CarlosRangel23/CutTagFlow/blob/main/data/example_samplesheet.csv). Importantly, modify all options available and references needed for the pipeline in [conf/params.config](https://github.com/CarlosRangel23/CutTagFlow/blob/main/conf/params.config)

#### Required columns
The samplesheet must contain the following columns:

```
sample,fastq1,fastq2,histone_mark
```


* The **`sample` column must follow the format**:
  ```
  sampleName_histoneMark
  ```

* The `fastq1` column must contain path to r1 file.
* The `fastq2` column must contain path to r2 file.
* The `histone_mark` column must match the corresponding modification.


#### Example

```csv
sample,fastq1,fastq2,histone_mark
BPES2_H3K27Ac,/path/to/file_r1.fastq.gz,/path/to/file_r2.fastq.gz,H3K27Ac
BPES2_H3K4me3,/path/to/file_r1.fastq.gz,/path/to/file_r2.fastq.gz,H3K4me3
BPES4_H3K27Ac,/path/to/file_r1.fastq.gz,/path/to/file_r2.fastq.gz,H3K27Ac
BPES4_H3K4me3,/path/to/file_r1.fastq.gz,/path/to/file_r2.fastq.gz,H3K4me3
```

### Running the pipeline
Before running, load Nextflow if you are in a HPC cluster. Be aware that your cluster may have some guidance or specific rules for using Nextflow.

```bash
module load apps/binapps/nextflow/25.10.4
```

```bash
NXF_VER=25.10.4 nextflow run main.nf --samplesheet data/example_samplesheet.csv [options]
```

### Available parameters

- You can execute the entire pipeline with `--all true`

#### First steps (preprocessing and alignment)

- **first_steps**  
  Performs all preprocessing and alignment

- **fastqc**  
  Runs quality control on raw FASTQ files using FastQC. Runs also MultiQC.

- **trim**  
  Performs adapter and quality trimming of raw reads.

- **fastqc_trim**  
  Runs FastQC again after trimming to assess read quality improvement. Runs also MultiQC.

- **alignment**  
  Aligns reads to the reference genome and generates BAM files using bowtie2.

#### Second steps (downstream analysis)

- **second_steps**  
  Performs all downstream analyses

- **filtering**  
  Filters aligned reads (e.g. low quality, mitochondrial reads, etc.) to improve signal-to-noise ratio.

- **peaks**  
  Performs peak calling to identify enriched regions (supports both narrow and broad marks) using macs3.

- **diffbind**  
  Generates consensus peaks using DiffBind and according to the most wanted minimum overlap.

- **coverage**  
  Generates coverage tracks (e.g. bigWig) for visualization in genome browsers.

For this step, you need to add `--meta_spike` parameter. Example is defined here: [data/sample_metadata.txt](https://github.com/CarlosRangel23/CutTagFlow/blob/main/data/sample_metadata.txt).

- **variant**  
  Performs GATK variant calling

#### Example run

```bash
NXF_VER=25.10.4 nextflow run main.nf \
  --samplesheet data/example_samplesheet.csv \
  --first_steps true \
  -profile slurm
```

### Output

The pipeline generates:

- Quality control reports (FastQC and MultiQC)
- Trimmed FASTQ files
- Alignment files (BAM + indexes)
- Filtered BAM files (Duplicated and Deduplicated)
- Peak files (narrowPeak / broadPeak)
- Differential consensus peaks (DiffBind)
- Coverage tracks (bigWig)
- VCF file for each histone mark group
