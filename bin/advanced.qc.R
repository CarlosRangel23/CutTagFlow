#!/usr/bin/env Rscript

# =========================================================================
# 1. LIBRARY LOADING & CONFIGURATION
# =========================================================================
# Suppress startup messages to keep the Nextflow log clean
suppressPackageStartupMessages({
  library(ATACseqQC)
  library(Rsamtools)
  library(GenomicAlignments)
  library(GenomicRanges)
  library(rtracklayer)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene) # Adjust depending on your reference genome
})

# Parse command-line arguments passed by the Nextflow module
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Missing arguments. Usage: advanced.qc.R <sample> <histone_mark> <label> <bam> <peak>")
}

sample_id    <- args[1]
histone_mark <- args[2]
label        <- args[3]
bam_file     <- args[4]
peak_file    <- args[5]

# Define the reference genome annotation database for TSSE
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# =========================================================================
# 2. FRiP SCORE CALCULATION (Fraction of Reads in Peaks)
# =========================================================================
message(paste("[QC]", sample_id, "-", label, ": Calculating FRiP score..."))

# Import peak regions (rtracklayer smoothly handles MACS3 narrowPeak/broadPeak formats)
peaks <- import(peak_file)

# Extract total mapped reads directly from the BAM index for high-speed processing
idx <- idxstatsBam(bam_file)
total_reads <- sum(idx$mapped)

# Count how many reads overlap with the defined peak regions
# Filter out secondary alignments and technical duplicates for accuracy
param <- ScanBamParam(which = peaks, flag = scanBamFlag(isSecondaryAlignment = FALSE, isDuplicate = FALSE))
reads_in_peaks <- countBam(bam_file, param = param)$records

# Compute the final FRiP ratio
frip_score <- ifelse(total_reads > 0, reads_in_peaks / total_reads, 0)

# =========================================================================
# 3. TSSE SCORE CALCULATION (TSS Enrichment Score)
# =========================================================================
message(paste("[QC]", sample_id, "-", label, ": Calculating TSSE score..."))

# Extract Transcription Start Sites (TSS) and resize regions to the exact 1bp coordinate
TSS <- genes(txdb)
TSS <- resize(TSS, width = 1, fix = "start")

# Load paired-end reads into the GalignmentPairs structure required by ATACseqQC
gal <- readGAlignmentPairs(bam_file, param = ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = FALSE)))

# Compute the flanking enrichment profile around the TSS regions
tsse_result <- TSSeScore(gal, TSS)

# Extract the single global numeric score for this sample
tsse_score <- tsse_result$TSSeScore

# =========================================================================
# 4. EXPORT METRICS TO CSV
# =========================================================================
# Structure the metric dataframe with clean English column headers
qc_data <- data.frame(
  sample       = sample_id,
  histone_mark = histone_mark,
  label        = label,
  total_reads  = total_reads,
  peaks_count  = length(peaks),
  frip         = frip_score,
  tsse         = tsse_score,
  stringsAsFactors = FALSE
)

# Generate the dynamic output filename matching the Nextflow output block pattern
output_name <- paste0(sample_id, ".", label, ".qc_metrics.csv")
write.csv(qc_data, file = output_name, row.names = FALSE)

message(paste("[QC] Successfully completed for:", sample_id))