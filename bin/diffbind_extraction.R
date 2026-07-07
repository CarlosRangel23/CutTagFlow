#!/usr/bin/env Rscript

library(DiffBind)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: diffbind_extraction.R <label> <histone_mark> <comma_separated_peaks> <comma_separated_bams> <min_overlap>")

lbl         <- args[1]
mark        <- args[2]
peak_files  <- unlist(strsplit(args[3], ","))
bam_files   <- unlist(strsplit(args[4], ","))
min_overlap <- as.integer(args[5])

n_samples <- length(peak_files)
if (n_samples < 1) {
  stop("Error: No peak files provided to DiffBind.")
}

message(paste("Processing DiffBind for Mark:", mark, "| Label:", lbl, "| Min Overlap:", min_overlap))

# 1. Dynamically build the Sample Sheet
data_for_diffbind <- lapply(1:n_samples, function(i) {
  pk_file   <- peak_files[i]
  bam_file  <- bam_files[i]
  sample_id <- gsub("\\.(narrowPeak|broadPeak)$", "", basename(pk_file))
  caller_val <- if (grepl("\\.narrowPeak$", pk_file)) "narrowpeak" else "broadpeak"
  
  data.frame(
    SampleID   = sample_id,
    bamReads   = bam_file,
    Peaks      = pk_file,
    PeakCaller = caller_val,
    stringsAsFactors = FALSE
  )
})

df_diffbind <- do.call(rbind, data_for_diffbind)

# 2. Initialize DiffBind Object with the user-defined minOverlap
dba_obj <- dba(sampleSheet = df_diffbind, minOverlap = min_overlap)
dba_obj <- dba.count(dba_obj, summits = FALSE)

# 3. Retrieve the Consensus Peakset
consensus_peaks <- dba.peakset(dba_obj, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)

consensus_df <- consensus_peaks %>%
  mutate(peak_id = paste0(CHR, ":", START, "-", END))
 
consensus_bed <- data.frame(
        CHR   = consensus_df$CHR,
        START = as.integer(consensus_df$START - 1), 
        END   = as.integer(consensus_df$END),
        NAME  = consensus_df$peak_id,
        stringsAsFactors = FALSE
    )
# 4. Save Outputs using structured naming conventions
out_prefix <- paste0("DiffBind_", mark, "_", lbl, "_minOverlap", min_overlap)

# Save the full RData environment/object for downstream differential analysis
save(dba_obj, consensus_peaks, file = paste0(out_prefix, ".RData"))
write.table(consensus_bed, file = paste0(out_prefix, "consensus_bed.bed"), quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = FALSE)

if (!is.null(consensus_peaks) && nrow(consensus_peaks) > 0) {
  write_csv(consensus_peaks, file = paste0(out_prefix, "_peakset.csv"))
  message(paste("Successfully saved peakset with", nrow(consensus_peaks), "consensus peaks."))
} else {
  message("Warning: No consensus peaks passed the minOverlap threshold. CSV not generated.")
}