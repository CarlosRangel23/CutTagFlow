#!/usr/bin/env Rscript

library(DiffBind)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: consensus_decay.R <label> <histone_mark> <comma_separated_peaks> <comma_separated_bams>")

lbl          <- args[1]
mark         <- args[2]
peak_files   <- unlist(strsplit(args[3], ","))
bam_files    <- unlist(strsplit(args[4], ","))

n_samples <- length(peak_files)
if (n_samples < 2) {
  message(paste("Skipping:", lbl, "-", mark, "(Less than 2 samples)"))
  quit(status = 0)
}

# Dynamically build the Sample Sheet for DiffBind including bamReads
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

df_diffbind  <- do.call(rbind, data_for_diffbind)
overlap_rate <- numeric(n_samples)

# Loop through minOverlap thresholds to calculate the decay curve
for (i in 1:n_samples) {
  dba_obj_i <- dba(sampleSheet = df_diffbind, minOverlap = i)
  consensus_peaks <- dba.peakset(dba_obj_i, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
  overlap_rate[i] <- if (is.null(consensus_peaks)) 0 else nrow(consensus_peaks)
}

# Prepare data frame for plotting (Filtering for N >= 2)
n_seq <- 2:n_samples
overlap_rate_filtered <- overlap_rate[n_seq]

plot_df <- data.frame(
  N = factor(n_seq),
  Peaks = overlap_rate_filtered
)

fill_color <- if(lbl == "withDups") "#2980b9" else "#c0392b"
edge_color <- if(lbl == "withDups") "#1f5f8a" else "#8e281e"

# Generate the bar plot using ggplot2
p <- ggplot(plot_df, aes(x = N, y = Peaks)) +
  theme_minimal(base_size = 12) +
  geom_bar(stat = "identity", fill = fill_color, color = edge_color, width = 0.7, alpha = 0.85) +
  geom_text(
    aes(label = format(Peaks, big.mark = ",")),
    vjust = -0.5, fontface = "bold", color = "#2c3e50", size = 2.6
  ) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = paste("Consensus Peaks Decay:", mark),
    subtitle = paste("Condition:", lbl, "| Displaying Overlaps from N = 2 to", n_samples),
    x = "Minimum Overlapping Samples (N)",
    y = "Total Consensus Peaks"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13, color = "#2c3e50", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10, color = "#7f8c8d", margin = margin(b = 12)),
    axis.title = element_text(face = "bold", size = 11, color = "#2c3e50"),
    axis.text = element_text(color = "#34495e", size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray93"),
    panel.grid.minor = element_blank()
  )

# Output filename is relative for Nextflow's working directory architecture
png_name <- paste0("Decay_Plot_", mark, "_", lbl, ".png")
ggsave(filename = png_name, plot = p, width = 8, height = 5.5, dpi = 150)