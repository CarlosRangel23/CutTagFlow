#!/usr/bin/env Rscript

# Load required bioinformatic and plotting libraries
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(patchwork)

# 1. Discover all pipeline input files dynamically in the Nextflow execution directory
picard_files   <- list.files(pattern = "_picard")
frag_files     <- list.files(pattern = "_fragmentLen")
idxstats_files <- list.files(pattern = "_idxstats")

# Initialize global empty data frames to accumulate structured metrics
dup_results   <- data.frame()
frag_data     <- data.frame()
idxstats_data <- data.frame()

# =========================================================================
# STEP 1: PARSE PICARD DUPLICATION METRICS
# =========================================================================
for (f in picard_files) {
  sample_name        <- str_remove(basename(f), "_picard.*")
  lines              <- readLines(f)
  metrics_header_idx <- grep("## METRICS CLASS", lines)
  
  metrics <- read.table(f, skip = metrics_header_idx, nrows = 1, header = TRUE, sep = "\t", fill = TRUE)
  
  mapped_reads <- as.numeric(metrics$READ_PAIRS_EXAMINED[1])
  dup_rate     <- as.numeric(metrics$PERCENT_DUPLICATION[1]) * 100
  est_lib_size <- as.numeric(metrics$ESTIMATED_LIBRARY_SIZE[1])
  unique_reads <- mapped_reads * (1 - (dup_rate / 100))
  
  name_split <- str_split(sample_name, "_")[[1]]
  sampleid   <- name_split[1]
  histone    <- name_split[2]
  
  dup_results <- rbind(dup_results, data.frame(
    Sample = sample_name,
    Sample_ID = sampleid,
    Histone = histone,
    Mapped_Fragments = mapped_reads,
    Duplication_Rate_Pct = dup_rate,
    Estimated_Library_Size = est_lib_size,
    Unique_Fragments = unique_reads
  ))
}

# =========================================================================
# STEP 2: PARSE SAMTOOLS FRAGMENT LENGTHS
# =========================================================================
for (f in frag_files) {
  sample_name <- str_remove(basename(f), "_fragmentLen.*")
  name_split  <- str_split(sample_name, "_")[[1]]
  sampleid    <- name_split[1]
  histone     <- name_split[2]
  
  df <- read.table(f, header = FALSE, col.names = c("fragLen", "fragCount"))
  
  df <- df %>% mutate(
    Weight = fragCount / sum(fragCount),
    Sample = sample_name,
    Sample_ID = sampleid,
    Histone = histone
  )
  
  frag_data <- rbind(frag_data, df)
}

# =========================================================================
# STEP 3: PARSE SAMTOOLS IDXSTATS
# =========================================================================
for (f in idxstats_files) {
  sample_name <- str_remove(basename(f), "_idxstats.*")
  name_split  <- str_split(sample_name, "_")[[1]]
  sampleid    <- name_split[1]
  histone     <- name_split[2]
  
  df <- read.table(f, sep = "\t", stringsAsFactors = FALSE, col.names = c("chrom", "length", "mapped", "unmapped"))
  
  total_mapped  <- sum(df$mapped[df$chrom != "*"])
  chrM_mapped   <- sum(df$mapped[df$chrom %in% c("chrM", "MT", "M", "chrM_spikein")])
  chrX_mapped   <- sum(df$mapped[df$chrom == "chrX"])
  chrY_mapped   <- sum(df$mapped[df$chrom == "chrY"])
  chrEBV_mapped <- sum(df$mapped[df$chrom == "chrEBV"])
  
  idxstats_data <- rbind(idxstats_data, data.frame(
    Sample = sample_name,
    Sample_ID = sampleid,
    Histone = histone,
    Total_Mapped_Reads = total_mapped,
    ChrM_Reads = chrM_mapped,
    ChrM_Fraction = ifelse(total_mapped > 0, chrM_mapped / total_mapped, 0),
    ChrX_Fraction = ifelse(total_mapped > 0, chrX_mapped / total_mapped, 0),
    ChrY_Fraction = ifelse(total_mapped > 0, chrY_mapped / total_mapped, 0),
    ChrEBV_Fraction = ifelse(total_mapped > 0, chrEBV_mapped / total_mapped, 0)
  ))
}

# Consolidate metrics tables into a comprehensive flat text report
qc_summary_df <- merge(dup_results, idxstats_data, by = c("Sample", "Sample_ID", "Histone"))
write.table(qc_summary_df, file = "alignment_dup_qc_summary.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# =========================================================================
# STEP 4: GENERATE dup_plots.pdf
# =========================================================================
pdf("dup_plots.pdf", width = 12, height = 6)

# Page 1: Global barplot of duplication rates per sample (ATAC legacy style)
p_dup_bar <- ggplot(qc_summary_df, aes(Sample, Duplication_Rate_Pct)) +
  geom_col(fill = "steelblue") +
  theme_bw() +
  labs(title = "Duplication Rate Levels per Sample", x = "Sample", y = "Percent duplication (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_dup_bar)

# Page 2: Combined comparative panels by Histone mark (CUT&Tag style)
figA <- ggplot(qc_summary_df, aes(x = Histone, y = Duplication_Rate_Pct, fill = Histone)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = Sample_ID), width = 0.15, size = 3, alpha = 0.8) +
  theme_bw(base_size = 14) + labs(y = "Duplication Rate (%)", x = "", color = "Sample ID")

figB <- ggplot(qc_summary_df, aes(x = Histone, y = Estimated_Library_Size, fill = Histone)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = Sample_ID), width = 0.15, size = 3, alpha = 0.8) +
  theme_bw(base_size = 14) + labs(y = "Estimated Library Size", x = "", color = "Sample ID")

figC <- ggplot(qc_summary_df, aes(x = Histone, y = Unique_Fragments, fill = Histone)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = Sample_ID), width = 0.15, size = 3, alpha = 0.8) +
  theme_bw(base_size = 14) + labs(y = "# of Unique Fragments", x = "", color = "Sample ID")

dup_panels <- ggarrange(figA, figB, figC, ncol = 3, common.legend = TRUE, legend = "bottom")
print(dup_panels)

dev.off()


# =========================================================================
# STEP 5: GENERATE mito_sex_plots.pdf
# =========================================================================
pdf("mito_sex_plots.pdf", width = 11, height = 6)

# -------------------------------------------------------------------------
# PANEL A: Global Mitochondrial Check
# -------------------------------------------------------------------------
p_mito_dot <- ggplot(qc_summary_df, aes(Sample, ChrM_Fraction)) +
  geom_point(size = 3, color = "darkred") +
  theme_bw() +
  labs(
    title = "Mitochondrial Read Fraction per Sample", 
    x = "Sample", 
    y = "chrM / total mapped reads"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_mito_dot)

# -------------------------------------------------------------------------
# PREPARATION: Pivot Sexual Chromosome Data to Long Format
# -------------------------------------------------------------------------
chrXY_long <- qc_summary_df %>%
  select(Sample, Sample_ID, Histone, ChrX_Fraction, ChrY_Fraction) %>%
  pivot_longer(cols = c(ChrX_Fraction, ChrY_Fraction), names_to = "chromosome", values_to = "fraction") %>%
  mutate(
    chromosome = recode(chromosome, ChrX_Fraction = "X", ChrY_Fraction = "Y"),
    percentage = fraction * 100
  )

# -------------------------------------------------------------------------
# PANEL B: Generate Sex Chromosome Percentage Pages per Histone Mark
# -------------------------------------------------------------------------
unique_histones <- unique(qc_summary_df$Histone)

for (current_histone in unique_histones) {
  
  # Filter long-format data for the specific histone subset
  target_long_df <- chrXY_long %>% filter(Histone == current_histone)
  
  p_sex_histone <- ggplot(target_long_df, aes(x = Sample_ID, y = percentage, fill = chromosome)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    theme_bw() +
    labs(
      title = paste("Chromosome X and Y Read Percentages - Target:", current_histone),
      x = "Biological Sample ID", 
      y = "Percentage of total mapped reads (%)", 
      fill = "Chromosome"
    ) +
    scale_fill_manual(values = c("X" = "#F8766D", "Y" = "#00BFC4")) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Print the clean page for the current histone mark into the PDF stream
  print(p_sex_histone)
}

dev.off()

# =========================================================================
# STEP 6: GENERATE fragmentlengthplots.pdf
# =========================================================================
pdf("fragmentlengthplots.pdf", width = 11, height = 5.5)

# Extract and iterate iteratively over unique patient/donor Sample_IDs
unique_ids <- unique(frag_data$Sample_ID)

for (current_id in unique_ids) {
  
  # Filter dataset containing BOTH histone marks for the current biological individual
  sample_subset <- frag_data %>% filter(Sample_ID == current_id)
  
  # -----------------------------------------------------------------------
  # PANEL A: Unified Violin + Line Chart Page for the Current Patient
  # -----------------------------------------------------------------------
  
  # Left sub-plot: Side-by-side Histone mark comparison violins
  p_violin <- ggplot(sample_subset, aes(x = Histone, y = fragLen, weight = Weight, fill = Histone)) +
    geom_violin(alpha = 0.7, scale = "width", color = "black") +
    scale_y_continuous(breaks = seq(0, 800, 100)) +
    coord_cartesian(ylim = c(0, 600)) +
    scale_fill_manual(values = c("H3K4me3" = "#4682B4", "H3K27Ac" = "#E17055")) +
    theme_bw(base_size = 12) + 
    labs(y = "Fragment Length (bp)", x = "") +
    theme(legend.position = "none")
  
  # Right sub-plot: Overlay nucleosomal periodicity profiles
  p_line <- ggplot(sample_subset, aes(x = fragLen, y = fragCount, color = Histone)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("H3K4me3" = "#4682B4", "H3K27Ac" = "#E17055")) +
    theme_bw(base_size = 12) + 
    coord_cartesian(xlim = c(0, 500)) +
    labs(x = "Fragment Length (bp)", y = "Count", color = "Histone Mark")
  
  # Combine left and right panels with a common legend at the bottom
  page_combo <- ggarrange(
    p_violin, p_line, 
    ncol = 2, 
    widths = c(1, 1.5), 
    common.legend = TRUE, 
    legend = "bottom"
  )
  
  page_combo <- annotate_figure(
    page_combo, 
    top = text_grob(paste("QC Fragment Size Summary - Biological ID:", current_id), face = "bold", size = 16)
  )
  
  # Print the consolidated Page 1 for this donor
  print(page_combo)
  
  # -----------------------------------------------------------------------
  # PANEL B: Linear vs Log Scale Comparison (Both Histones Combined)
  # -----------------------------------------------------------------------
  # Restrict evaluation range to 600bp for clear visualization profiles
  df_scale <- sample_subset %>% filter(fragLen <= 600)
  
  p_linear <- ggplot(df_scale, aes(x = fragLen, y = Weight, color = Histone)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("H3K4me3" = "#4682B4", "H3K27Ac" = "#E17055")) +
    theme_bw() +
    labs(
      title = "Linear scale profile", 
      x = "Fragment length (bp)", 
      y = "Normalized read density"
    ) +
    theme(legend.position = "none") # Remove to avoid redundant legends on layout
  
  p_log <- ggplot(df_scale, aes(x = fragLen, y = Weight, color = Histone)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = c(200, 400, 600), linetype = "dashed", color = "grey40", linewidth = 0.6) +
    scale_color_manual(values = c("H3K4me3" = "#4682B4", "H3K27Ac" = "#E17055")) +
    scale_y_log10() +
    scale_x_continuous(breaks = seq(0, 600, by = 200)) +
    theme_bw() +
    labs(
      title = "Log scale profile", 
      x = "Fragment length (bp)", 
      y = "Normalized read density (log)"
    ) +
    theme(legend.position = "none")
  
  # Use patchwork to bind the traces together and generate a common bottom legend
  page_scales <- (p_linear | p_log) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom")
  
  # Add the main header for the scales page
  page_scales <- page_scales + 
    plot_annotation(title = paste("Fragment Density Distribution Scales - Biological ID:", current_id),
                    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)))
  
  # Print Page 2 for this donor
  print(page_scales)
}

dev.off()


# =========================================================================
# STEP 7: GENERATE fragment_length_highlight_all_samples.pdf
# =========================================================================

# Define the function to highlight all histone marks belonging to a specific Sample_ID
plot_frag_highlight <- function(frag_dataframe, highlight_id) {
  
  # Create a logical flag based on Sample_ID instead of individual sub-samples
  plot_df <- frag_dataframe %>% 
    mutate(highlight = (Sample_ID == highlight_id)) %>% 
    filter(fragLen <= 600)
  
  ggplot(plot_df, aes(x = fragLen, y = Weight, group = Sample)) +
    # Background: Traces from all OTHER sample IDs in light grey
    geom_line(data = filter(plot_df, !highlight), color = "grey70", alpha = 0.20, linewidth = 0.4) +
    
    # Foreground: Highlight ALL histone marks for the CURRENT target Sample_ID
    geom_line(data = filter(plot_df, highlight), aes(color = Histone), linewidth = 1.3, alpha = 0.85) +
    
    # Structural nucleosomal thresholds
    geom_vline(xintercept = c(200, 400), linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    
    # Maintain strict consistent color mapping per Histone mark
    scale_color_manual(values = c("H3K4me3" = "#4682B4", "H3K27Ac" = "#E17055")) +
    scale_y_log10() +
    coord_cartesian(xlim = c(0, 600)) +
    labs(
      title = "Fragment Length Distribution Profile Overview",
      subtitle = paste0("Highlighted Biological Sample: ", highlight_id),
      x = "Fragment length (bp)", 
      y = "Relative frequency (log scale)", 
      color = "Histone Mark"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", hjust = 0),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
    )
}

# Open the PDF device to deploy the multi-page sample-ID report
pdf("fragment_length_highlight_all_samples.pdf", width = 8, height = 6)

# Loop iteratively over unique biological Sample_IDs instead of sub-samples
unique_ids <- unique(frag_data$Sample_ID)

for (id in unique_ids) {
  print(plot_frag_highlight(frag_data, id))
}

dev.off()