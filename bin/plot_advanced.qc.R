#!/usr/bin/env Rscript

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# ----------------------------------------------------------------
# PARSE COMMAND LINE ARGUMENTS
# ----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Find the indices where our flags are declared
frip_flag_idx   <- which(args == "--frip")
cutoff_flag_idx <- which(args == "--cutoffs")

if (length(frip_flag_idx) == 0 || length(cutoff_flag_idx) == 0) {
  stop("Missing parameters. Usage: plot_advanced.qc.R --frip [files...] --cutoffs [files...]")
}

# Extract file paths sandwiched between flags or until the end of args array
if (frip_flag_idx < cutoff_flag_idx) {
  frip_paths   <- args[(frip_flag_idx + 1):(cutoff_flag_idx - 1)]
  cutoff_paths <- args[(cutoff_flag_idx + 1):length(args)]
} else {
  cutoff_paths <- args[(cutoff_flag_idx + 1):(frip_flag_idx - 1)]
  frip_paths   <- args[(frip_flag_idx + 1):length(args)]
}

# ----------------------------------------------------------------
# 1. READ AND PREPARE FRiP DATA 
# ----------------------------------------------------------------
if (length(frip_paths) == 0) stop("No FRiP files were passed to the script.")

frip_list <- lapply(frip_paths, read.csv)
frip_all  <- do.call(rbind, frip_list)

# Coerce structural columns to factors and standardize names cleanly
frip_all <- frip_all %>%
  mutate(
    sample       = as.character(sample),
    histone_mark = as.character(histone_mark),
    label        = case_when(
      grepl("noDups", label)   ~ "Deduplicated",
      grepl("withDups", label) ~ "Duplicated",
      TRUE                     ~ "Unknown"
    )
  )

frip_all$label        <- as.factor(frip_all$label)
frip_all$histone_mark <- as.factor(frip_all$histone_mark)

# --- AUTOMATIC SAMPLE CLEANING ---
frip_all <- frip_all %>%
  mutate(sample_clean = mapply(function(s, m) {
    gsub(paste0("_", m), "", s)
  }, sample, histone_mark)) %>%
  mutate(sample_clean = as.factor(sample_clean))

# ----------------------------------------------------------------
# 2. READ AND PREPARE CUTOFF FILES 
# ----------------------------------------------------------------
if (length(cutoff_paths) == 0) stop("No cutoff summary files were passed to the script.")

cutoff_list <- lapply(cutoff_paths, function(path) {
  filename <- basename(path) # Ej: "BPES2_H3K27Ac.withDups_cutoff_analysis.txt"
  
  label_val <- "Unknown"
  if (grepl("withDups", filename)) label_val <- "Duplicated"
  if (grepl("noDups", filename))   label_val <- "Deduplicated"
  
  mark_val <- "Unknown"
  if (grepl("H3K27Ac", filename))   mark_val <- "H3K27Ac"
  if (grepl("H3K4me3", filename))   mark_val <- "H3K4me3"
  if (grepl("H3K27me3", filename))  mark_val <- "H3K27me3"
  if (grepl("H3K9me3", filename))   mark_val <- "H3K9me3"
  if (grepl("H3K36me3", filename))  mark_val <- "H3K36me3"
  if (grepl("H3K4me1", filename))   mark_val <- "H3K4me1"
  if (grepl("H3K4me2", filename))   mark_val <- "H3K4me2"
  
  sample_clean_val <- gsub(paste0("_", mark_val, ".*"), "", filename)
  
  df <- read.table(path, header = TRUE, comment.char = "#")
  df$sample_clean <- sample_clean_val
  df$label        <- label_val
  df$histone_mark <- mark_val
  
  return(df)
})

cutoff_all <- do.call(rbind, cutoff_list)

# ----------------------------------------------------------------
# COLOR PALETTE HELPERS (Shared across blocks)
# ----------------------------------------------------------------
metric_colors <- c("frip_peaks" = "#2c3e50", "frip_tss_2kb" = "#16a085")
metric_labels <- c("frip_peaks" = "FRiP Peaks", "frip_tss_2kb" = "FRiP TSS (2kb)")

# ----------------------------------------------------------------
# BLOCK 0: FULLY AUTOMATED SUBplots (Faceted & Sorted dynamically)
# ----------------------------------------------------------------
available_marks  <- unique(frip_all$histone_mark)
available_labels <- unique(frip_all$label)

for (mark in available_marks) {
  for (lbl in available_labels) {
    
    df_sub <- frip_all %>% 
      filter(histone_mark == mark & label == lbl) %>%
      pivot_longer(cols = c(frip_peaks, frip_tss_2kb), names_to = "metric", values_to = "score") %>%
      mutate(sample_order = reorder(sample_clean, -score * (metric == "frip_peaks")))
    
    if (nrow(df_sub) == 0) next
    
    p_sub <- ggplot(df_sub, aes(x = sample_order, y = score, fill = metric)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.9) +
      theme_minimal() +
      scale_fill_manual(values = metric_colors, labels = metric_labels) +
      labs(title = paste(mark, "FRiP Metrics (", lbl, ")"), 
           x = "Sample", y = "FRiP Score") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.position = "bottom", 
            legend.title = element_blank())
    
    file_name_sub <- paste0("FRiPs_", mark, "_", lbl, ".png")
    ggsave(file_name_sub, plot = p_sub, width = 6, height = 5, dpi = 150)
  }
}

# ----------------------------------------------------------------
# PLOT 1: Multi-metric Faceted Comparison per Histone Mark
# ----------------------------------------------------------------
for (mark in available_marks) {
  
  df_mark <- frip_all %>%
    filter(histone_mark == mark) %>%
    pivot_longer(cols = c(frip_peaks, frip_tss_2kb), 
                 names_to = "metric", 
                 values_to = "score") %>%
    mutate(sample_order = reorder(sample_clean, -score * (metric == "frip_peaks")))
  
  p_bars <- ggplot(df_mark, aes(x = sample_order, y = score, fill = metric)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.9) +
    facet_wrap(~ label, scales = "free_x") + 
    theme_minimal() +
    labs(
      title = paste("FRiP Metrics for", mark),
      x = "Sample", 
      y = "Proportion Score"
    ) +
    scale_fill_manual(values = metric_colors, labels = metric_labels) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  file_name_p1 <- paste0("FRiP_comparison_", mark, ".png")
  ggsave(file_name_p1, plot = p_bars, width = 9, height = 6, dpi = 150)
}

# ----------------------------------------------------------------
# PLOT 2: Boxplot comparing label (Dups vs noDups) split by Mark
# ----------------------------------------------------------------
p2 <- ggplot(frip_all, aes(x = label, y = frip_peaks, fill = label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2, aes(color = label)) +
  facet_wrap(~ histone_mark) +
  theme_bw() +
  labs(title = "FRiP Peaks Distribution: Duplicated vs Deduplicated",
       x = "Condition", y = "FRiP Peaks Score") +
  scale_fill_manual(values = c(
    "Duplicated"   = "#3498db", 
    "Deduplicated" = "#e74c3c"
  )) +
  scale_color_manual(values = c(
    "Duplicated"   = "#2980b9", 
    "Deduplicated" = "#c0392b"
  )) + 
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 11)
  )

ggsave("FRiP_boxplot.png", plot = p2, width = 8, height = 6, dpi = 150)

# ----------------------------------------------------------------
# PLOT 3: MACS3 Cutoff Impact Curve (Filtered by qscore >= 1)
# ----------------------------------------------------------------
cutoff_filtered <- cutoff_all %>%
  filter(qscore >= 1)

p3 <- ggplot(cutoff_filtered, aes(x = qscore, y = npeaks, group = interaction(sample_clean, label), color = label)) +
  geom_line(alpha = 0.7, size = 1) +
  facet_wrap(~ histone_mark, scales = "free_y") +
  theme_minimal() +
  labs(title = "MACS3 Cutoff Plot",
       subtitle = "Displaying high-confidence thresholds (qscore >= 1)",
       x = "Significance Threshold (qscore)", y = "Total Peaks (npeaks)") +
  scale_color_manual(values = c(
    "Duplicated"   = "#2980b9", 
    "Deduplicated" = "#c0392b"
  )) + 
  geom_vline(
    xintercept = 5,
    linetype = "dashed",
    color = "grey40"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

ggsave("MACS3_cutoff_peaks_impact_filtered.png", plot = p3, width = 9, height = 6, dpi = 150)

# ----------------------------------------------------------------
# 4. GLOBAL SYNTHESIS TABLE (.CSV)
# ----------------------------------------------------------------
table_dups <- frip_all %>% 
  filter(label == "Duplicated") %>% 
  select(sample = sample_clean, histone_mark, total_reads_withDups = total_reads, frip_tss_2kb_withDups = frip_tss_2kb, frip_peaks_withDups = frip_peaks)

table_nodups <- frip_all %>% 
  filter(label == "Deduplicated") %>% 
  select(sample = sample_clean, histone_mark, total_reads_noDups = total_reads, frip_tss_2kb_noDups = frip_tss_2kb, frip_peaks_noDups = frip_peaks)

final_table <- full_join(table_dups, table_nodups, by = c("sample", "histone_mark"))
write.csv(final_table, "Experiment_QC_Summary.csv", row.names = FALSE)