#!/usr/bin/env Rscript

# =========================================================================
# 1. LOAD VISUALIZATION LIBRARIES
# =========================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggpubr)
})

# =========================================================================
# 2. CONSOLIDATE FRIP AND TSSE METRICS
# =========================================================================
message("[GLOBAL QC] Gathering and merging all per-sample CSV files...")

frip_files <- list.files(pattern = "\\\\.frip\\\\.csv$")
tsse_files <- list.files(pattern = "\\\\.tsse\\\\.csv$")

if (length(frip_files) == 0 || length(tsse_files) == 0) {
  stop("CRITICAL ERROR: No .frip.csv or .tsse.csv files found in the working directory.")
}

# Read and bind tables natively
df_frip <- frip_files %>% map_df(~read.csv(.x))
df_tsse <- tsse_files %>% map_df(~read.csv(.x))

# Master join: combine FRiP and TSSE data into a single long table
df_global <- full_join(df_frip, df_tsse, by = c("sample", "histone_mark", "label"))

# Export the absolute complete master table of the experiment
write.csv(df_global, "Experiment_QC_Summary.csv", row.names = FALSE)

# =========================================================================
# 3. GENERATE THE PLOTS
# =========================================================================
message("[GLOBAL QC] Generating comparative plots...")

# Theme template for clean aesthetics
pipeline_theme <- theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "#f2f2f2"),
        panel.grid.minor = element_blank())

# Plot 1: FRiP in Peaks
p1 <- ggplot(df_global, aes(x = histone_mark, y = frip_peaks, fill = label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), size = 1.5, alpha = 0.6) +
  labs(title = "Fraction of Reads in Peaks (FRiP)", x = "Histone Mark", y = "FRiP Score", fill = "Condition") +
  scale_fill_manual(values = c("withDups" = "#e41a1c", "noDups" = "#377eb8")) +
  pipeline_theme
ggsave("FRiP_peaks_comparison.png", plot = p1, width = 7, height = 5, dpi = 300)

# Plot 2: FRiP in TSS 2kb regions
p2 <- ggplot(df_global, aes(x = histone_mark, y = frip_tss_2kb, fill = label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), size = 1.5, alpha = 0.6) +
  labs(title = "FRiP around TSS (±2 kb)", x = "Histone Mark", y = "FRiP TSS Score", fill = "Condition") +
  scale_fill_manual(values = c("withDups" = "#e41a1c", "noDups" = "#377eb8")) +
  pipeline_theme
ggsave("FRiP_TSS_comparison.png", plot = p2, width = 7, height = 5, dpi = 300)

# Plot 3: TSSE Scores
p3 <- ggplot(df_global, aes(x = histone_mark, y = tsse_score, fill = label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), size = 1.5, alpha = 0.6) +
  labs(title = "TSS Enrichment Score (TSSE)", x = "Histone Mark", y = "TSSE Value", fill = "Condition") +
  scale_fill_manual(values = c("withDups" = "#e41a1c", "noDups" = "#377eb8")) +
  pipeline_theme
ggsave("TSSE_scores_comparison.png", plot = p3, width = 7, height = 5, dpi = 300)

# =========================================================================
# 4. PARSE MACS3 CUTOFF ANALYSIS LOGS
# =========================================================================
message("[GLOBAL QC] Parsing MACS3 cutoff analysis logs...")

cutoff_files_list <- list.files(pattern = "_cutoff_analysis\\\\.txt$")

if (length(cutoff_files_list) > 0) {
  df_cutoff_master <- cutoff_files_list %>% map_df(function(f) {
    clean_name <- gsub("_cutoff_analysis\\\\.txt", "", f)
    name_parts <- unlist(strsplit(clean_name, "\\\\."))
    s_id  <- name_parts[1]
    lbl   <- name_parts[2]
    
    lines <- readLines(f, warn = FALSE)
    data_lines <- lines[!grepl("^#", lines)]
    if(length(data_lines) > 1) {
      dt <- read.table(text = data_lines, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      dt$sample <- s_id
      dt$label  <- lbl
      return(dt)
    }
    return(NULL)
  })

  if (nrow(df_cutoff_master) > 0) {
    score_col <- colnames(df_cutoff_master)[1]
    
    p4 <- ggplot(df_cutoff_master, aes_string(x = score_col, y = "npeaks", color = "label", group = "interaction(sample, label)")) +
      geom_line(alpha = 0.5, size = 0.8) +
      labs(title = "MACS3 Cutoff Score Impact on Peak Calling", x = paste("Cutoff Score Threshold (", score_col, ")"), y = "Number of Called Peaks", color = "Condition") +
      scale_color_manual(values = c("withDups" = "#e41a1c", "noDups" = "#377eb8")) +
      facet_wrap(~sample, scales = "free_y") +
      theme_bw() +
      theme(strip.background = element_rect(fill = "#f2f2f2"))
      
    ggsave("MACS3_cutoff_peaks_impact.png", plot = p4, width = 10, height = 6, dpi = 300)
  }
} else {
  message("[GLOBAL QC] WARNING: No MACS3 cutoff analysis files found to plot.")
}

message("[GLOBAL QC] All comparison charts successfully exported!")