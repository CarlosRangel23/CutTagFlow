#!/usr/bin/env Rscript

# Load required libraries
library(stringr)
library(ggplot2)

# Initialize data frames for parsing
data <- data.frame(Sample = character(), Metric = character(), Reads = numeric(), Percentage = numeric(), stringsAsFactors = FALSE)
overall_data <- data.frame(Sample = character(), TotalReads = numeric(), OverallRate = numeric(), stringsAsFactors = FALSE)
log_files <- commandArgs(trailingOnly = TRUE)

for (f in log_files) {
  lines <- readLines(f)
  sample_name <- str_remove(basename(f), "\\.bowtie2\\.txt")
  
  # 1. Target specific lines in the Bowtie2 output
  line_total <- lines[grep("reads; of these", lines)]
  line_0     <- lines[grep("aligned concordantly 0 times", lines)]
  line_1     <- lines[grep("aligned concordantly exactly 1 time", lines)]
  line_mult  <- lines[grep("aligned concordantly >1 times", lines)]
  rate_line  <- lines[grep("overall alignment rate", lines)]
  
  # 2. Extract Absolute Reads (first integer in the line)
  total_reads <- as.numeric(str_extract(line_total, "[0-9]+"))
  reads_0     <- as.numeric(str_extract(line_0, "[0-9]+"))
  reads_1     <- as.numeric(str_extract(line_1, "[0-9]+"))
  reads_mult  <- as.numeric(str_extract(line_mult, "[0-9]+"))
  
  # 3. Extract Percentages (numbers immediately followed by '%')
  pct_0    <- as.numeric(str_extract(line_0, "[0-9.]+(?=%)"))
  pct_1    <- as.numeric(str_extract(line_1, "[0-9.]+(?=%)"))
  pct_mult <- as.numeric(str_extract(line_mult, "[0-9.]+(?=%)"))
  rate     <- as.numeric(str_extract(rate_line, "[0-9.]+(?=%)"))
  
  # 4. Append to data frames
  data <- rbind(data, data.frame(Sample = sample_name, Metric = "0 times", Reads = reads_0, Percentage = pct_0))
  data <- rbind(data, data.frame(Sample = sample_name, Metric = "Exactly 1 time", Reads = reads_1, Percentage = pct_1))
  data <- rbind(data, data.frame(Sample = sample_name, Metric = ">1 times", Reads = reads_mult, Percentage = pct_mult))
  
  overall_data <- rbind(overall_data, data.frame(Sample = sample_name, TotalReads = total_reads, OverallRate = rate))
}

# =========================================================================
# SECTION 1: CREATE A CONSOLIDATED DATA FRAME AND WRITE TO A SUMMARY TXT
# =========================================================================

# Build a comprehensive flat dataframe containing all information per sample
summary_df <- data.frame(Sample = character(), stringsAsFactors = FALSE)

for (s in unique(data$Sample)) {
  s_overall <- overall_data[overall_data$Sample == s, ]
  s_0 <- data[data$Sample == s & data$Metric == "0 times", ]
  s_1 <- data[data$Sample == s & data$Metric == "Exactly 1 time", ]
  s_mult <- data[data$Sample == s & data$Metric == ">1 times", ]
  
  row_data <- data.frame(
    Sample                  = s,
    Total_Input_Reads       = s_overall$TotalReads,
    Overall_Alignment_Rate_Pct = s_overall$OverallRate,
    Aligned_0_Times_Reads   = s_0$Reads,
    Aligned_0_Times_Pct     = s_0$Percentage,
    Aligned_1_Time_Reads    = s_1$Reads,
    Aligned_1_Time_Pct      = s_1$Percentage,
    Aligned_Mult_Times_Reads = s_mult$Reads,
    Aligned_Mult_Times_Pct   = s_mult$Percentage,
    stringsAsFactors        = FALSE
  )
  
  summary_df <- rbind(summary_df, row_data)
}

# Export the final consolidated dataframe to a tab-separated text file (.txt)
write.table(summary_df, 
            file = "alignment_summary_report.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# =========================================================================
# SECTION 2: GENERATE GGPLOT (STACKED BAR WITH OVERALL RATE ON TOP)
# =========================================================================
p <- ggplot(data, aes(x = Sample, y = Percentage, fill = Metric)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  
  geom_text(data = overall_data, 
            aes(x = Sample, y = 102, label = paste0(OverallRate, "%")), 
            inherit.aes = FALSE, 
            vjust = 0, fontface = "bold", size = 4.5, color = "black") +
  
  scale_fill_manual(values = c("0 times" = "#e41a1c", "Exactly 1 time" = "#4daf4a", ">1 times" = "#377eb8")) +
  
  labs(title = "Bowtie2 Concordant Alignment Distribution",
       subtitle = "Overall alignment rate displayed on top of bars",
       x = "Samples",
       y = "Percentage of Reads (%)",
       fill = "Concordant Alignment") +
  
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 110)

ggsave("alignment_summary_plot.png", plot = p, width = 8, height = 6, dpi = 300)