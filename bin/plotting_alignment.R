#!/usr/bin/env Rscript
library(ggplot2)
library(stringr)

log_files <- list.files(pattern = "\\\\.bowtie2\\\\.txt$")
data <- data.frame()
for (f in log_files) {
    lines <- readLines(f)
    sample_name <- str_remove(f, "\\\\.bowtie2\\\\.txt")
    
    rate_line <- lines[grep("overall alignment rate", lines)]
    rate <- as.numeric(str_extract(rate_line, "[0-9.]+"))
    data <- rbind(data, data.frame(Sample = sample_name, Rate = rate))
}
p <- ggplot(data, aes(x = Sample, y = Rate, fill = Sample)) +
     geom_bar(stat = "identity", width = 0.5, color = "black") +
     theme_minimal() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1),
           legend.position = "none") +
     labs(title = "Global Alignment Efficiency (Bowtie2)",
          x = "Sample",
          y = "Alignment (%)") +
     ylim(0, 100)

ggsave("alignment_efficiency.pdf", plot = p, width = 7, height = 5)
