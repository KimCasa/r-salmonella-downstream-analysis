# Load libraries
library(tidyverse)

# Read QC summary
qc <- read_csv("sample_qc_summary.csv", show_col_types = FALSE)

# Histogram of N50 (kb)
ggplot(qc, aes(x = N50_kb, fill = QC_PASS)) +
  geom_histogram(bins = 40, color = "white", alpha = 0.8) +
  scale_fill_manual(values = c("TRUE" = "#2b8cbe", "FALSE" = "#de2d26")) +
  labs(
    title = "Distribution of N50 per Sample",
    x = "N50 (kb)",
    y = "Number of Samples",
    fill = "QC Pass"
  ) +
  theme_minimal(base_size = 13)
