args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 05_figures/01_plot_early_treatment.R config.R")
source(args[[1]])
suppressPackageStartupMessages(library(ggplot2))
summary_table <- read.csv("outputs/tables/02_early_treatment_by_exposure.csv")
plot <- ggplot(summary_table, aes(endpoint, proportion, fill = factor(A))) +
  geom_col(position = position_dodge(width = .8), width = .7) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = NULL, y = "Observed proportion", fill = "Baseline resistance (A)",
    title = "Early appropriate active treatment by baseline resistance") + theme_minimal(base_size = 12)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("outputs/figures/01_early_treatment_by_exposure.png", plot, width = 8, height = 5, dpi = 300)
