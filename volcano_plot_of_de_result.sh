#!/bin/bash

# Check for input file
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path_to_diff_expr_results>"
    exit 1
fi

# Assign input file
DIFF_EXPR_FILE=$1

# R script to create a volcano plot
Rscript -e "
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load the differential expression results
deg_results <- read.csv('$DIFF_EXPR_FILE', sep='\t', header=TRUE, row.names=1)

# Create a new column for -log10(p_value)
deg_results\$neg_log10_pval <- -log10(deg_results\$p_val)

# Identify top 10 upregulated and downregulated genes
top10_up <- deg_results %>%
  filter(avg_log2FC > 0) %>%
  arrange(desc(avg_log2FC)) %>%
  head(10)

top10_down <- deg_results %>%
  filter(avg_log2FC < 0) %>%
  arrange(avg_log2FC) %>%
  head(10)

# Create the volcano plot
volcano_plot <- ggplot(deg_results, aes(x = avg_log2FC, y = neg_log10_pval)) +
  geom_point(aes(color = (p_val < 0.05 & abs(avg_log2FC) > 1)), alpha=0.5, size=1) + # smaller and more transparent points
  geom_point(data = top10_up, aes(x = avg_log2FC, y = neg_log10_pval), color = 'blue', size = 2) +
  geom_point(data = top10_down, aes(x = avg_log2FC, y = neg_log10_pval), color = 'orange', size = 2) +
  geom_hline(yintercept = -log10(0.01), linetype = 'dashed', color = 'red') +
  geom_vline(xintercept = c(-1, 1), linetype = 'dashed', color = 'red') +
  scale_color_manual(values = c('grey', 'red')) +
  labs(title = 'Volcano Plot of Differential Expression',
       x = 'Log2 Fold Change',
       y = '-Log10 P-Value') +
  theme_minimal() +
  theme(legend.position = 'none')

# Add labels for top 10 genes with background
volcano_plot <- volcano_plot +
  geom_label(data = top10_up, aes(label = rownames(top10_up)), vjust = -0.5, color = 'blue', size = 3, fill = 'white', label.size = 0.5) +
  geom_label(data = top10_down, aes(label = rownames(top10_down)), vjust = 1.5, color = 'orange', size = 3, fill = 'white', label.size = 0.5)

# Save the plot
ggsave('volcano_plot.pdf', plot = volcano_plot, width = 8, height = 6)
"


