#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <seurat_rds_file> <barcodes_fileA> <barcodes_fileB> <de_genes_file>"
    exit 1
fi

# Assign input arguments to variables
SEURAT_RDS=$1
BARCODES_FILE_A=$2
BARCODES_FILE_B=$3
DE_GENES_FILE=$4

# Run R script
Rscript -e "
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)

# Load the Seurat object
seurat_obj <- readRDS('${SEURAT_RDS}')

# Load barcodes for group A and group B
barcodes_A <- readLines('${BARCODES_FILE_A}')
barcodes_B <- readLines('${BARCODES_FILE_B}')

# Combine barcodes into a single data frame
cell_data <- data.frame(barcode = c(barcodes_A, barcodes_B),
                         group = rep(c('Group A', 'Group B'), times = c(length(barcodes_A), length(barcodes_B))))

# Load DE genes
de_genes <- read.table('${DE_GENES_FILE}', header = TRUE, sep = '\t', row.names=1)

# Extract top genes based on avg_log2FC
top_genes_up <- rownames(de_genes[de_genes$avg_log2FC > 0, ])[1:10]
top_genes_down <- rownames(de_genes[de_genes$avg_log2FC < 0, ])[1:10]
top_genes <- c(top_genes_up, top_genes_down)

# Check if top_genes is empty
if (length(top_genes) == 0) {
    stop('No genes found in DE genes data.')
}

# Extract data for top genes
data_to_plot <- as.data.frame(seurat_obj[['RNA']]@data[top_genes, ])  # Updated to use [['RNA']]
data_to_plot$cell_type <- ifelse(rownames(data_to_plot) %in% barcodes_A, 'Group A', 'Group B')

# Melt the data for plotting
data_melted <- melt(data_to_plot, id.vars = 'cell_type', variable.name = 'gene', value.name = 'expression', na.rm = TRUE)

# Check the melted data structure
print(str(data_melted))

# Create violin plot
ggplot(data_melted, aes(x = gene, y = expression, fill = cell_type)) +
    geom_violin(position = position_dodge(0.9), alpha = 0.7) +
    labs(title = 'Violin Plots of Top DE Genes',
         x = 'Gene',
         y = 'Expression') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme_minimal() +
    geom_jitter(position = position_jitterdodge(), size = 0.5, alpha = 0.5)

# Save the plot
ggsave('violin_plot_top_genes.png', width = 10, height = 6)
"

