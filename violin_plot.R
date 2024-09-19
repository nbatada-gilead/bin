# Load necessary libraries
library(Seurat)
library(reshape2)
library(ggplot2)
library(dplyr)

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
seurat_file <- args[1]
barcodesA_file <- args[2]
barcodesB_file <- args[3]
de_genes_file <- args[4]

# Load Seurat object
seurat_obj <- readRDS(seurat_file)

# Load barcodes
barcodesA <- readLines(barcodesA_file)
barcodesB <- readLines(barcodesB_file)

# Load DE genes
de_genes <- read.table(de_genes_file, header = TRUE, sep = "\t")
top_genes <- c(as.character(de_genes[de_genes$log2FC > 0, "gene_name"][1:10]),
               as.character(de_genes[de_genes$log2FC < 0, "gene_name"][1:10]))

# Fetch data for top genes
data_to_plot <- FetchData(seurat_obj, vars = c(top_genes, "orig.ident"))

# Convert orig.ident to character
data_to_plot$orig.ident <- as.character(data_to_plot$orig.ident)

# Add metadata for grouping
data_to_plot$cell_type <- ifelse(rownames(data_to_plot) %in% barcodesA, "Group A", "Group B")

# Melt data for ggplot
data_melted <- melt(data_to_plot, id.vars = "cell_type", variable.name = "gene", value.name = "expression", na.rm = TRUE)

print(head(data_melted,10))
# Generate violin plot
plot <- ggplot(data_melted, aes(x = cell_type, y = expression, fill = cell_type)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +
  facet_wrap(~ gene, scales = "free_y") +  # Create separate plots for each gene
  theme_minimal() +
  labs(title = "Violin Plot of Top DE Genes",
       x = "Cell Type",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 0))

# Save the plot
ggsave("violin_plot.pdf", plot, width = 12, height = 8)

