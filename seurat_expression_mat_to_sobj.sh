#!/bin/bash

# Usage information when no arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <gene_expression_tsv> <metadata_csv> <output_seurat_rds>"
    echo "Optional: Provide a specific cluster number to export (e.g., Tregs)"
    exit 1
fi

# Get arguments
EXPRESSION_FILE="$1"
METADATA_FILE="$2"
OUTPUT_FILE="$3"

# Set a default output file if the user does not provide one
if [ -z "$OUTPUT_FILE" ]; then
    OUTPUT_FILE="output_seurat_object.rds"
fi

# Create the R script
R_SCRIPT=$(cat <<EOF
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)

# Read arguments
expression_data <- fread("$EXPRESSION_FILE", sep = "\t", header = TRUE, row.names = 1)
metadata <- read.csv("$METADATA_FILE", row.names = 1)

# Create Seurat object
seurat_object <- CreateSeuratObject(counts = as.matrix(expression_data))
seurat_object <- AddMetaData(seurat_object, metadata)

# Seurat processing pipeline
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object, npcs = 30)
seurat_object <- FindNeighbors(seurat_object, dims = 1:30)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
seurat_object <- RunUMAP(seurat_object, dims = 1:30)

# Display the columns available in meta.data
cat("Here are the available meta.data columns:\n")
print(colnames(seurat_object@meta.data))

# Ask the user whether to export specific clusters
cat("Do you want to export a specific cluster? (yes/no): ")
export_cluster <- tolower(readline())

if (export_cluster == "yes") {
    # Ask user for the column to use as 'Ident'
    cat("Enter the column from meta.data to use as Ident (e.g., 'seurat_clusters'): ")
    ident_column <- readline()
    
    # Set the column as the Ident
    Idents(seurat_object) <- seurat_object@meta.data[[ident_column]]
    
    # Ask which cluster(s) to export
    cat("Enter the cluster number(s) to export (comma-separated for multiple): ")
    cluster_nums <- readline()
    
    # Convert the user input into a vector of cluster numbers
    cluster_nums_vector <- strsplit(cluster_nums, ",")[[1]] %>% as.numeric()
    
    # Subset the Seurat object by the selected clusters
    seurat_subset <- subset(seurat_object, idents = cluster_nums_vector)
    
    # Save the subset as an RDS file
    saveRDS(seurat_subset, file = paste0("$OUTPUT_FILE", ".rds"))
    
    cat("Clusters", cluster_nums, "exported as", "$OUTPUT_FILE", ".rds\n")
    
} else {
    # Save the full Seurat object
    saveRDS(seurat_object, file = paste0("$OUTPUT_FILE", ".rds"))
    cat("Full Seurat object exported as", "$OUTPUT_FILE", ".rds\n")
}

EOF
)

# Run the R script
echo "$R_SCRIPT" | R --no-save
