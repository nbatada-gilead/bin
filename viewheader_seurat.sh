#!/bin/bash

# Usage: seurat_show_meta.sh sobj.rds

# Check if the file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 sobj.rds"
    exit 1
fi

# Store the RDS object path
SOBJ_PATH=$1

# R script to list meta.data columns and their unique values in a cleaner format
Rscript - <<EOF
suppressPackageStartupMessages(library("Seurat"))
library(Seurat)

# Load Seurat object
sobj <- readRDS("$SOBJ_PATH")

# Access the meta.data
meta_data <- sobj@meta.data

# Function to format and print unique values with frequencies
print_top_values <- function(col, values, top_n=10) {
    cat( col, "\n")
    cat("--------------------------------------------------------\n")
    
    # Format unique values and frequencies
    unique_values <- sort(table(values), decreasing = TRUE)
    top_values <- head(unique_values, top_n)
    
    # Print top values with padding for readability
    for (i in seq_along(top_values)) {
        cat(sprintf("%-30s %d\n", names(top_values)[i], top_values[i]))
    }
    
    if (length(unique_values) > top_n) {
        cat("... (showing top", top_n, "of", length(unique_values), "unique values)\n")
    }
    cat("\n")
}

# Iterate over each column and print top 10 unique values
for (col in colnames(meta_data)) {
    print_top_values(col, meta_data[[col]])
}
EOF

