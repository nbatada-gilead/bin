#!/bin/bash

# Print usage if -h or --help is provided
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  echo "Usage: find_diff_expr_genes.sh seurat_rds barcodeA.txt barcodeB.txt output_file"
  echo "Example: ./find_diff_expr_genes.sh luoma2020_tregs_cluster_bioturing.rds barcodes_tregs_healthy.txt barcodes_tregs_colitis.txt output_diff_expr_genes.tsv"
  exit 0
fi

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
  echo "Error: Missing arguments."
  echo "Usage: find_diff_expr_genes.sh seurat_rds barcodeA.txt barcodeB.txt output_file"
  exit 1
fi

# Assign input arguments to variables
SEURAT_RDS=$1
BARCODE_A=$2
BARCODE_B=$3
OUTPUT_FILE=$4

# R script to run differential expression analysis
Rscript - <<EOF
# Load necessary libraries
library(Seurat)
packageVersion("Seurat")

# Load Seurat object
sobj <- readRDS("${SEURAT_RDS}")

# Read barcode files
barcodesA <- readLines("${BARCODE_A}")
barcodesB <- readLines("${BARCODE_B}")

# Subset Seurat object for each cell type
sobj_A <- subset(sobj, cells = barcodesA)
sobj_B <- subset(sobj, cells = barcodesB)

print(dim(sobj_A))
print(dim(sobj_B))

# Merge Seurat objects
sobj_combined <- merge(sobj_A, sobj_B)

# Add 'cell_type' metadata to distinguish between A and B cells
sobj_combined\$cell_type_temp <- ifelse(Cells(sobj_combined) %in% barcodesA, "A", "B")

# Check if 'data' layer exists and normalize if it doesn't
if (!"data" %in% Assays(sobj_combined)) {
  sobj_combined <- NormalizeData(sobj_combined)
}

#print(colnames(sobj_combined@meta.data))

# Set cell identity using 'cell_type' metadata column
Idents(sobj_combined) <- "cell_type_temp"

# Perform differential expression analysis
markers <- FindMarkers(sobj_combined, ident.1 = "A", ident.2 = "B", logfc.threshold = 0.25, test.use="wilcox", min.pct=0.25, min.diff.pct=0.25, min.cells.group=8)

# Write results to the specified output file
write.table(markers, file = "${OUTPUT_FILE}", sep = ",", quote = FALSE, row.names = TRUE)
EOF
