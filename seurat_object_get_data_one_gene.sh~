#!/bin/bash

# Help function
show_help() {
  echo "Usage: $0 <path_to_rds_file> <export_data (True/False)> <expression_type (counts/normalized)> <gene_names...>"
  echo "Options:"
  echo "  <path_to_rds_file> : Path to the Seurat object (.rds file)."
  echo "  <export_data>      : Whether to export the data (default: False)."
  echo "  <expression_type>  : Type of expression data to add to metadata (counts or normalized)."
  echo "  <gene_names>       : Space-separated list of gene names."
  exit 0
}

# Check for help argument
if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]] || [[ "$1" == "-H" ]] || [[ "$1" == "--HELP" ]]; then
  show_help
fi

# Check if RDS file is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <path_to_rds_file> <export_data (True/False)> <expression_type (counts/normalized)> <gene_names...>"
  exit 1
fi

RDS_FILE="$1"
EXPORT_DATA=${2:-False}
EXPRESSION_TYPE=${3:-counts}
GENES="${@:4}"

# Check if file exists
if [ ! -f "$RDS_FILE" ]; then
  echo "Error: File '$RDS_FILE' not found!"
  exit 1
fi

# Prepare gene list for R script
GENES_R=$(printf "'%s'," $(echo $GENES) | sed 's/,$//')

# Run R script with the provided RDS file and options
Rscript --vanilla -e "
suppressPackageStartupMessages(library(Seurat))

# Message before reading the file
cat('Reading Seurat object from:', '$RDS_FILE', '\n')
cat('--------------------\n')

# Read the Seurat object
sobj <- readRDS('$RDS_FILE')

# Message after reading the file
cat('Successfully read the Seurat object.\n')
print(dim(sobj))
cat('--------------------\n')

# Print Seurat version
cat('Seurat version:', as.character(packageVersion('Seurat')), '\n')
cat('--------------------\n')

# Print Seurat object summary
print(sobj)
cat('--------------------\n')

# Get assay name
assay_name <- names(sobj@assays)[1]

# Add gene expression data to metadata if genes are provided
genes_list <- c($GENES_R)
if (length(genes_list) > 0) {
  for (gene in genes_list) {
    if (gene %in% rownames(GetAssayData(sobj, assay = assay_name))) {
      if ('$EXPRESSION_TYPE' == 'counts') {
        expr_data <- GetAssayData(sobj, assay = assay_name, layer = 'counts')[gene, ]
      } else {
        expr_data <- GetAssayData(sobj, assay = assay_name, layer = 'data')[gene, ]
      }
      sobj@meta.data[[gene]] <- as.vector(expr_data)
      cat('Gene expression data added to metadata for:', gene, '\n')
    } else {
      cat('Warning: Gene', gene, 'not found in the dataset.\n')
    }
  }
  cat('--------------------\n')
}

# Create a separate file to store the barcode presence data
barcodes <- rownames(sobj@meta.data)
rds_barcodes <- colnames(GetAssayData(sobj, assay = assay_name))
barcode_presence <- ifelse(rds_barcodes %in% barcodes, 1, 0)
output_file <- paste0(substr('$RDS_FILE', 1, nchar('$RDS_FILE') - 4), '_barcode_presence.tsv')
write.table(data.frame(Barcode = rds_barcodes, Presence = barcode_presence), 
            file = output_file, sep = '\t', quote = FALSE, row.names = FALSE)
cat(sprintf('Barcode presence data written to %s\n', output_file))

# Export RNA counts and scaled data if requested
if ('$EXPORT_DATA' == 'True') {
   # Export meta.data as TSV
   meta_file <- paste0(substr('$RDS_FILE', 1, nchar('$RDS_FILE') - 4), '_meta.tsv')
   write.table(sobj@meta.data, file=meta_file, sep='\t', quote=FALSE, col.names=NA)
   cat(sprintf('Meta data written to %s\n', meta_file))

} else {
  cat('Meta data export skipped. Give argument True if you want to export meta data.\n')
}
"

