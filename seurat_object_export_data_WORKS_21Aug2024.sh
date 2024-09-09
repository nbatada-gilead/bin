#!/bin/bash

# Function to display help message
function show_help() {
  echo "Usage: $0 <path_to_rds_file> <format> <data_type> <filename_prefix>"
  echo "  <path_to_rds_file>  : Path to the Seurat object (RDS file)."
  echo "  <format>            : 'full' (default), '10x', or 'mm'."
  echo "                        - 'full' exports the meta data and full assay data."
  echo "                        - '10x' exports in 10x format (matrix, features, barcodes)."
  echo "                        - 'mm' exports in Matrix Market format."
  echo "  <data_type>         : 'counts' or 'scaled' (default: counts)."
  echo "  <filename_prefix>   : Prefix for all output files (optional)."
  echo "  -h, --help          : Display this help message."
  exit 0
}

# Display help if requested
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  show_help
fi

# Check if RDS file is provided
if [ -z "$1" ]; then
  echo "Error: <path_to_rds_file> is required."
  show_help
fi

RDS_FILE="$1"
FORMAT="${2:-full}"
DATA_TYPE="${3:-counts}"
PREFIX="${4:-}"

# Check if file exists
if [ ! -f "$RDS_FILE" ]; then
  echo "Error: File '$RDS_FILE' not found!"
  exit 1
fi

# Run R script with the provided RDS file
Rscript --vanilla -e "
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))

# Read the Seurat object
sobj <- readRDS('$RDS_FILE')
cat('Successfully read the Seurat object.\n')

# Export meta data
meta_file <- paste0('${PREFIX}${RDS_FILE%.rds}', '_meta.tsv')
write.table(sobj@meta.data, file=meta_file, sep='\t', quote=FALSE, col.names=NA, row.names=TRUE)
cat('Meta data written to', meta_file, '\n')

# Export data based on the format and data type
if ('$FORMAT' == 'full') {
  if ('$DATA_TYPE' == 'counts') {
    data_matrix <- GetAssayData(sobj, assay='RNA', layer='counts')
    data_file <- paste0('${PREFIX}${RDS_FILE%.rds}', '_counts.tsv')
  } else {
    data_matrix <- GetAssayData(sobj, assay='RNA', layer='scale.data')
    data_file <- paste0('${PREFIX}${RDS_FILE%.rds}', '_scaled.tsv')
  }
  write.table(as.matrix(data_matrix), file=data_file, sep='\t', quote=FALSE, col.names=NA, row.names=TRUE)
  cat('Data matrix written to', data_file, '\n')

} else if ('$FORMAT' == '10x') {
  if ('$DATA_TYPE' == 'counts') {
    data_matrix <- GetAssayData(sobj, assay='RNA', layer='counts')
  } else {
    data_matrix <- GetAssayData(sobj, assay='RNA', layer='scale.data')
  }
  
  # Export matrix
  matrix_file <- paste0('${PREFIX}${RDS_FILE%.rds}', '_matrix.mtx')
  cat('Attempting to write matrix file to', matrix_file, '\n')
  writeMM(data_matrix, file=matrix_file)
  cat('Matrix file written to', matrix_file, '\n')
  
  # Export features (genes)
  features_file <- paste0('${PREFIX}${RDS_FILE%.rds}', '_features.tsv')
  write.table(data.frame(gene_id = rownames(data_matrix),
                         gene_name = rownames(data_matrix),
                         feature_type = rep('Gene Expression', nrow(data_matrix))),
              file = features_file, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
  cat('Features file written to', features_file, '\n')

  # Export barcodes (cells)
  barcodes_file <- paste0('${PREFIX}${RDS_FILE%.rds}', '_barcodes.tsv')
  write.table(colnames(data_matrix), file=barcodes_file, sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
  cat('Barcodes file written to', barcodes_file, '\n')

} else if ('$FORMAT' == 'mm') {
  if ('$DATA_TYPE' == 'counts') {
    data_matrix <- GetAssayData(sobj, assay='RNA', layer='counts')
  } else {
    data_matrix <- GetAssayData(sobj, assay='RNA', layer='scale.data')
  }
  mm_file <- paste0('${PREFIX}${RDS_FILE%.rds}', '_matrix.mtx')
  cat('Attempting to write Matrix Market file to', mm_file, '\n')
  writeMM(data_matrix, file=mm_file)
  cat('Matrix Market file written to', mm_file, '\n')
}
"
