#!/bin/bash

# Check if RDS file is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <path_to_rds_file> [export_data (default: False)]"
  exit 1
fi

RDS_FILE="$1"
EXPORT_DATA=${2:-False}

# Check if file exists
if [ ! -f "$RDS_FILE" ]; then
  echo "Error: File '$RDS_FILE' not found!"
  exit 1
fi

# Run R script with the provided RDS file
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

# Get assay name
assay_name <- names(sobj@assays)[1]

# Print active identity
cat('Active Identity:\n')
print(head(sobj@active.ident, 1))
cat('--------------------\n')


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
