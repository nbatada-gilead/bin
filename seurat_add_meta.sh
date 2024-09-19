#!/bin/bash

# Usage: seurat_add_meta.sh sobj.rds col_orig=Tissue col_new=Tissue_rds value_mapping=first_value::new_first_value value_mapping=second_value::new_second_value

# Parse input arguments
sobj=$1
col_orig=$(echo "$2" | cut -d'=' -f2)
col_new=$(echo "$3" | cut -d'=' -f2)
shift 3
value_mapping=""

# Process value_mapping pairs
for arg in "$@"; do
    old_value=$(echo "$arg" | cut -d'::' -f1)
    new_value=$(echo "$arg" | cut -d'::' -f2)
    value_mapping="$value_mapping\"$old_value\" = \"$new_value\", "
done

# Remove the last comma and space
value_mapping=$(echo "$value_mapping" | sed 's/, $//')

# Create an R script to update the Seurat object
cat <<EOF > update_seurat.R
suppressMessages(library(Seurat))

# Load the Seurat object
sobj <- readRDS("$sobj")

# Create value mapping
value_mapping <- list($value_mapping)

# Update the new column based on mapping
sobj@meta.data[["$col_new"]] <- with(sobj@meta.data, ifelse($col_orig %in% names(value_mapping), value_mapping[[as.character($col_orig)]], as.character($col_orig)))

# Save the updated Seurat object
saveRDS(sobj, sub(".rds", "_nb_meta.rds", "$sobj"))

# Print friendly feedback
old_values <- unique(sobj@meta.data[[col_orig]])
new_values <- unique(sobj@meta.data[[col_new]])

cat("Added new column '$col_new'.\\nOld values in '$col_orig':", old_values, "\\nMapped to new values in '$col_new':", new_values, "\\n")
EOF

# Run the R script
Rscript update_seurat.R

# Clean up
rm update_seurat.R

