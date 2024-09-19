#!/usr/bin/env python3
import argparse
import rpy2.robjects as ro
from rpy2.robjects import r
import os

def update_meta_data(sobj_file, col_orig, col_new, value_mapping, rds_suffix):
    """
    Updates the Seurat object's meta.data with a new column based on the provided value mappings.
    Adds support for '_rest_' to handle all unspecified values.
    """
    # Load the Seurat object in R
    r(f'sobj <- readRDS("{sobj_file}")')

    # Parse the value mapping into a dictionary
    mapping_dict = dict(item.split(':') for item in value_mapping.split(','))

    # Check if _rest_ is specified in the mapping dictionary
    rest_value = mapping_dict.pop('_rest_', None)

    # Convert the mapping dictionary to an R list format
    mapping_r_list = ', '.join([f'"{k}" = "{v}"' for k, v in mapping_dict.items()])

    # Apply the mapping in R and handle potential indexing issues
    r(f'''
    # Create mapping list in R
    mapping <- list({mapping_r_list})

    # Print the mapping list for debugging
    cat("Mapping list:\\n")
    print(mapping)

    # Ensure the column is treated as a character vector
    sobj@meta.data[["{col_orig}"]] <- as.character(sobj@meta.data[["{col_orig}"]])
    
    # Apply the mapping with support for '_rest_' value
    sobj@meta.data[["{col_new}"]] <- sapply(sobj@meta.data[["{col_orig}"]], function(x) {{
        if (x %in% names(mapping)) {{
            return(mapping[[x]])
        }} else {{
            # Use '_rest_' value if specified, otherwise return original value
            return({f'"{rest_value}"' if rest_value else 'x'})
        }}
    }})

    # Save the updated Seurat object with the specified suffix
    saveRDS(sobj, "{sobj_file[:-4]}_{rds_suffix}.rds")

    # Print friendly feedback
    old_values <- unique(sobj@meta.data[["{col_orig}"]])
    new_values <- unique(sobj@meta.data[["{col_new}"]])
    cat("Added new column '{col_new}'.\\nOld values in '{col_orig}':", old_values, "\\nMapped to new values in '{col_new}':", new_values, "\\n")
    ''')

def main():
    # Initialize argument parser
    parser = argparse.ArgumentParser(
        description='Add a new column to Seurat meta.data based on value mappings with optional "_rest_" for unspecified values.',
        epilog='Example usage: python seurat_add_meta.py -rds Zhang2024_CRC_GSE236581_tregs.rds -orig Tissue -new Tissue_rds -map Tumor:colon_tumor,Normal:colon,Blood:pbmc,_rest_:other -s v2'
    )

    # Define named arguments
    parser.add_argument('-r', '--rds_file', type=str, required=True, help='Path to the Seurat RDS file.')
    parser.add_argument('-o', '--col_orig', type=str, required=True, help='The original column name in meta.data.')
    parser.add_argument('-n', '--col_new', type=str, required=True, help='The new column name to create in meta.data.')
    parser.add_argument('-m', '--value_mapping', type=str, required=True, help='Value mappings in the format old_value:new_value separated by commas. Use "_rest_" to handle other values.')
    parser.add_argument('-s', '--rds_suffix', type=str, required=True, help='Suffix to append to the new RDS file name (rds suffix optional).')

    # Parse arguments
    args = parser.parse_args()

    # Call the update_meta_data function
    update_meta_data(args.rds_file, args.col_orig, args.col_new, args.value_mapping, args.rds_suffix)

if __name__ == '__main__':
    main()

