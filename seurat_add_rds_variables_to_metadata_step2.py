#!/usr/bin/env python3

import argparse
import os
import rpy2.robjects as ro
from rpy2.robjects import r

def parse_mapping_list(mapping_file):
    mappings = {}
    with open(mapping_file, 'r') as f:
        lines = f.readlines()

    for line in lines:
        line = line.strip()
        if not line:
            continue  # Skip empty lines
        key, value = line.split('=', 1)
        var_name, attr = key.split('.', 1)
        if var_name not in mappings:
            mappings[var_name] = {}
        mappings[var_name][attr] = value
    return mappings

def main():
    parser = argparse.ArgumentParser(description='Apply mappings to Seurat object and save updated object')
    parser.add_argument('-r', '--rds', required=True, help='Path to original Seurat RDS file')
    parser.add_argument('-m', '--mapping', required=True, help='Path to mapping list file')
    parser.add_argument('-o', '--output', required=True, help='Path to save updated Seurat RDS file')
    args = parser.parse_args()

    # Load the Seurat object using rpy2
    try:
        readRDS = ro.r['readRDS']
        sobj = readRDS(args.rds)
        # Assign 'sobj' to R global environment
        ro.globalenv['sobj'] = sobj
    except Exception as e:
        print(f"Error loading RDS file: {e}")
        return

    # Parse the mapping_list.txt
    mappings = parse_mapping_list(args.mapping)

    # For each mapping, apply it using R code
    for dest_col, mapping_info in mappings.items():
        src_col = mapping_info.get('src', '')
        if not src_col or src_col == '':
            print(f"No source column specified for '{dest_col}'. Skipping.")
            continue

        # Prepare mapping in R
        from_values = mapping_info.get('from', '')
        to_values = mapping_info.get('to', '')
        default_value = mapping_info.get('default', '')

        from_list = from_values.split(',') if from_values else []
        to_list = to_values.split(',') if to_values else []

        if len(from_list) != len(to_list):
            print(f"Mismatch in 'from' and 'to' lengths for '{dest_col}'. Skipping.")
            continue

        if not from_list and not to_list:
            # Copy values as is
            r_code = f'''
            # Ensure the column is character
            sobj@meta.data[["{src_col}"]] <- as.character(sobj@meta.data[["{src_col}"]])

            # Copy the source column to the destination column
            sobj@meta.data[["{dest_col}"]] <- sobj@meta.data[["{src_col}"]]

            # Convert to lower case and replace spaces with underscores
            sobj@meta.data[["{dest_col}"]] <- tolower(sobj@meta.data[["{dest_col}"]])
            sobj@meta.data[["{dest_col}"]] <- gsub(" ", "_", sobj@meta.data[["{dest_col}"]])
            '''
        else:
            # Create mapping dictionary in R
            mapping_r_list = ', '.join([f'"{k}" = "{v}"' for k, v in zip(from_list, to_list)])

            # Prepare R code to apply the mapping
            r_code = f'''
            # Ensure the column is character
            sobj@meta.data[["{src_col}"]] <- as.character(sobj@meta.data[["{src_col}"]])

            # Create mapping list in R
            mapping <- list({mapping_r_list})

            # Apply the mapping, handling NA values
            sobj@meta.data[["{dest_col}"]] <- sapply(sobj@meta.data[["{src_col}"]], function(x) {{
                if (is.na(x)) {{
                    # Handle NA values
                    return({f'"{default_value}"' if default_value else 'NA'})
                }} else if (x %in% names(mapping)) {{
                    return(mapping[[x]])
                }} else {{
                    # Use default value if specified, otherwise return original value
                    return({f'"{default_value}"' if default_value else 'x'})
                }}
            }})

            # Convert to lower case and replace spaces with underscores
            sobj@meta.data[["{dest_col}"]] <- tolower(sobj@meta.data[["{dest_col}"]])
            sobj@meta.data[["{dest_col}"]] <- gsub(" ", "_", sobj@meta.data[["{dest_col}"]])
            '''

        # Execute the R code
        try:
            ro.r(r_code)
            print(f"Applied mapping for '{dest_col}'.")
        except Exception as e:
            print(f"Error applying mapping for '{dest_col}': {e}")
            continue

    # Save the updated Seurat object
    try:
        saveRDS = ro.r['saveRDS']
        saveRDS(ro.globalenv['sobj'], args.output)
        print(f"Updated Seurat object saved to {args.output}")
    except Exception as e:
        print(f"Error saving updated RDS file: {e}")
        return

if __name__ == '__main__':
    main()

    
