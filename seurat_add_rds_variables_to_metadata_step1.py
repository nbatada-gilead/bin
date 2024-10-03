#!/usr/bin/env python3

import argparse
import re
from collections import defaultdict
import json
import rpy2.robjects as ro
from rpy2.robjects import r, NA_Character

def clean_column_name(name):
    return name.lower().replace('.', '_').replace(' ', '_')

def load_dictionary(dict_file):
    with open(dict_file, 'r') as f:
        dictionary = json.load(f)
    return dictionary

def main():
    parser = argparse.ArgumentParser(description='Generate mapping_list.txt from Seurat RDS file with optional source columns and dictionaries')
    parser.add_argument('-r', '--rds', required=True, help='Path to Seurat RDS file')
    parser.add_argument('-o', '--output', default='mapping_list.txt', help='Output mapping list file')

    # Optional arguments for source columns
    parser.add_argument('--diseasename_src', help='Source column for diseasename_rds')
    parser.add_argument('--group_src', help='Source column for group_rds')
    parser.add_argument('--celltype_src', help='Source column for celltype_rds')
    parser.add_argument('--tissuesite_src', help='Source column for tissuesite_rds')
    parser.add_argument('--tissue_src', help='Source column for tissue_rds')

    # Optional arguments for dictionaries
    parser.add_argument('--disease_dict', help='Path to disease dictionary JSON file')
    parser.add_argument('--tissue_dict', help='Path to tissue dictionary JSON file')

    args = parser.parse_args()

    # Load dictionaries if provided
    disease_dict = {}
    if args.disease_dict:
        disease_dict = load_dictionary(args.disease_dict)
    else:
        # Default disease dictionary
        disease_dict = {
            "CRC": "coad",
            "Lung": "luad",
            "Breast Cancer": "brca",
            "Prostate Cancer": "prad",
            # Add more mappings as needed
        }

    tissue_dict = {}
    if args.tissue_dict:
        tissue_dict = load_dictionary(args.tissue_dict)
    else:
        # Default tissue dictionary
        tissue_dict = {
            "Blood": "pbmc",
            "Peripheral blood": "pbmc",
            "Umbilical cord blood": "pbmc",
            "Colon": "colon",
            "Colorectum": "colon",
            "Intestine": "colon",
            "rectum": "colon",
            "ileum": "colon",
            "small intestine": "colon",
            # Add more mappings as needed
        }

    # Load the Seurat object using rpy2
    try:
        readRDS = ro.r['readRDS']
        sobj = readRDS(args.rds)
        # Assign 'sobj' to R global environment
        ro.globalenv['sobj'] = sobj
    except Exception as e:
        print(f"Error loading RDS file: {e}")
        return

    # Extract meta.data column names and unique values using R code
    try:
        r('meta_data <- sobj@meta.data')
        colnames = r('colnames(meta_data)')
        meta_data_cols = list(colnames)

        # Get unique values for each column
        meta_data_uniques = {}
        for col in meta_data_cols:
            uniques = r(f'unique(as.character(meta_data${col}))')
            meta_data_uniques[col] = [str(val) for val in uniques if val != NA_Character]
    except Exception as e:
        print(f"Error extracting meta.data: {e}")
        return

    # Define the variables we need to find and their optional source columns
    variables = {
        'diseasename_rds': args.diseasename_src,
        'group_rds': args.group_src,
        'celltype_rds': args.celltype_src,
        'tissuesite_rds': args.tissuesite_src,
        'tissue_rds': args.tissue_src
    }

    # Predefined keywords for each variable
    keywords = {
        'tissue_rds': ['tissue', 'organ', 'sample_site'],
        'tissuesite_rds': ['site', 'location', 'anatomy', 'anatomical'],
        'celltype_rds': ['celltype', 'cell_type', 'cell', 'lineage', 'cell type'],
        'group_rds': ['group', 'condition', 'status', 'disease_state', 'batch'],
        'diseasename_rds': ['disease', 'diagnosis', 'pathology', 'tumor_type', 'tumor type', 'tumortype']
    }

    # Predefined value keywords for each variable
    value_keywords = {
        'tissue_rds': ['blood', 'lung', 'colon', 'tissue', 'liver', 'brain', 'pbmc'],
        'tissuesite_rds': ['adjacent', 'normal', 'tumor', 'metastasis'],
        'celltype_rds': ['t cell', 'b cell', 'macrophage', 'neuron', 'fibroblast', 'cd4', 'cd8', 'treg'],
        'group_rds': ['tumor', 'normal', 'control', 'disease', 'healthy'],
        'diseasename_rds': ['cancer', 'carcinoma', 'disease', 'diabetes', 'luad', 'coad']
    }

    # Initialize a dictionary to store possible matches
    possible_mappings = defaultdict(list)

    # Prepare the mappings dictionary
    mappings_dict = {}

    for var, provided_src_col in variables.items():
        dest_col = var
        src_col = provided_src_col  # Use provided source column if given

        if src_col and src_col in meta_data_cols:
            # Source column is provided and exists in meta.data
            src_values = meta_data_uniques[src_col]
        else:
            # Attempt to guess the source column
            # Analyze column names and values to find possible matches
            for col in meta_data_cols:
                col_clean = clean_column_name(col)
                unique_values = meta_data_uniques[col]
                unique_values_clean = [val.lower() for val in unique_values]

                # Check if column name contains any keywords
                for kw in keywords[var]:
                    if re.search(r'\b' + re.escape(kw.lower()) + r'\b', col_clean):
                        possible_mappings[var].append((col, 'name_match'))
                        break  # No need to check other keywords

                # Check if values contain any keywords
                for val_kw in value_keywords[var]:
                    if any(re.search(r'\b' + re.escape(val_kw.lower()) + r'\b', val.lower()) for val in unique_values):
                        possible_mappings[var].append((col, 'value_match'))
                        break  # No need to check other value keywords

            candidates = possible_mappings.get(var, [])
            if candidates:
                # Prioritize name matches over value matches
                for candidate in candidates:
                    if candidate[1] == 'name_match':
                        src_col = candidate[0]
                        break
                else:
                    src_col = candidates[0][0]  # Take the first value_match
                src_values = meta_data_uniques[src_col]
            else:
                print(f"No valid source column found for '{dest_col}'. Skipping.")
                continue

        # Remove NA values
        src_values = [val for val in src_values if val != 'NA']

        # Prepare mapping_from and mapping_to
        mapping_from = src_values.copy()
        mapping_to = [val.lower().replace(' ', '_') for val in src_values]

        # Apply dictionaries for diseasename_rds and tissue_rds
        if dest_col == 'diseasename_rds':
            mapping_to = [disease_dict.get(val, val.lower().replace(' ', '_')) for val in src_values]
            default = ''  # Use original value if not in dictionary
        elif dest_col == 'tissue_rds':
            mapping_to = [tissue_dict.get(val, val.lower().replace(' ', '_')) for val in src_values]
            default = ''  # Use original value if not in dictionary
        elif dest_col == 'group_rds':
            # Map 'normal' to 'control', and 'tumor' to 'tumor', others to 'other'
            group_mapping = {'normal': 'control', 'tumor': 'tumor'}
            mapping_to = [group_mapping.get(val.lower(), 'other') for val in src_values]
            default = 'other'
        elif dest_col in ['celltype_rds', 'tissuesite_rds']:
            # Copy values as is
            mapping_from = []
            mapping_to = []
            default = ''
        else:
            default = ''

        # Store the mapping information
        mappings_dict[var] = {
            'src': src_col,
            'from': mapping_from,
            'to': mapping_to,
            'default': default
        }

    # Generate meta_mappings.txt
    with open(args.output, 'w') as f:
        for var, mapping_info in mappings_dict.items():
            src_col = mapping_info['src']
            mapping_from = mapping_info['from']
            mapping_to = mapping_info['to']
            default = mapping_info['default']
            f.write(f'{var}.src={src_col}\n')
            f.write(f'{var}.from={",".join(mapping_from)}\n')
            f.write(f'{var}.to={",".join(mapping_to)}\n')
            f.write(f'{var}.default={default}\n\n')

    print(f"Mapping list saved to {args.output}. Please review and edit it as necessary.")

if __name__ == '__main__':
    main()
