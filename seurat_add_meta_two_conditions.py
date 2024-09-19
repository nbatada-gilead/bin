#!/usr/bin/env python3

import argparse
import rpy2.robjects as ro
from rpy2.robjects import r
import os

def update_meta_data_using_operator(sobj_file, orig1, orig2, val1, val2, operator, new_col, new_val, rds_suffix):
    """
    Updates the Seurat object's meta.data with a new column based on a condition using multiple columns and an operator (AND/OR).
    If the new column already exists, only the rows that meet the condition are updated.
    """
    # Construct the query string based on the operator (AND or OR)
    if operator == "AND":
        query_str = f'{orig1} == "{val1}" & {orig2} == "{val2}"'
    elif operator == "OR":
        query_str = f'{orig1} == "{val1}" | {orig2} == "{val2}"'
    else:
        raise ValueError("Operator must be either 'AND' or 'OR'")

    # Load the Seurat object in R
    r(f'sobj <- readRDS("{sobj_file}")')

    # Apply the query in R
    r(f'''
    # Ensure the columns are treated as character vectors
    sobj@meta.data[["{orig1}"]] <- as.character(sobj@meta.data[["{orig1}"]])
    sobj@meta.data[["{orig2}"]] <- as.character(sobj@meta.data[["{orig2}"]])

    # Check if the new column already exists
    if (!("{new_col}" %in% colnames(sobj@meta.data))) {{
        # If it doesn't exist, create the new column and initialize it with NA or some default value
        sobj@meta.data[["{new_col}"]] <- NA
    }}

    # Select rows based on the query
    selected_rows <- with(sobj@meta.data, {query_str})

    # Update the new column for the selected rows
    sobj@meta.data[selected_rows, "{new_col}"] <- "{new_val}"

    # Save the updated Seurat object with the specified suffix
    saveRDS(sobj, "{sobj_file[:-4]}_{rds_suffix}.rds")

    # Print friendly feedback
    cat("Updated column '{new_col}' with value '{new_val}' for rows satisfying the condition: {query_str}")
    ''')

def main():
    # Initialize argument parser
    parser = argparse.ArgumentParser(
        description='Update Seurat meta.data based on multiple column conditions using a specified operator (AND/OR). If the new column exists, only update rows that meet the condition.',
        epilog='Example usage: python seurat_add_meta_using_operator.py -rds sobj.rds --orig1 "tissue" --orig2 "group" --val1 "value1" --val2 "value2" --operator "AND" --new_col new_colname --new_val "new_value" --suffix v1'
    )

    # Define named arguments
    parser.add_argument('-r', '--rds_file', type=str, required=True, help='Path to the Seurat RDS file.')
    parser.add_argument('--orig1', type=str, required=True, help='First original column to use in query.')
    parser.add_argument('--orig2', type=str, required=True, help='Second original column to use in query.')
    parser.add_argument('--val1', type=str, required=True, help='First value to filter on (in orig1 column).')
    parser.add_argument('--val2', type=str, required=True, help='Second value to filter on (in orig2 column).')
    parser.add_argument('--operator', type=str, required=True, choices=["AND", "OR"], help='Logical operator to combine conditions (AND/OR).')
    parser.add_argument('--new_col', type=str, required=True, help='The new column name to create in meta.data or update if it already exists.')
    parser.add_argument('--new_val', type=str, required=True, help='The new value to assign for selected rows.')
    parser.add_argument('--suffix', type=str, required=True, help='Suffix to append to the new RDS file name.')

    # Parse arguments
    args = parser.parse_args()

    # Call the function to update the meta.data
    update_meta_data_using_operator(args.rds_file, args.orig1, args.orig2, args.val1, args.val2, args.operator, args.new_col, args.new_val, args.suffix)

if __name__ == '__main__':
    main()

