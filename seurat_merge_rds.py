#!/usr/bin/env Rscript

# Load required libraries
library(Seurat)
library(argparse)
library(Matrix)

# Create parser object
parser <- ArgumentParser(description = 'Merge multiple Seurat RDS files into a single Seurat object with a single counts layer, optionally merging data and scale.data layers, merging meta.data, and adding batch information.')

parser$add_argument('-f', '--files', nargs = '+', required = TRUE, help = 'List of RDS files to merge')
parser$add_argument('-o', '--output', required = FALSE, default = 'merged_counts_seurat.rds', help = 'Output RDS file name')
parser$add_argument('-d', '--include_data', action = 'store_true', help = 'Include data and scale.data layers in the merged Seurat object')

args <- parser$parse_args()

# Initialize lists to hold counts matrices, data matrices, scale.data matrices, meta.data, and sample names
counts_list <- list()
data_list <- list()
scale_data_list <- list()
meta_data_list <- list()
sample_names <- c()

# Load each Seurat object and extract counts data and meta.data
cat("Loading Seurat objects and extracting data...\n")
for (i in seq_along(args$files)) {
    filepath <- args$files[i]
    sample_name <- basename(filepath)
    cat("Reading", filepath, "...\n")
    seurat_obj <- readRDS(filepath)
    
    # Extract counts matrix
    counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
    counts_list[[sample_name]] <- counts

    # Extract data and scale.data if required
    if (args$include_data) {
        data <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
        data_list[[sample_name]] <- data
        
        scale_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "scale.data")
        scale_data_list[[sample_name]] <- scale_data
    }
    
    # Extract meta.data
    meta_data <- seurat_obj@meta.data
    # Add a column batch_rds with the name of the source input RDS file
    meta_data$batch_rds <- sample_name
    meta_data_list[[sample_name]] <- meta_data

    sample_names <- c(sample_names, sample_name)
    
    # Remove variables to free memory
    rm(seurat_obj)
    gc()
}

# Get the union of all genes across all counts matrices
cat("Computing the union of genes across all datasets...\n")
all_genes <- Reduce(union, lapply(counts_list, rownames))

# Align counts matrices
cat("Aligning counts matrices to have the same genes...\n")
counts_aligned_list <- list()
for (sample_name in names(counts_list)) {
    counts <- counts_list[[sample_name]]
    # Identify missing genes in the current counts matrix
    missing_genes <- setdiff(all_genes, rownames(counts))
    if (length(missing_genes) > 0) {
        # Create a zero matrix for missing genes
        zero_matrix <- Matrix(0, nrow = length(missing_genes), ncol = ncol(counts), sparse = TRUE)
        rownames(zero_matrix) <- missing_genes
        colnames(zero_matrix) <- colnames(counts)
        # Combine the existing counts matrix with the zero matrix
        counts <- rbind(counts, zero_matrix)
    }
    # Ensure the counts matrix has all genes in the same order
    counts <- counts[all_genes, , drop = FALSE]
    counts_aligned_list[[sample_name]] <- counts

    # Remove the original counts to save memory
    counts_list[[sample_name]] <- NULL
    gc()
}

# Merge counts matrices
cat("Merging counts matrices...\n")
merged_counts <- do.call(cbind, counts_aligned_list)
rm(counts_aligned_list)
gc()

# If include_data is TRUE, also merge data and scale.data
if (args$include_data) {
    # Merge data matrices
    cat("Merging data matrices...\n")
    data_aligned_list <- list()
    for (sample_name in names(data_list)) {
        data <- data_list[[sample_name]]
        missing_genes <- setdiff(all_genes, rownames(data))
        if (length(missing_genes) > 0) {
            zero_matrix <- Matrix(0, nrow = length(missing_genes), ncol = ncol(data), sparse = TRUE)
            rownames(zero_matrix) <- missing_genes
            colnames(zero_matrix) <- colnames(data)
            data <- rbind(data, zero_matrix)
        }
        data <- data[all_genes, , drop = FALSE]
        data_aligned_list[[sample_name]] <- data

        # Remove the original data to save memory
        data_list[[sample_name]] <- NULL
        gc()
    }
    merged_data <- do.call(cbind, data_aligned_list)
    rm(data_aligned_list)
    gc()

    # Merge scale.data matrices
    cat("Merging scale.data matrices...\n")
    scale_data_aligned_list <- list()
    for (sample_name in names(scale_data_list)) {
        scale_data <- scale_data_list[[sample_name]]
        missing_genes <- setdiff(all_genes, rownames(scale_data))
        if (length(missing_genes) > 0) {
            zero_matrix <- matrix(0, nrow = length(missing_genes), ncol = ncol(scale_data))
            rownames(zero_matrix) <- missing_genes
            colnames(zero_matrix) <- colnames(scale_data)
            scale_data <- rbind(scale_data, zero_matrix)
        }
        scale_data <- scale_data[all_genes, , drop = FALSE]
        scale_data_aligned_list[[sample_name]] <- scale_data

        # Remove the original scale.data to save memory
        scale_data_list[[sample_name]] <- NULL
        gc()
    }
    merged_scale_data <- do.call(cbind, scale_data_aligned_list)
    rm(scale_data_aligned_list)
    gc()
}

# Merge meta.data
cat("Merging meta.data...\n")
# Get all unique column names across all meta.data
all_meta_cols <- unique(unlist(lapply(meta_data_list, colnames)))
# Align meta.data frames
aligned_meta_data_list <- list()
for (sample_name in names(meta_data_list)) {
    meta_data <- meta_data_list[[sample_name]]
    missing_cols <- setdiff(all_meta_cols, colnames(meta_data))
    if (length(missing_cols) > 0) {
        for (col in missing_cols) {
            meta_data[[col]] <- NA
        }
    }
    meta_data <- meta_data[, all_meta_cols, drop = FALSE]
    aligned_meta_data_list[[sample_name]] <- meta_data

    # Remove the original meta_data to save memory
    meta_data_list[[sample_name]] <- NULL
    gc()
}
# Combine meta.data frames
merged_meta_data <- do.call(rbind, aligned_meta_data_list)
rm(aligned_meta_data_list)
gc()

# Ensure that rownames of merged_meta_data match colnames of merged_counts
rownames(merged_meta_data) <- colnames(merged_counts)

# Create new Seurat object with the merged counts matrix
cat("Creating new Seurat object with merged counts...\n")
sobj <- CreateSeuratObject(counts = merged_counts, meta.data = merged_meta_data)
rm(merged_counts)
rm(merged_meta_data)
gc()

# If include_data is TRUE, add merged data and scale.data to the Seurat object
if (args$include_data) {
    cat("Adding merged data and scale.data to the Seurat object...\n")
    sobj[["RNA"]] <- SetAssayData(sobj[["RNA"]], layer = "data", new.data = merged_data)
    rm(merged_data)
    gc()
    sobj[["RNA"]] <- SetAssayData(sobj[["RNA"]], layer = "scale.data", new.data = merged_scale_data)
    rm(merged_scale_data)
    gc()
}

# If include_data is FALSE, remove data and scale.data layers
if (!args$include_data) {
    cat("Removing data and scale.data layers...\n")
    DefaultAssay(sobj) <- "RNA"
    sobj[["RNA"]]@data <- new("dgCMatrix")
    sobj[["RNA"]]@scale.data <- matrix()
}

# Save the new Seurat object
cat("Saving merged Seurat object to", args$output, "...\n")
saveRDS(sobj, file = args$output)
cat("Merging complete.\n")
