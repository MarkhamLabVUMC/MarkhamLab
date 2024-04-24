# R script for producing volcano plots
# Note preliminary steps remain same to run differential expression analysis
# However pre-filtering of data is not done in order for volcano plots to work

################################################################################

library("DESeq2")
library(readr)  # for read_csv

# Load counts matrix
counts_file_path <- "D:/Analysis/counts_matrix.csv"
counts <- read_csv(counts_file_path, col_names = TRUE)
rownames(counts) <- counts$...1 # Replace '...1' with the actual name of the gene ID column if necessary
# Create copy of counts matrix with gene names
counts_with_gene_names <- counts
counts <- counts[,-1] # Drop the gene ID column so that only numeric count data remains

# Load metadata
metadata_file_path <- "D:/Analysis/10051_NM_metadata.csv"
metadata <- read_csv(metadata_file_path)
rownames(metadata) <- metadata$Sample
metadata <- metadata[,-1] # Remove the 'Sample' column after setting it as row names

# Rename rownames in metadata to match the column names in counts
rownames(metadata) <- colnames(counts)

# Check that the sample names in the counts matrix and metadata match
stopifnot(all(colnames(counts) == rownames(metadata)))

# which contains the conditions for each sample
metadata$Condition <- factor(metadata$Condition,
                             levels = c("Untreated", "TcdB", "TcsL", "H2O2", "Vehicle", "TcdA", "TcdB+A"))

# Now 'Untreated' is the first level and thus the reference level by default
# If you want to explicitly set it as the reference you can use relevel
metadata$Condition <- relevel(metadata$Condition, ref = "Untreated")


# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = DataFrame(metadata),
                              design = ~ Condition) # Replace 'Condition' with the actual condition variable name if different

# Differential Expression Analysis
# Run the DESeq pipeline
dds <- DESeq(dds)
res <- results(dds)
res

################################################################################

# Need to add gene names to res to label volcano plots

# Assigning row names from counts_with_gene_names to res
rownames(res) <- rownames(counts_with_gene_names)

# Now, create the gene mapping again
res$gene <- counts_with_gene_names$...1  # Adjust the column name as needed

# Check if any gene names are missing after assignment
if (any(is.na(res$gene))) {
  stop("Some gene names could not be found. Please check the data.")
} else {
  message("Gene names successfully merged!")
}

# Check the head of res to ensure correctness
head(res)

# Three Different types of volcano plots are produced

################################################################################

# 1. Unlabeled volcano plot using EnhancedVolcano

# Load EnhancedVolcano if not already loaded
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)

# Create the volcano plot using EnhancedVolcano
png("D:\\Analysis\\volcano_plot_unlabeled.png", width = 1200, height = 600)
EnhancedVolcano(res,
                lab = res$gene,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano plot',
                pCutoff = 0.05,  # General p-value cutoff for significance highlighting
                FCcutoff = 2,    # General fold change cutoff for significance highlighting
                pointSize = 2.0,
                labSize = 3.0,
                labCol = 'black',
                colAlpha = 0.8,
                drawConnectors = TRUE,
                selectLab = res$pvalue < 0.01 & abs(res$log2FoldChange) > 2.5)  # More stringent criteria for labeling
dev.off()

################################################################################

# 2. Labeled volcano plot using EnhancedVolcano

png("D:\\Analysis\\volcano_plot_labeled.png", width = 1200, height = 600)
EnhancedVolcano(res,
                lab = res$gene,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano plot',
                pCutoff = 0.05,  # General p-value cutoff for significance highlighting
                FCcutoff = 2,    # General fold change cutoff for significance highlighting
                pointSize = 2.0,
                labSize = 3.0,
                labCol = 'black',
                colAlpha = 0.8,
                drawConnectors = TRUE)
dev.off()

################################################################################

# 3. Unlabeled volcano plot using gg-Plot

# Install and load ggplot2
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
library(ggplot2)

png("D:\\Analysis\\volcano_plot_ggplot.png", width = 1200, height = 600)
# Prepare data for plotting
volcano_data <- as.data.frame(res)
volcano_data$gene <- rownames(volcano_data)

# Basic volcano plot
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue), color = padj < 0.05)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10 p-value", title = "Volcano Plot") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue")
dev.off()


# Install and load ggrepel
if (!requireNamespace("ggrepel", quietly = TRUE))
  install.packages("ggrepel")

# Load necessary libraries
library(ggplot2)
library(ggrepel)  # For better label placement

# Prepare the data for plotting
res$logP <- -log10(res$pvalue)  # Create a new column for the transformed p-values
res$significant <- res$pvalue < 0.05 & abs(res$log2FoldChange) > 1  # Define significance