# This follows the DESeq2 tutorial: 
# https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input

# 1. Download DESeq2

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# 2. 
library("DESeq2")
library(readr)  # for read_csv

# Load counts matrix
counts_file_path <- "D:/MarkhamLab/Analysis/counts_matrix.csv"
counts <- read_csv(counts_file_path, col_names = TRUE)
rownames(counts) <- counts$...1 # Replace '...1' with the actual name of the gene ID column if necessary
# Create copy of counts matrix with gene names
counts_with_gene_names <- counts
counts <- counts[,-1] # Drop the gene ID column so that only numeric count data remains

# Load metadata
metadata_file_path <- "D:/MarkhamLab/Analysis/10051_NM_metadata.csv"
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


# Histogram of counts per condition

# Set row names of 'metadata' to the 'Sample' column values
rownames(metadata) <- paste0("sample", sprintf("%04d", as.integer(rownames(metadata))))

# Now retry checking for missing samples
#missing_samples <- setdiff(rownames(metadata), colnames(raw_counts))
#print(missing_samples)

library(ggplot2)

# Get raw counts from DESeqDataSet object
raw_counts <- counts(dds)

# Iterate over each unique condition
for (condition in unique(metadata$Condition)) {
  # Subset sample names for the current condition
  samples_in_condition <- rownames(metadata)[metadata$Condition == condition]
  
  # Subset counts for samples in the current condition
  counts_in_condition <- raw_counts[, samples_in_condition]
  
  # Calculate the sum of counts per gene across all samples in the condition
  gene_counts_sum <- rowSums(counts_in_condition)
  
  # Plot histogram of gene counts sum for the current condition
  p <- ggplot(data = data.frame(gene_counts_sum), aes(x = gene_counts_sum)) +
    geom_histogram(bins = 3, fill = "blue", color = "black") +
    ggtitle(paste("Histogram of gene counts for", condition)) +
    xlab("Sum of counts per gene") +
    ylab("Number of genes")
  
  print(p)
}

# Line Graph

library(ggplot2)
library(dplyr)

# Calculate average count per gene for each condition
average_counts_per_condition <- data.frame()
for (condition in unique(metadata$Condition)) {
  samples_in_condition <- rownames(metadata)[metadata$Condition == condition]
  counts_in_condition <- raw_counts[, samples_in_condition]
  
  # Calculate the average count per gene
  average_counts <- rowMeans(counts_in_condition)
  
  # Prepare a data frame for plotting
  condition_data <- data.frame(Gene = 1:nrow(raw_counts), 
                               AverageCount = average_counts, 
                               Condition = condition)
  average_counts_per_condition <- rbind(average_counts_per_condition, condition_data)
}
# Open a PNG device
png("D:\\MarkhamLab\\Analysis\\average_gene_counts_by_condition.png", width = 800, height = 600)

# Plot
ggplot(average_counts_per_condition, aes(x = Gene, y = AverageCount, color = Condition)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Average Gene Counts by Condition", x = "Gene (ordered by appearance)", y = "Average Count") +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Close the device
dev.off()

# 2. Pre-filtering
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# 3. Differential Expression Analysis
# Run the DESeq pipeline
dds <- DESeq(dds)
res <- results(dds)
res

write.csv(as.data.frame(res), 
          file="D:\\MarkhamLab\\Analysis\\results_table.csv")

# Install apeglm for log fold change below
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm")

library(apeglm)

# 4. LFC
# Log fold change shrinkage for visualization and ranking
resultsNames(dds)
# First, find the position of "Condition_TcdB.A_vs_Untreated" in the resultsNames(dds)
coefPosition <- which(resultsNames(dds) == "Condition_TcdB.A_vs_Untreated")

# Now use this position with lfcShrink
resLFC <- lfcShrink(dds, coef=coefPosition, type="apeglm")
resLFC

write.csv(as.data.frame(resLFC), 
          file="D:\\MarkhamLab\\Analysis\\results_table_LFC.csv")

# 5. Speed-up and parallelization thoughts
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocParallel")

library("BiocParallel")
# For windows, use SnowParam, not MulticoreParam: register(MulticoreParam(4))
register(SnowParam(4))

# 6. p-values & adjusted p-values
resOrdered <- res[order(res$pvalue),]

summary(res)

write.csv(as.data.frame(res), 
          file="D:\\MarkhamLab\\Analysis\\summarized_results.csv")

# Shows how many adjusted p-values were less than 0.1
sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)

write.csv(as.data.frame(res05), 
          file="D:\\MarkhamLab\\Analysis\\summarized_results_pvalue_cutoff.csv")

sum(res05$padj < 0.05, na.rm=TRUE)

# 7. Independent hypothesis weighting
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("IHW")

library("IHW")
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)

write.csv(as.data.frame(resIHW), 
          file="D:\\MarkhamLab\\Analysis\\summarized_results_IHW.csv")

sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult

# 8. Exploring and exporting results
# log2 fold changes attributable to a given variable 
#over the mean of normalized counts for all the samples 
# in the DESeqDataSet
plotMA(res, ylim=c(-2,2))

# MA-plot for the shrunken log2 fold changes
plotMA(resLFC, ylim=c(-2,2))

# Skip it, not useful, Ran for hours without output, stopped it
#idx <- identify(res$baseMean, res$log2FoldChange)
#rownames(res)[idx]

resultsNames(dds)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ashr")

# because we are interested in treated vs untreated, we set 'coef=2'
# For our data, means Condition_TcdB_vs_Untreated
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")

par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

plotCounts(dds, gene=which.min(res$padj), intgroup="Condition")

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="Condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=Condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

mcols(res)$description

write.csv(as.data.frame(resOrdered), 
          file="D:\\MarkhamLab\\Analysis\\condition_treated_results.csv")

resSig <- subset(resOrdered, padj < 0.1)
resSig

write.csv(as.data.frame(resSig), 
          file="D:\\MarkhamLab\\Analysis\\condition_treated_results_pvalue.csv")

# Need Type in metadata
#colData(dds)
#ddsMF <- dds
#levels(ddsMF$type)
#levels(ddsMF$type) <- sub("-.*", "", levels(ddsMF$type))
#levels(ddsMF$type)

# Data Transformations & Visualization
# Count Data Transformations
# Extracting Transformed Values
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

write.csv(as.data.frame(assay(vsd)), 
          file="D:\\MarkhamLab\\Analysis\\assay_vsd.csv")

# Effects of Transformations on Variance
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("vsn")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("hexbin")

# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))

# Heatmap of Count Matrix
# Need Type again
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pheatmap")

library("pheatmap")
# Displays top 20 genes with highest normalized counts
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Condition","Round")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# Attempting to get gene names and their counts

# Get the normalized counts
norm_counts <- rowMeans(counts(dds, normalized=TRUE))

# Order the genes by their normalized counts
ordered_genes <- order(norm_counts, decreasing=TRUE)

# Get the names of the top 20 genes from the counts matrix with gene names
top_genes <- rownames(counts_with_gene_names)[ordered_genes[1:20]]

# Get the normalized counts of the top 20 genes
top_counts <- norm_counts[ordered_genes[1:20]]

# Create a dataframe with the gene names and their normalized counts
df_top_genes <- data.frame(Gene = top_genes, NormalizedCounts = top_counts)

# Print the dataframe
print(df_top_genes)

write.csv(as.data.frame(df_top_genes), 
          file="D:\\MarkhamLab\\Analysis\\heatmap_findings.csv")

# Heatmap with gene names depicted
library("pheatmap")

# Assuming 'counts_with_gene_names' exists and has gene names in the first column
# Extract gene names using the 'select' indices
gene_names <- rownames(counts_with_gene_names)[select]

# Use 'assay(ntd)' to get the normalized counts for the top 20 genes
normalized_counts_top_genes <- assay(ntd)[select,]

# Assign gene names as row names to the normalized counts matrix
rownames(normalized_counts_top_genes) <- gene_names

# Prepare the annotation data frame for the columns
df <- as.data.frame(colData(dds)[, c("Condition", "Round")])

# Generate the heatmap with gene names displayed
pheatmap(normalized_counts_top_genes, cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)


# Convert the matrix to a DataFrame
normalized_counts_df <- as.data.frame(normalized_counts_top_genes)

# Add the gene names as the first column of the DataFrame
normalized_counts_df$Gene <- rownames(normalized_counts_df)

# Rearrange the DataFrame to have 'Gene' as the first column
normalized_counts_df <- normalized_counts_df[, c("Gene", setdiff(names(normalized_counts_df), "Gene"))]
normalized_counts_df <- normalized_counts_df[,-1]
# Now 'normalized_counts_df' is a DataFrame with the desired structure
print(normalized_counts_df)

write.csv(as.data.frame(normalized_counts_df), 
          file="D:\\MarkhamLab\\Analysis\\heatmap_findings_per_sample.csv")


# Heatmap with values on it
library(ggplot2)
library(tidyr)

# Assuming 'normalized_counts_top_genes' contains the top 20 genes with row names as gene names
# Convert the matrix to a long format for ggplot
data_long <- as.data.frame(normalized_counts_top_genes) %>%
  tibble::rownames_to_column("Gene") %>%
  tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "Value")

# Plot heatmap with ggplot2
ggplot(data_long, aes(x = Sample, y = Gene, fill = Value)) +
  geom_tile() + # Create the heatmap tiles
  geom_text(aes(label = sprintf("%.2f", Value)), size = 3) + # Annotate with values
  scale_fill_gradient(low = "blue", high = "red") + # Color gradient
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Improve X labels readability
  labs(fill = "Normalized\nCount") # Legend title

# Note: Adjust 'sprintf("%.2f", Value)' to change the format of the numbers shown

# PCA Plot
# Saved as PCA_Plot_Samples
plotPCA(vsd, intgroup=c("Condition", "Cells"))

# Open a PNG device
png("D:\\MarkhamLab\\Analysis\\PCA_Plot_Samples_Customized.png", width = 800, height = 600)

pcaData <- plotPCA(vsd, intgroup=c("Condition", "Cells"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Cells)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

# Close the device
dev.off()


dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

# What is ctrlGenes?
dds <- estimateSizeFactors(dds, controlGenes=ctrlGenes)
dds <- DESeq(dds)

# Combining 'Condition' and 'Cells' into a single factor for the analysis
dds$group <- factor(paste0(dds$Condition, dds$Cells))

# Updating the design formula to use the new 'group' factor
design(dds) <- ~ group

# Running DESeq2 with the updated design
dds <- DESeq(dds)

# Checking the results names, which will help you specify contrasts
resultsNames(dds)

# Check results for particular group
results(dds, name="group_TcdANL.902_vs_H2O2NL.902")

# Likelihood Ratio Test
dds <- DESeq(dds, test="LRT", reduced=~1)
res <- results(dds)

write.csv(as.data.frame(res), 
          file="D:\\MarkhamLab\\Analysis\\results_table_LRT.csv")

# Boxplot
# y-axis is log-transformed Cook's distances of the RNA-seq count data.
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

# Dispersion plot
plotDispEsts(dds)



par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-2.5,2.5)
resGA <- results(dds, lfcThreshold = 0.5, altHypothesis = "greaterAbs", test="Wald")
resLA <- results(dds, lfcThreshold=.5, altHypothesis="lessAbs", test="Wald")
resG <- results(dds, lfcThreshold=.5, altHypothesis="greater", test="Wald")
resL <- results(dds, lfcThreshold=.5, altHypothesis="less", test="Wald")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()



use <- res$baseMean > metadata(res)$filterThreshold
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)

################################################################################
library(DESeq2)
library(readr)  # Optional, for read_csv function which is faster than read.csv

# Load your counts data
counts_file_path <- "D:/MarkhamLab/Analysis/counts_matrix.csv"
counts <- read_csv(counts_file_path)

# Load your metadata
metadata_file_path <- "D:/MarkhamLab/Analysis/10051_NM_metadata.csv"
metadata <- read_csv(metadata_file_path)

# Ensure that 'metadata' DataFrame has the necessary columns and conditions set correctly
# Typically, 'metadata' should have a column 'group' or similar that indicates the condition or group for each sample

# If the first column of 'counts' is gene IDs, set it as row names
rownames(counts) <- counts$GeneID  # Change 'GeneID' to the actual column name in your CSV
counts <- counts[,-1]  # Remove the gene ID column after setting it as row names

# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = DataFrame(metadata),
  design = ~ group  # Ensure this matches a column name in your metadata
)

# Run the DESeq analysis
dds <- DESeq(dds)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Glimma")
legend("topright", fill=rev(colori), legend=rev(names(colori)))

## pre-load to avoid load messages in report
library(Glimma)
library(edgeR)
library(DESeq2)

glimmaMDS(dds)


# Run DESeq and subsequent analysis
dds <- DESeq(dds, quiet = FALSE)

# Use glimma to create an MA plot
glimmaMA(dds)


glimmaVolcano(dds)
################################################################################
# For volcano plot
# 2. 
library("DESeq2")
library(readr)  # for read_csv

# Load counts matrix
counts_file_path <- "D:/MarkhamLab/Analysis/counts_matrix.csv"
counts <- read_csv(counts_file_path, col_names = TRUE)
rownames(counts) <- counts$...1 # Replace '...1' with the actual name of the gene ID column if necessary
# Create copy of counts matrix with gene names
counts_with_gene_names <- counts
counts <- counts[,-1] # Drop the gene ID column so that only numeric count data remains

# Load metadata
metadata_file_path <- "D:/MarkhamLab/Analysis/10051_NM_metadata.csv"
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

# 3. Differential Expression Analysis
# Run the DESeq pipeline
dds <- DESeq(dds)
res <- results(dds)
res

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


# Load EnhancedVolcano if not already loaded
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)

# Create the volcano plot using EnhancedVolcano
png("D:\\MarkhamLab\\Analysis\\volcano_plot_unlabeled.png", width = 1200, height = 600)
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


# Labeled volcano plot
png("D:\\MarkhamLab\\Analysis\\volcano_plot_labeled.png", width = 1200, height = 600)
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

# gg-Plot Volcano Plot

# Install and load ggplot2
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
library(ggplot2)

png("D:\\MarkhamLab\\Analysis\\volcano_plot_ggplot.png", width = 1200, height = 600)
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