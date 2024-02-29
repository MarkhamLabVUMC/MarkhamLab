#1. Mapping

#Need to install rtracklayer if you do not have
BiocManager::install("rtracklayer")

library(rtracklayer)

gtf_path <- "D:\kallisto_Homo_sapiens\Homo_sapiens.GRCh38.111.gtf.gz"

gtf <- rtracklayer::import(gtf_path)

tx2gene <- data.frame(
  TXNAME = gtf$transcript_id,
  GENEID = gtf$gene_name
)

# Remove potential duplicates
#tx2gene <- tx2gene[!duplicated(tx2gene$TXNAME), ]

# Set the proper column names for tximport
colnames(tx2gene) <- c("TXNAME", "GENENAME")

# Remove rows where GENENAME is NA
tx2gene <- tx2gene[!is.na(tx2gene$GENENAME), ]

head(tx2gene)

#2. Install tximport

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport")

# Load the tximport library

library(tximport)

# Set the base directory where all your "kallisto" directories are located
base_dir <- "D:\\kallisto_Homo_sapiens"

# Create the vector of file paths for your kallisto abundance.tsv files
# Assuming that the files are compressed as .tsv.gz
# If they are not compressed, remove the .gz extension from the file path
files <- file.path(base_dir, sprintf("sample_%04d", 1:30), "abundance.tsv")

# Set names for the files based on your sample identifiers
names(files) <- paste0("sample", sprintf("%04d", 1:30))

# Run tximport for kallisto TSV files with ignoreTxVersion set to TRUE
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# View the first few entries of the counts matrix
head(txi.kallisto.tsv$counts)

# Get the counts matrix from the txi.kallisto.tsv object
counts_matrix <- txi.kallisto.tsv$counts

counts_matrix <- round(counts_matrix, 0)
str(counts_matrix)

# Define the file names for the .matrix and .csv files
matrix_file <- "D:\\Analysis\\counts_matrix.matrix"
csv_file <- "D:\\Analysis\\counts_matrix.csv"

# Save as .matrix (plain text format)
write.table(counts_matrix, file = matrix_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# Save as .csv file
write.csv(counts_matrix, file = csv_file, quote = FALSE, row.names = TRUE)
