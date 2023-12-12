# loading the required libraries
library( "DESeq2" )
library(ggplot2)
library(tidyverse)
library(dplyr)

# Specifying the path to the gene count file
Gene_path <- "C:\\Users\\saini\\Downloads\\all_samples_gene_count.txt"

# reading the file into a dataframe
gene_df <- read.table(Gene_path, header = TRUE, sep = "\t")

# columns to be dropped
columns_to_drop <- c("Chr", "Start", "End", "Strand", "Length")

# Using square brackets to drop the specified columns
gene_df <- gene_df[, !colnames(gene_df) %in% columns_to_drop]

# trimming the column name to only have ascension number as sample name
for (i in 2:ncol(gene_df)) {  
  colnames(gene_df)[i] <- sub(".*\\.(ERR[0-9]+)_sorted\\.bam", "\\1", colnames(gene_df)[i])
}

# converting the first column to row names
gene_df <- gene_df %>% tibble::column_to_rownames(var = "Geneid")

View(gene_df)

# Store row names
rownames_deg_df <- rownames(gene_df)

# Convert the columns to numeric
gene_df <- sapply(gene_df, as.numeric)

# Restore row names
rownames(gene_df) <- rownames_deg_df

# sum of all the counts in each row
row_sums <- rowSums(gene_df)

# dropping gene for which the sum of counts is 0
gene_df <- gene_df[row_sums > 0, ]

# view
View(gene_df)

# Specifying the path to the miRNA count file
meta_path <- "C:\\Users\\saini\\Downloads\\E-MTAB-11917.txt"

# reading the file into a dataframe
metadata <- read.table(meta_path, header = TRUE, sep = "\t", fill = TRUE)

# Remove periods from column names
colnames(metadata) <- gsub("\\.", "", colnames(metadata))

# Drop columns with any NA values
metadata <- metadata[, colSums(is.na(metadata)) == 0]

# dropping every alternate row as the data is repeated twice for forward and reverse read
metadata <- metadata[c(TRUE, FALSE), ]

# listing noise words that need to be removed from the dataframe
words_to_remove <- c("Characteristics", "Comment", "Protocol", "Term")

# Remove specified words from column names
for (word in words_to_remove) {
  colnames(metadata) <- gsub(word, "", colnames(metadata), ignore.case = TRUE)
}

# replacing the Sample names with ERR ascension numbers as the fastq samples were named by ERR ascension number
metadata <- metadata %>%
  mutate(ENA_SAMPLE = case_when(
    ENA_SAMPLE == "ERS12375704" ~ "ERR9922199",
    ENA_SAMPLE == "ERS12375705" ~ "ERR9922200",
    ENA_SAMPLE == "ERS12375706" ~ "ERR9922201",
    ENA_SAMPLE == "ERS12375707" ~ "ERR9922202",
    ENA_SAMPLE == "ERS12375708" ~ "ERR9922203",
    ENA_SAMPLE == "ERS12375709" ~ "ERR9922204",
    ENA_SAMPLE == "ERS12375710" ~ "ERR9922205",
    ENA_SAMPLE == "ERS12375711" ~ "ERR9922206",
    TRUE ~ ENA_SAMPLE
  ))

# converting the first column to row names
rownames(metadata) <- metadata$ENA_SAMPLE

View(metadata)

# making sure the row names in deg_df matches to column names in metadata
all(colnames(gene_df) %in% rownames(metadata))

# are they in the same order?
all(colnames(gene_df) == rownames(metadata))

# Creating DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData = gene_df,
                              colData = metadata,
                              design = ~ disease)

# Running DESeq analysis
dds <- DESeq(dds)

# Check the levels of the "disease" factor
levels(dds$disease)

# Differential expression analysis
res<- results(dds, contrast = c("disease","normal", "primary sclerosing cholangitis"))

# Order the results by adjusted p-value
res <- res[order(res$padj),]

# Create a data frame with the ordered results
results_df <- data.frame(res)

results_df

# Print significant genes (padj < 0.05) with row names
significant_genes <- results_df[which(results_df$padj < 0.05), ]

View(significant_genes)

# getting DEGs with padj as 0.05 and log2foldchange 1 (positive and negative)
DEGS <- significant_genes[which(significant_genes$padj < 0.05 & significant_genes$log2FoldChange > 1 | significant_genes$log2FoldChange < -1), ]

View(DEGS)

# upregulated genes
upregulated_genes <- DEGS[which(DEGS$padj < 0.05 & DEGS$log2FoldChange > 1), ]

View(upregulated_genes)

# downregulated gene
downregulated_genes <- DEGS[which(DEGS$padj < 0.05 & DEGS$log2FoldChange < -1), ]

View(downregulated_genes)

# volcano plot
library(EnhancedVolcano)
EnhancedVolcano(res, lab = rownames(results_df), 
                x = 'log2FoldChange', 
                y = 'pvalue', 
                title = 'Volcano Plot',
                pCutoff = 0.001,  
                FCcutoff = 1.0,
                pointSize = 1.0,
                labSize = 3.0)

