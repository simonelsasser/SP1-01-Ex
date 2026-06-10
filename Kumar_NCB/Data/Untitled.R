# Install DESeq2 and ggplot2 if you haven't already
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages("ggplot2")

library(DESeq2)
library(ggplot2)

# 1. Load the count data and metadata
# Note: row.names = 1 sets the gene_id as the row names
counts_data <- read.table("rsem.merged.gene_counts_2021.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
meta_data <- read.csv("metaData.csv", header = TRUE, stringsAsFactors = FALSE)

# 2. Clean the counts data
# The second column 'transcript_id(s)' is character data and not needed for DESeq2, so we remove it
counts_data <- counts_data[, -1] 

# 3. Subset for the 'naive' condition
# We only want rows in meta_data where condition is 'naive'
meta_naive <- meta_data[meta_data$condition == "naive", ]
# Subset the columns in counts_data to match these naive samples
counts_naive <- counts_data[, meta_naive$Name]

# Ensure meta_naive row order exactly matches the column order of counts_naive
meta_naive <- meta_naive[match(colnames(counts_naive), meta_naive$Name), ]

# 4. Create the DESeq2 dataset
# The design formula specifies we are testing the effect of 'treatment'
dds <- DESeqDataSetFromMatrix(countData = counts_naive,
                              colData = meta_naive,
                              design = ~ treatment)

# 5. Run the DESeq2 pipeline
dds <- DESeq(dds)

# 6. Extract the results
# We specify the contrast to compare 'EZH2i' (treated) against 'none' (control)
res <- results(dds, contrast = c("treatment", "EZH2i", "none"))

# 7. Generate a Volcano Plot
# Convert results to a dataframe for plotting
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Remove rows with NA values for pvalue or padj to avoid plotting errors
res_df <- na.omit(res_df)

# Define significance thresholds
log2FC_threshold <- 1
padj_threshold <- 0.05

# Add a category column for coloring points in the plot
res_df$significance <- "Not Significant"
res_df$significance[res_df$log2FoldChange > log2FC_threshold & res_df$padj < padj_threshold] <- "Upregulated in EZH2i"
res_df$significance[res_df$log2FoldChange < -log2FC_threshold & res_df$padj < padj_threshold] <- "Downregulated in EZH2i"

# Plot using ggplot2
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Downregulated in EZH2i" = "blue", 
                                "Not Significant" = "grey", 
                                "Upregulated in EZH2i" = "red")) +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Differential Expression: Naive + EZH2i vs Naive",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Significance") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Display the plot
print(volcano_plot)

# Save the plot to a file
ggsave("volcano_plot_naive_EZH2i_vs_naive.png", plot = volcano_plot, width = 8, height = 6)

# Optional: Save the full results table
write.csv(as.data.frame(res), file = "deseq2_results_naive_EZH2i_vs_naive.csv")
