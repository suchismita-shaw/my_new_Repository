library(DESeq2)

df <- read.csv("merged_expression_data.csv", row.names = 1)

# Separate gene names and Ensembl IDs
rownames(df) <- df$GeneSymbol_EnsemblGeneID
df <- df[, -1]  # Remove gene ID column after setting as rownames
sample_info <- data.frame(
  condition = c("Smoked", "Control", "Smoked", "Control", "Smoked", "Control", "Smoked"),
  genotype = c("alpha7G", "alpha7E260A:G", "alpha7E260A:G", "alpha7G", "alpha7G", "alpha7E260A:G", "alpha7E260A:G"),
  sex = c("Male", "Male", "Male", "Female", "Female", "Female", "Female")
)
rownames(sample_info) <- colnames(df)

dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = sample_info,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
head(res)

res_sig <- res[which(res$padj < 0.05), ]
res_sig <- na.omit(res_sig)
head(res_sig)

res_ordered <- res[order(res$padj), ]
write.csv(as.data.frame(res_ordered), file = "DEG_results.csv")

plotMA(res, ylim = c(-5, 5), main = "DESeq2 MA-plot")
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("condition", "genotype"))

library(pheatmap)

# Select top 30 DE genes
topgenes <- head(order(res$padj), 30)

# Extract normalized counts
norm_counts <- assay(vsd)[topgenes, ]

# Add gene names
rownames(norm_counts) <- rownames(res)[topgenes]

pheatmap(norm_counts, cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, annotation_col = sample_info)


library(EnhancedVolcano)

df <- read.csv("merged_expression_data.csv", row.names = 1)

rownames(res) <- rownames(df)
head(rownames(res))

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Differential Expression',
                pCutoff = 0.05,
                FCcutoff = 1
)

EnhancedVolcano(res,
                lab = rownames(res),                      # Gene labels
                x = 'log2FoldChange',                     # X-axis: fold change
                y = 'pvalue',                             # Y-axis: raw p-value (you can use 'padj' too)
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~ "p-value"),
                title = 'Differential Expression',
                subtitle = 'Smoked vs Control',
                pCutoff = 0.05,                           # p-value threshold
                FCcutoff = 1,                             # log2FC threshold
                pointSize = 3.0,
                labSize = 3.5,
                colAlpha = 0.7,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2') # Color scheme
)



