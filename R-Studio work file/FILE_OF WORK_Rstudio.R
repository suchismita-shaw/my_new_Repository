library(tidyverse)
library(reshape2)

# Load data
expr_data <- read.csv("C:\\Users\\User\\Downloads\\merged_expression_data.csv")  

# Convert wide to long format
expr_long <- melt(expr_data, id.vars = "GeneSymbol_EnsemblGeneID",
                  variable.name = "Sample", value.name = "Expression")

# Extract metadata from sample names
expr_long <- expr_long %>%
  separate(Sample, into = c("Genotype", "Condition", "Sex"), sep = "_")

# Filter one gene
gene_line <- expr_long %>% 
  filter(GeneSymbol_EnsemblGeneID == "0610007C21Rik_ENSMUSG00000013622.8")

# Line plot (mean across genotype-condition-sex)
ggplot(gene_line, aes(x = interaction(Genotype, Condition, Sex), y = Expression, group = 1)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(size = 2) +
  labs(title = "Expression Trend of 0610007C21Rik",
       x = "Sample Group", y = "Expression Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate mean expression per gene
top_genes <- expr_long %>%
  group_by(GeneSymbol_EnsemblGeneID) %>%
  summarise(MeanExpr = mean(Expression)) %>%
  arrange(desc(MeanExpr)) %>%
  slice(1:5)

# Bar plot
ggplot(top_genes, aes(x = reorder(GeneSymbol_EnsemblGeneID, -MeanExpr), y = MeanExpr, fill = GeneSymbol_EnsemblGeneID)) +
  geom_bar(stat = "identity") +
  labs(title = "Top 5 Expressed Genes",
       x = "Gene", y = "Mean Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Create gene-wise expression matrix
expr_wide <- expr_data
rownames(expr_wide) <- expr_wide$GeneSymbol_EnsemblGeneID
expr_wide <- expr_wide[,-1]
expr_wide <- t(expr_wide) %>% as.data.frame()

# Scatter for two genes
ggplot(expr_wide, aes(x = `0610007C21Rik_ENSMUSG00000013622.8`, 
                      y = `0610007L01Rik_ENSMUSG00000053094.8`)) +
  geom_point(color = "darkred", size = 3) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  labs(title = "Scatter: Expression Correlation between 0610007C21Rik and 0610007L01Rik",
       x = "0610007C21Rik Expression", y = "0610007L01Rik Expression") +
  theme_minimal()

