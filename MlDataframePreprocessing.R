# Load necessary libraries
library(TCGAbiolinks)
library(dplyr)
library(randomForest)
library(ggplot2)

# ---------------- Step 1: Download TCGA-LUAD RNA-Seq Data ----------------
# Query for RNA-Seq data (FPKM format)
rna_query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Download and prepare RNA-Seq data
GDCdownload(rna_query)
rna_data <- GDCprepare(rna_query)

# ---------------- Step 2: Preprocess RNA-Seq Data ----------------
# Extract the "FPKM-UQ" assay matrix from the SummarizedExperiment object
rna_expression <- assays(rna_data)$fpkm_uq_unstrand  

# Log-transform FPKM data to stabilize variance
rna_expression_log <- log2(rna_expression + 1)

# Get row names of dataframe
row_names <- rownames(rna_expression_log)

# Extract Ensembl IDs from row names (remove decimal versions)
row_names_clean <- gsub("\\..*", "", row_names)

# Check if any of the Ensembl IDs are present in the row names
matching_rows <- row_names[which(row_names_clean %in% ensembl_ids)]

# Ensure rownames are gene IDs and colnames are sample IDs
rownames(rna_expression_log) <- rownames(rna_data)
colnames(rna_expression_log) <- colnames(rna_data)

ensembl_ids_to_keep <- c("ENSG00000075624.17", "ENSG00000090339.9", "ENSG00000100197.22", 
                         "ENSG00000100320.23", "ENSG00000103024.7", "ENSG00000109971.14", 
                         "ENSG00000110651.12", "ENSG00000111640.15", "ENSG00000120217.14", 
                         "ENSG00000130234.12", "ENSG00000132589.16", "ENSG00000134460.18", 
                         "ENSG00000135318.12", "ENSG00000135404.12", "ENSG00000137166.16", 
                         "ENSG00000137845.15", "ENSG00000139350.12", "ENSG00000146648.19", 
                         "ENSG00000150093.20", "ENSG00000159640.17", "ENSG00000164111.15", 
                         "ENSG00000165029.17", "ENSG00000167552.14", "ENSG00000169174.11", 
                         "ENSG00000171680.23")

filtered_data_ML <- rna_expression_log[rownames(rna_expression_log) %in% ensembl_ids_to_keep, ]
filtered_data_ML_geneExp <- geneExp[rownames(geneExp) %in% ensembl_ids_to_keep, ]


mostVarValuesLUSC <- apply(rna_expression_log, 2, max)
mostCountsValues <- apply(geneExp, 2, max)

max_values <- apply(rna_expression_log_LUSC, 2, max)  # Find max value per column
df_max <- data.frame(max_values) 

# Set row names same as original dataframe
rownames(df_max) <- rownames(filtered_data_ML)  

ggplot(mostVarValues, aes(x = rownames(mostVarValues), y = mostVarValues)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Top 5 RNA Expression Levels",
    x = "Sample ID",
    y = "Most Variable Values"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability

# ---------------- Step 3: Identify Gene with Most Drastic Changes ----------------
# Compute the variance of expression for each gene across all samples
gene_variance <- apply(rna_expression_log_LUSC, 1, var)

# Find the gene with the highest variance
most_variable_gene <- rownames(rna_expression_log_LUSC)[which.max(gene_variance)]
cat("Gene with the most drastic changes in expression:", most_variable_gene, "\n")


# Identify the most variable gene for each sample
max_diff_genes <- apply(rna_expression_log_LUSC, 2, function(sample) {
  which.max(abs(sample - mean(sample))) # Gene with the largest deviation in expression
})

# Get gene names for each sample
max_diff_gene_names <- rownames(rna_expression_log_LUSC)[max_diff_genes]

# Create a summary table of results
results <- data.frame(
  SampleID = colnames(rna_expression_log_LUSC),
  MostVariableGene = max_diff_gene_names
)

get_important_gene <- function(sample_id, expression_data) {
  # Create binary labels: 1 for the current sample, 0 for others
  sample_labels <- as.factor(ifelse(colnames(expression_data) == sample_id, "Sample", "Others"))
  
  # Fit random forest model
  rf_model <- randomForest(
    x = t(expression_data),  # Features: genes (transposed to samples as rows)
    y = sample_labels,       # Labels: binary for the target sample
    importance = TRUE,
    ntree = 500
  )
  
  # Extract the importance matrix
  importance_data <- as.data.frame(importance(rf_model))
  importance_data$Gene <- rownames(importance_data)
  
  # Sort by Mean Decrease in Accuracy and select the top gene
  top_gene <- importance_data[order(-importance_data$MeanDecreaseAccuracy), ][1, c("Gene", "MeanDecreaseAccuracy")]
  
  # Add sample ID for reference
  top_gene <- data.frame(SampleID = sample_id, Gene = top_gene$Gene, MeanDecreaseAccuracy = top_gene$MeanDecreaseAccuracy)
  
  return(top_gene)
}


# Run random forest for each sample to find the most important gene
important_genes <- lapply(colnames(rna_expression_log_LUSC), function(sample_id) {
  get_important_gene(sample_id, rna_expression_log_LUSC)
})

# Combine results into a single data frame
important_genes_df <- do.call(rbind, important_genes)



sample_id <- colnames(rna_expression_log)[4]  # Choose the first sample as a representative
sample_labels <- as.factor(ifelse(colnames(rna_expression_log) == sample_id, "Sample", "Others"))

# ---------------- Step 4: Identify Mean Decrease in Accuracy for Genes in Results ----------------

# Filter the expression data to include only genes in the results
filtered_genes <- unique(results$MostVariableGene)
rna_expression_filtered <- rna_expression_log_LUSC[filtered_genes, ]

# Create binary labels for the first sample
sample_id <- results$SampleID[1]
sample_labels <- as.factor(ifelse(colnames(rna_expression_filtered) == sample_id, "Sample", "Others"))

# Fit a random forest model for the first sample using the filtered genes
rf_model_filtered <- randomForest(
  x = t(rna_expression_filtered),  # Features: filtered genes (transposed to samples as rows)
  y = sample_labels,               # Labels: binary for the target sample
  importance = TRUE,
  ntree = 500
)

# Extract the importance matrix for the filtered genes
importance_filtered <- as.data.frame(importance(rf_model_filtered))
importance_filtered$Gene <- rownames(importance_filtered)

# Merge the results to include Mean Decrease in Accuracy for each gene
results_with_mda <- results %>%
  left_join(importance_filtered, by = c("MostVariableGene" = "Gene")) %>%
  select(SampleID, MostVariableGene, MeanDecreaseAccuracy)

# ---------------- Step 5: Save and Display Results ----------------
# Save the results to a file
write.csv(results_with_mda, "results_with_mda.csv", row.names = FALSE)

# Display the results
print(results_with_mda)


results_with_mda <- results %>%
  rowwise() %>%
  mutate(
    MeanDecreaseAccuracy = {
      # Create binary labels for the current sample
      sample_labels <- as.factor(ifelse(colnames(rna_expression_filtered) == SampleID, "Sample", "Others"))
      
      # Fit a random forest model for the current sample
      rf_model <- randomForest(
        x = t(rna_expression_filtered),  # Features: filtered genes (transposed to samples as rows)
        y = sample_labels,               # Labels: binary for the target sample
        importance = TRUE,
        ntree = 500
      )
      
      # Extract the Mean Decrease in Accuracy for the specific gene
      importance_data <- as.data.frame(importance(rf_model))
      importance_data$Gene <- rownames(importance_data)
      
      # Get the MDA for the gene of interest
      mda <- importance_data %>%
        filter(Gene == MostVariableGene) %>%
        pull(MeanDecreaseAccuracy)
      
      # Return the MDA value (NA if not found)
      ifelse(length(mda) > 0, mda, NA)
    }
  )


machineLearningDataframe <- inner_join(adjresults_odd, results_with_mda, 
                                       by = c("cases" = "SampleID"))




***************************
  
#Graph Design
  

  library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)  # Load tibble for rownames_to_column()

# Ensure the dataframe has Ensembl IDs as row names
rna_filtered <- rna_filtered %>%
  rownames_to_column(var = "EnsemblID")  # Convert row names to a column

# Convert to long format for ggplot
rna_long <- rna_filtered %>%
  tidyr::pivot_longer(cols = -EnsemblID, names_to = "SampleID", values_to = "ExpressionValue")

# Get the top 5 highest expression values across all samples
top5_values <- rna_long %>%
  arrange(desc(ExpressionValue)) %>%
  slice_head(n = 5)

# Plot the bar chart
ggplot(top5_values, aes(x = SampleID, y = ExpressionValue, fill = SampleID)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = EnsemblID), vjust = -0.5, size = 4) +  # Add Ensembl ID labels
  labs(
    title = "Top 5 Expression Values Across Samples",
    x = "Sample ID",
    y = "Expression Value"
  ) +
  theme_minimal()  
  
  ggplot(top5_values, aes(x = SampleID, y = ExpressionValue)) +
  geom_bar(stat = "identity", fill = "grey") +  # All bars in grey
  geom_text(aes(label = EnsemblID), vjust = -0.5, size = 4) +  # Add Ensembl ID labels
    labs(
       title = "Top 5 Expression Values Across Samples",
       x = "Sample ID",
       y = "Expression Value"
         ) 
       theme_minimal()
       
       ensembl_ids <- c("ENSG00000159640", "ENSG00000090339", "ENSMUSG00000000442", "ENSMUSG00000027150",
                        "ENSG00000110651", "ENSG00000169174", "ENSG00000135318", "ENSG00000146648",
                        "ENSMUSG00000024401", "ENSG00000139350", "ENSMUSG00000020053", "ENSG00000164111",
                        "ENSG00000137845", "ENSG00000075624", "ENSG00000120217", "ENSG00000137166",
                        "ENSG00000165029", "ENSG00000171680", "ENSG00000134460", "ENSG00000103024",
                        "ENSG00000167552", "ENSG00000150093", "ENSG00000130234", "ENSG00000100197",
                        "ENSG00000100320", "ENSMUSG00000029724", "ENSG00000135404", "ENSG00000132589",
                        "ENSG00000111640", "ENSG00000109971")