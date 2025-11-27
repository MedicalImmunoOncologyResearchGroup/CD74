
library(DESeq2)
library(dplyr)
library(ggplot2)


n_repeats<-5

for (i in 1:n_repeats) {
    # Read the data. This example is for spleen samples
    countfile <- paste0('/Users/oipulk/Documents/scRNASeq/data/Eleftheria_Maranou_Mar2024/analysis/figures/pseudobulk/spleen/counts/pseudobulk_counts_spleen_aggregatedDC_', i,'.csv')
    counts <- read.csv(countfile, row.names=1)

    meta_file <- paste0('/Users/oipulk/Documents/scRNASeq/data/Eleftheria_Maranou_Mar2024/analysis/figures/pseudobulk/spleen/counts/pseudobulk_metadata_spleen_aggregatedDC_', i,'.csv')
    metadata <- read.csv(meta_file, row.names=1)

                                        # Ensure counts and metadata are aligned
    common_samples <- intersect(colnames(counts), rownames(metadata))
    counts <- counts[, common_samples]
    metadata <- metadata[common_samples, ]

                                        # Set up condition, genotype, and pathogenicity as factors
    metadata$condition <- factor(metadata$condition)
    metadata$genotype <- factor(metadata$wt.ko, levels=c("wt", "ko"))
    metadata$pathogenicity <- factor(metadata$pathogenicity, levels=c("naive", "pathogenic"))

                                        # Create DESeqDataSet object
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = metadata,
                                  design = ~ genotype + pathogenicity + genotype:pathogenicity)

                                        # Run DESeq2
    dds <- DESeq(dds)

                                        # Get results for the interaction term
    res_interaction <- results(dds, name="genotypeko.pathogenicitypathogenic")

                                        # Convert results to data frame and add gene names
    res_df <- as.data.frame(res_interaction) %>%
        mutate(gene = rownames(.))

                                        # Order by adjusted p-value
    res_df_ordered <- res_df %>% arrange(padj)

                                        # Export results
    outfile <- paste0('/Users/oipulk/Documents/scRNASeq/data/Eleftheria_Maranou_Mar2024/analysis/figures/pseudobulk/spleen/deseq2_interaction_results/spleen_aggregatedDC_deseq2_interaction_results_', i,'.csv')
    write.csv(res_df_ordered, outfile, row.names = FALSE)
    }


### Extras:

## Repeats are dealt with later in Python

## Function to perform DESeq2 analysis
#perform_deseq2 <- function(counts, metadata) {
#  dds <- DESeqDataSetFromMatrix(countData = counts,
#                                colData = metadata,
#                                design = ~ genotype + pathogenicity + genotype:pathogenicity)
#  dds <- DESeq(dds)
#  res_interaction <- results(dds, name="genotypeko.pathogenicitypathogenic")
#  return(res_interaction)
#}

## If you need to handle repeats (e.g., for DCs with small cell numbers):
#handle_repeats <- function(counts, metadata, n_repeats=10) {
#  results_list <- list()
#  for (i in 1:n_repeats) {
#    # Randomly select one replicate for each unique sample
#    unique_samples <- unique(sub("_\\d+$", "", rownames(metadata)))
#    selected_samples <- sapply(unique_samples, function(s) {
#      sample(grep(paste0("^", s, "_"), rownames(metadata), value = TRUE), 1)
#    })
#    
#    subset_counts <- counts[, selected_samples]
#    subset_metadata <- metadata[selected_samples, ]
    
#    results_list[[i]] <- perform_deseq2(subset_counts, subset_metadata)
#  }
  
#  # Combine results (you may want to adjust this based on your needs)
#  combined_results <- do.call(cbind, lapply(results_list, function(x) x$log2FoldChange))
#  mean_log2FC <- rowMeans(combined_results)
#  combined_pvalues <- do.call(cbind, lapply(results_list, function(x) x$pvalue))
#  mean_pvalue <- apply(combined_pvalues, 1, function(x) exp(mean(log(x))))
  
#  final_results <- data.frame(
#    log2FoldChange = mean_log2FC,
#    pvalue = mean_pvalue,
#    padj = p.adjust(mean_pvalue, method = "BH")
#  )
#  return(final_results)
#}

## For B cells (no repeats needed)
#print("Performing DESeq2 analysis for B cells...")
#b_cell_results <- perform_deseq2(counts, metadata)
#write.csv(as.data.frame(b_cell_results), "b_cell_deseq2_results.csv")

## Uncomment the following lines if you need to analyze DCs with repeats
#print("Performing DESeq2 analysis for DCs with repeats...")
#dc_results <- handle_repeats(counts, metadata, n_repeats=10)
#write.csv(dc_results, "spleen_dc_deseq2_results_with_repeats.csv")

print("DESeq2 analysis complete. Results exported.")
