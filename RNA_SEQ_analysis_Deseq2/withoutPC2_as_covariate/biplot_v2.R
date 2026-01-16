#!/usr/bin/env Rscript

# Biplot visualization following strict best practices
# Shows sample distribution and gene expression patterns

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(patchwork)
  library(matrixStats)
  library(RColorBrewer)
})

#-------------------------------------------------------------------------------
# HELPER FUNCTION: Convert Ensembl to Gene Symbols
#-------------------------------------------------------------------------------
get_gene_symbol_mapping <- function() {
  # Collect all successful gene mappings from any model
  all_mappings <- data.frame()
  
  for (model_name in names(full_results)) {
    results <- full_results[[model_name]]$results
    if ("gene_symbol" %in% names(results)) {
      # Check if this model has actual gene symbols (not just ensembl IDs)
      test_symbols <- head(results$gene_symbol, 10)
      test_ids <- head(results$ensembl_id, 10)
      
      # If gene_symbol is different from ensembl_id, it's a valid mapping
      if (!all(test_symbols == test_ids)) {
        mapping <- results %>%
          dplyr::select(ensembl_id, gene_symbol) %>%
          mutate(ensembl_clean = gsub("\\.\\d+$", "", ensembl_id)) %>%
          filter(gene_symbol != ensembl_id) %>%
          distinct()
        
        all_mappings <- rbind(all_mappings, mapping)
      }
    }
  }
  
  # Remove duplicates, keeping first occurrence
  all_mappings <- all_mappings %>%
    group_by(ensembl_clean) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  return(all_mappings)
}

# Fix the convert_ensembl_to_genes function
convert_ensembl_to_genes <- function(ensembl_ids) {
  # Clean version numbers
  clean_ids <- gsub("\\.\\d+$", "", ensembl_ids)
  
  # Get the mapping from successful models
  mapping <- get_gene_symbol_mapping()
  
  if (nrow(mapping) > 0) {
    # Map the symbols
    symbols <- mapping$gene_symbol[match(clean_ids, mapping$ensembl_clean)]
    
    # For any unmapped, keep the cleaned ensembl ID
    symbols[is.na(symbols)] <- clean_ids[is.na(symbols)]
    
    message(sprintf("  Mapped %d/%d gene symbols", 
                    sum(!is.na(match(clean_ids, mapping$ensembl_clean))), 
                    length(clean_ids)))
    
    return(symbols)
  } else {
    message("  No gene symbol mapping available, using Ensembl IDs")
    return(clean_ids)
  }
}


better_treatment <- function() {
  # Get mapping from a successful model 
  mapping <- get_gene_symbol_mapping()
  
  if (nrow(mapping) == 0) {
    message("No successful gene mappings found")
    return()
  }
  
  # List of models that likely have ensembl IDs instead of symbols #man
  treatment_models <- c("treatment_main",
                        "treatment_sex_interaction",
                        "treatment_strain_interaction",
                        "three_way_interaction",
                        "treatment_effect_in_all_females",
                        "treatment_strain_interaction_in_all_females",
                        "treatment_effect_in_all_males",
                        "treatment_strain_interaction_in_all_males",
                        "treatment_effect_in_all_C57BL6",
                        "treatment_sex_interaction_in_all_C57BL6",
                        "treatment_effect_in_all_DBA2J",
                        "treatment_sex_interaction_in_all_DBA2J",
                        "treatment_effect_in_female_C57BL6",
                        "treatment_effect_in_female_DBA2J",
                        "treatment_effect_in_male_C57BL6",
                        "treatment_effect_in_male_DBA2J")
  for (model_name in treatment_models) {
    if (model_name %in% names(full_results)) {
      results <- full_results[[model_name]]$results
      
      # Check if this model needs fixing
      if (all(results$gene_symbol == results$ensembl_id)) {
        message(sprintf("Fixing gene symbols for: %s", model_name))
        
        # Clean and map
        results$ensembl_clean <- gsub("\\.\\d+$", "", results$ensembl_id)
        
        # Map symbols
        results$gene_symbol <- mapping$gene_symbol[match(results$ensembl_clean, 
                                                         mapping$ensembl_clean)]
        
        # Keep ensembl ID for unmapped genes
        unmapped <- is.na(results$gene_symbol)
        results$gene_symbol[unmapped] <- results$ensembl_clean[unmapped]
        
        # Update the results
        full_results[[model_name]]$results <<- results
        
        # Count successful mappings
        n_mapped <- sum(!unmapped)
        message(sprintf("  Mapped %d/%d genes", n_mapped, nrow(results)))
      }
    }
  }
}

better_treatment()


#-------------------------------------------------------------------------------
# MAIN BIPLOT FUNCTION
#-------------------------------------------------------------------------------

create_biplot_with_samples <- function(pca_results, full_results, model_name, 
                                       n_genes = 15, output_dir = "wald_test_results") {
  
  # Check if model exists in pca_results
  if (!(model_name %in% names(pca_results))) {
    message(sprintf("  Skipping %s: No PCA results available", model_name))
    return(NULL)
  }
  
  # Check if model has enough significant genes
  n_sig <- sum(full_results[[model_name]]$results$padj < 0.05, na.rm = TRUE)
  if (n_sig < 10 ) {
    message(sprintf("  Skipping %s: Only %d significant genes (threshold: 10)", 
                    model_name, n_sig))
    return(NULL)
  }
  
  message(sprintf("Creating biplot for: %s (%d significant genes)", model_name, n_sig))
  
  # Extract PCA components
  pca_obj <- pca_results[[model_name]]$pca_object
  pca_data <- pca_results[[model_name]]$pca_data
  var_explained <- pca_results[[model_name]]$var_explained
  
  # Get loadings for PC1 and PC2
  loadings <- as.data.frame(pca_obj$rotation[, 1:2])
  loadings$gene <- rownames(loadings)
  
  # Convert to gene symbols
  loadings$gene_symbol <- convert_ensembl_to_genes(loadings$gene)
  
  # Get significance information from results
  results_df <- full_results[[model_name]]$results
  clean_results <- results_df %>%
    mutate(ensembl_clean = gsub("\\.\\d+$", "", ensembl_id)) %>%
    dplyr::select(ensembl_clean, padj, log2FoldChange, significant)
  
  # Match loadings with significance
  loadings$ensembl_clean <- gsub("\\.\\d+$", "", loadings$gene)
  loadings <- loadings %>%
    left_join(clean_results, by = "ensembl_clean")
  
  # Calculate contribution to variance
  loadings$contribution <- sqrt(loadings$PC1^2 + loadings$PC2^2)
  
  # Prioritize significant genes, then by contribution
  loadings <- loadings %>%
    mutate(
      is_significant = !is.na(padj) & padj < 0.05,
      direction = case_when(
        log2FoldChange > 0 ~ "Up",
        log2FoldChange < 0 ~ "Down",
        TRUE ~ "NS"
      )
    ) %>%
    arrange(desc(is_significant), desc(contribution))
  
  # Select top genes
  top_genes <- head(loadings, n_genes)
  
  # Scale loadings for visualization (best practice scaling)
  pc1_range <- range(pca_data$PC1)
  pc2_range <- range(pca_data$PC2)
  
  scale_factor <- min(
    diff(pc1_range) / (2 * max(abs(top_genes$PC1))),
    diff(pc2_range) / (2 * max(abs(top_genes$PC2)))
  ) * 0.7
  
  top_genes$PC1_scaled <- top_genes$PC1 * scale_factor
  top_genes$PC2_scaled <- top_genes$PC2 * scale_factor
  
  
  # Define color schemes
  treatment_colors <- c("vehicle" = "#0000FF", "THC" = "#FF0000")
  sex_shapes <- c("female" = 16, "male" = 17)
  direction_colors <- c("Up" = "#FFA500", "Down" = "#00441b", "NS" = "#999999")
  
  # Determine which factors are present in this model
  n_treatments <- length(unique(pca_data$treatment))
  n_sexes <- length(unique(pca_data$sex))
  n_strains <- length(unique(pca_data$strain))
  
  # Create the biplot
  p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
    # Add confidence ellipses if enough samples
    {if (nrow(pca_data) >= 6 && n_treatments > 1)
      stat_ellipse(aes(color = treatment), level = 0.95, type = "t", 
                   linetype = "dashed", alpha = 0.3)} +
    # Sample points
    geom_point(aes(color = treatment, shape = sex), 
               size = 4, alpha = 0.7) +
    # Sample labels
    geom_text_repel(
      aes(label = sample_name, color = treatment),
      size = 2.5,
      max.overlaps = 20,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.size = 0.2,
      segment.alpha = 0.5,
      segment.color = "grey70",
      show.legend = FALSE
    ) +
    # Gene loading vectors
    geom_segment(
      data = top_genes,
      aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled),
      arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
      alpha = 0.6,
      color = "darkred",
      size = 0.7
    ) +
    # Gene labels with background
    geom_label_repel(
      data = top_genes,
      aes(x = PC1_scaled, y = PC2_scaled, 
          label = gene_symbol, fill = direction),
      size = 2.8,
      max.overlaps = 30,
      color = "white",
      fontface = "bold",
      alpha = 0.85,
      box.padding = 0.2,
      point.padding = 0.1,
      segment.size = 0.3,
      segment.color = "darkred",
      segment.alpha = 0.5,
      min.segment.length = 0
    ) +
    # Scales
    scale_color_manual(values = treatment_colors, name = "Treatment") +
    scale_shape_manual(values = sex_shapes, name = "Sex") +
    scale_fill_manual(values = direction_colors, name = "Gene Direction") +
    # Labels
    labs(
      title = paste("Biplot:", gsub("_", " ", model_name)),
      subtitle = sprintf(
        "%d significant genes | Top %d loadings | %d samples | PC1: %.1f%% var | PC2: %.1f%% var",
        n_sig, nrow(top_genes), nrow(pca_data), var_explained[1], var_explained[2]
      ),
      x = paste0("PC1 (", var_explained[1], "% variance)"),
      y = paste0("PC2 (", var_explained[2], "% variance)")
    ) +
    # Theme
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey50", size = 1)
    ) +
    # Add origin lines
    geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.3) +
    geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.3)
  
  return(p)
}

#-------------------------------------------------------------------------------
# MAIN EXECUTION
#-------------------------------------------------------------------------------

message("\n==== BIPLOT GENERATION ====")
message("Creating biplots for models with >5 significant genes")
message("==========================================\n")

# Create output directory
biplot_dir <- file.path(output_dir, "biplots")
dir.create(biplot_dir, showWarnings = FALSE, recursive = TRUE)

# Store all biplots
all_biplots <- list()

# Generate biplots for each model in pca_results
for (model_name in names(pca_results)) {
  biplot <- create_biplot_with_samples(
    pca_results = pca_results,
    full_results = full_results,
    model_name = model_name,
    n_genes = 15,
    output_dir = biplot_dir
  )
  
  if (!is.null(biplot)) {
    # Save individual biplot
    output_file <- file.path(biplot_dir, paste0("biplot_", model_name, "_with_samples.png"))
    ggsave(output_file, biplot, width = 12, height = 10, dpi = 300)
    
    output_file_pdf <- file.path(biplot_dir, paste0("biplot_", model_name, "_with_samples.pdf"))
    ggsave(output_file_pdf, biplot, width = 12, height = 10, dpi = 300)
    
    all_biplots[[model_name]] <- biplot
    message(sprintf("  Saved: %s", basename(output_file)))
  }
}

# Create combined biplot if multiple plots exist
if (length(all_biplots) > 1) {
  message("\nCreating combined biplot figure...")
  
  n_plots <- length(all_biplots)
  ncol <- min(3, ceiling(sqrt(n_plots)))
  nrow <- ceiling(n_plots / ncol)
  
  # Adjust for combined view
  all_biplots_adjusted <- lapply(all_biplots, function(p) {
    p + theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      plot.title = element_text(size = 11),
      plot.subtitle = element_text(size = 8)
    )
  })
  
  combined_biplot <- wrap_plots(all_biplots_adjusted, ncol = ncol, nrow = nrow) +
    plot_annotation(
      title = "Biplot Analysis: Sample Distribution and Gene Expression Patterns",
      subtitle = "Models with >5 significant genes | Arrows: gene loadings | Labels: samples & genes",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )
  
  combined_file <- file.path(output_dir, "biplots_combined_with_samples.pdf")
  ggsave(combined_file, combined_biplot,
         width = min(24, ncol * 8),
         height = min(20, nrow * 7),
         dpi = 300, limitsize = FALSE)
  
  message(sprintf("Combined biplot saved: %s", basename(combined_file)))
}

# Summary
message("\n==========================================")
message(sprintf("Biplots generated: %d models", length(all_biplots)))
message(sprintf("Output directory: %s", biplot_dir))
message("==========================================\n")