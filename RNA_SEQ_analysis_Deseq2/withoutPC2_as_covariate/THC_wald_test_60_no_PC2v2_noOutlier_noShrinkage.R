# Run Wald test, with main, interaction and stratified analysis
# WITHOUT PC2 AS A COVARIATE
# USING RAW P-VALUES FOR ALL DOWNSTREAM ANALYSES

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm) 
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(ggrepel)
  library(patchwork)
  library(biomaRt)
  library(matrixStats)
  library(scales)
  library(ashr)
})

# ---------------------------------------------------------------
# 1. LOAD DATA 
# ---------------------------------------------------------------

message("=== Loading Data ===")

# Load counts
counts <- read.csv("counts_timepoint60no_outlier.csv", row.names = 1, check.names = FALSE)

metadata <- read.csv("metv2_timepoint60_no_PC2_no_outlier.csv", row.names = 1, check.names = FALSE)
all(row.names(metadata) %in% colnames(counts))

# Remove non-count columns
count_cols <- sapply(counts, is.numeric)
counts <- counts[, count_cols]

# Convert to integer matrix
counts <- round(as.matrix(counts))
mode(counts) <- "integer"

# Match samples between counts and metadata
common_samples <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, common_samples]
metadata <- metadata[common_samples, ]

message(sprintf("Samples: %d", ncol(counts)))
message(sprintf("Genes: %d", nrow(counts)))

# Set factor levels (reference levels: vehicle, female, C57BL_6)
metadata$treatment <- factor(metadata$Treatment, levels = c("vehicle", "THC"))
metadata$sex <- factor(metadata$Sex, levels = c("female", "male"))
metadata$strain <- factor(metadata$Strain, levels = c("C57BL_6", "DBA_2J"))

message("=== USING RAW P-VALUES (pvalue < 0.05) FOR ALL ANALYSES ===\n")

# ---------------------------------------------------------------
# 2. GENE NAME CONVERSION FUNCTION
# ---------------------------------------------------------------

convert_to_symbols <- function(ensembl_ids, batch_size = 3000, max_retries = 3) {
  message("  Converting Ensembl IDs to gene symbols...")
  
  # Clean version numbers
  clean_ids <- gsub("\\.\\d+$", "", ensembl_ids)
  unique_ids <- unique(clean_ids)
  
  if (length(unique_ids) == 0) {
    return(ensembl_ids)
  }
  
  # Try to load cached gene symbols if available
  cache_file <- "gene_symbols_cache.rds"
  if (file.exists(cache_file)) {
    message("  Loading gene symbols from cache...")
    cached_symbols <- readRDS(cache_file)
    symbols <- cached_symbols[clean_ids]
    symbols[is.na(symbols)] <- clean_ids[is.na(symbols)]
    return(symbols)
  }
  
  for (attempt in 1:max_retries) {
    result <- tryCatch({
      message(sprintf("  Attempt %d/%d to connect to BioMart...", attempt, max_retries))
      
      # Connect to BioMart with timeout
      mart <- useMart("ensembl", 
                      dataset = "mmusculus_gene_ensembl", 
                      host = "https://sep2025.archive.ensembl.org")
      
      # Initialize mapping
      symbols_map <- character(0)
      n_batches <- ceiling(length(unique_ids) / batch_size)
      
      for (i in seq_len(n_batches)) {
        start_idx <- (i - 1) * batch_size + 1
        end_idx <- min(i * batch_size, length(unique_ids))
        batch_ids <- unique_ids[start_idx:end_idx]
        
        message(sprintf("    Processing batch %d/%d (%d IDs)...", i, n_batches, length(batch_ids)))
        
        gene_map_batch <- getBM(
          attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'),
          filters = 'ensembl_gene_id',
          values = batch_ids,
          mart = mart
        )
        
        batch_symbols <- gene_map_batch$external_gene_name[match(batch_ids, gene_map_batch$ensembl_gene_id)]
        batch_symbols[is.na(batch_symbols) | batch_symbols == ""] <- batch_ids[is.na(batch_symbols) | batch_symbols == ""]
        names(batch_symbols) <- batch_ids
        symbols_map <- c(symbols_map, batch_symbols)
      }
      
      # Map back to original order
      symbols <- rep(NA_character_, length(clean_ids))
      names(symbols) <- clean_ids
      matching_indices <- match(names(symbols_map), names(symbols))
      symbols[matching_indices[!is.na(matching_indices)]] <- symbols_map[!is.na(matching_indices)]
      symbols[is.na(symbols)] <- clean_ids[is.na(symbols)]
      
      # Save successful mapping to cache
      all_symbols <- setNames(symbols, clean_ids)
      saveRDS(all_symbols, cache_file)
      message("  Gene symbols cached for future use.")
      
      message("  Conversion complete.")
      return(symbols)
      
    }, error = function(e) {
      if (attempt < max_retries) {
        message(sprintf("    Connection failed: %s. Retrying in 5 seconds...", e$message))
        Sys.sleep(5)
        NULL
      } else {
        message(sprintf("    All attempts failed: %s", e$message))
        "failed"
      }
    })
    
    if (!is.null(result) && result[1] != "failed") {
      return(result)
    }
  }
  
  # Final fallback
  message("  Using Ensembl IDs as fallback")
  return(ensembl_ids)
}

#-------------------------------------------------------
# 3. DESEQ2 ANALYSIS FUNCTION WITH LFC SHRINKAGE 
# MODIFIED TO USE RAW P-VALUES
#-------------------------------------------------------

run_wald_test <- function(counts, metadata, design_formula,
                          contrast = NULL, test_name = "test",
                          alpha = 0.05, lfc_threshold = 0,
                          apply_shrinkage = FALSE,
                          shrinkage_type = "apeglm") {
  
  message(sprintf("\n=== Running Wald Test: %s ===", test_name))
  message(sprintf(" Design: %s", deparse(design_formula)))
  message(" Using RAW P-VALUES (pvalue < 0.05) for significance")
  
  # Create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = design_formula
  )
  
  # Keep genes with at least 10 counts in at least 3 samples
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  message(" Genes after filtering: ", nrow(dds))
  
  # Run DESeq2 with error handling
  dds <- tryCatch({
    DESeq(dds,
          test = "Wald",
          fitType = "parametric",
          minReplicatesForReplace = 7,
          parallel = FALSE,
          quiet = TRUE)
  }, error = function(e) {
    message(" Parametric fit failed, using local regression")
    DESeq(dds,
          test = "Wald",
          fitType = "local",
          quiet = TRUE)
  })
  
  # Get the appropriate coefficient for results and shrinkage
  available_coefs <- resultsNames(dds)
  message(" Available coefficients: ", paste(available_coefs, collapse = ", "))
  
  # Determine which coefficient to use
  coef_name <- NULL
  use_contrast <- FALSE
  res <- NULL
  
  # Check if contrast is specified
  if (!is.null(contrast)) {
    res <- results(dds,
                   contrast = contrast,
                   alpha = alpha,
                   lfcThreshold = lfc_threshold,
                   independentFiltering = TRUE,
                   cooksCutoff = TRUE,
                   pAdjustMethod = "BH")
    use_contrast <- TRUE
    message("  Using contrast: ", paste(contrast, collapse = " "))
  } else {
    # For interaction terms, use the last coefficient
    coef_name <- available_coefs[length(available_coefs)]
    message("  Using coefficient: ", coef_name)
    res <- results(dds,
                   name = coef_name,
                   alpha = alpha,
                   lfcThreshold = lfc_threshold,
                   independentFiltering = TRUE,
                   cooksCutoff = TRUE,
                   pAdjustMethod = "BH")
  }
  
  # Apply LFC shrinkage properly
  res_shrunk <- NULL
  if (apply_shrinkage) {
    
    if (use_contrast) {
      message("  Applying ashr shrinkage (contrast detected)...")
      res_shrunk <- tryCatch({
        lfcShrink(dds,
                  contrast = contrast,
                  res = res,
                  type = "ashr",
                  quiet = TRUE)
      }, error = function(e) {
        message("    ashr shrinkage failed: ", e$message)
        NULL
      })
    } else if (!is.null(coef_name)) {
      
      # Check if it's an interaction term
      if (grepl(":", coef_name)) {
        message("  Applying ashr shrinkage (interaction term detected)...")
        res_shrunk <- tryCatch({
          lfcShrink(dds,
                    coef = coef_name,
                    res = res,
                    type = "ashr",
                    quiet = TRUE)
        }, error = function(e) {
          message("    ashr shrinkage failed: ", e$message)
          NULL
        })
      } else {
        # Simple coefficients: try apeglm first (preferred), then ashr
        message("  Attempting apeglm shrinkage for coefficient: ", coef_name)
        res_shrunk <- tryCatch({
          lfcShrink(dds,
                    coef = coef_name,
                    type = "apeglm",
                    quiet = TRUE)
        }, error = function(e) {
          message("    apeglm failed, trying ashr...")
          tryCatch({
            lfcShrink(dds,
                      coef = coef_name,
                      res = res,
                      type = "ashr",
                      quiet = TRUE)
          }, error = function(e2) {
            message("    ashr also failed: ", e2$message)
            NULL
          })
        })
      }
    }
    
    # Update results if shrinkage was successful
    if (!is.null(res_shrunk)) {
      message("  Shrinkage applied successfully")
      if (!"pvalue" %in% names(res_shrunk)) {
        res_shrunk$pvalue <- res$pvalue
        res_shrunk$padj <- res$padj
      }
      res <- res_shrunk
    } else {
      message("  Could not apply shrinkage, proceeding with unshrunken estimates")
    }
  }
  
  # Convert to data frame and add annotations
  res_df <- as.data.frame(res)
  res_df$ensembl_id <- rownames(res_df)
  res_df$gene_symbol <- convert_to_symbols(res_df$ensembl_id)
  
  # Add significance categories BASED ON RAW P-VALUE
  res_df <- res_df %>%
    mutate(
      significant = pvalue < alpha,
      direction = case_when(
        log2FoldChange > 0 ~ "up",
        log2FoldChange < 0 ~ "down",
        TRUE ~ "none"
      ),
      sig_category = case_when(
        is.na(pvalue) ~ "ns",
        pvalue == 0 ~ "***", 
        pvalue < 0.001 ~ "***",
        pvalue < 0.01 ~ "**",
        pvalue < 0.05 ~ "*",
        pvalue < 0.1 ~ ".",
        TRUE ~ "ns"
      )
    ) %>%
    arrange(pvalue)
  
  # Summary statistics BASED ON RAW P-VALUE
  sig_genes <- sum(res_df$significant, na.rm = TRUE)
  up_genes <- sum(res_df$significant & res_df$direction == "up", na.rm = TRUE)
  down_genes <- sum(res_df$significant & res_df$direction == "down", na.rm = TRUE)
  
  # Count genes with pvalue = 0
  zero_pvalue <- sum(res_df$pvalue == 0, na.rm = TRUE)
  if (zero_pvalue > 0) {
    message(sprintf(" Note: %d genes with pvalue = 0 (extremely significant)", zero_pvalue))
  }
  
  message(sprintf(" Significant genes (pvalue < 0.05): %d (up: %d, down: %d)",
                  sig_genes, up_genes, down_genes))
  
  # Also report padj for comparison
  sig_genes_padj <- sum(res_df$padj < 0.05, na.rm = TRUE)
  message(sprintf(" For comparison, padj < 0.05: %d genes", sig_genes_padj))
  
  # Save results to CSV
  write.csv(res_df,
            file.path(output_dir, paste0(test_name, "_results.csv")),
            row.names = FALSE)
  
  return(list(
    dds = dds,
    results = res_df,
    shrinkage_applied = !is.null(res_shrunk),
    summary = data.frame(
      test = test_name,
      total_genes = sum(!is.na(res_df$pvalue)),
      significant_pvalue = sig_genes,
      upregulated_pvalue = up_genes,
      downregulated_pvalue = down_genes,
      significant_padj = sig_genes_padj,
      zero_pvalue = zero_pvalue,
      shrinkage = ifelse(!is.null(res_shrunk), "Yes", "No")
    )
  ))
}

#--------------------------------------------------
# 4. MAIN ANALYSIS WORKFLOW
#--------------------------------------------------

# Create output directory
output_dir <- "wald_test_results_rawPvalue_noPC2_no_outlier_without_shrinkage"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Initialize results storage
full_results <- list()
full_summaries <- list()

message("\n")
message("========================================================")
message("  NO PC2 COVARIATE IN MODELS")
message("  USING RAW P-VALUES (pvalue < 0.05) FOR ALL ANALYSES")
message("========================================================")
message("\n")

# ===============================================================
# a. MAIN EFFECTS ANALYSIS
# ===============================================================

message("\n==== MAIN EFFECTS ANALYSIS ====")

# Treatment effect
full_results$treatment_main <- run_wald_test(
  counts, metadata,
  design_formula = ~ treatment,
  contrast = c("treatment", "THC", "vehicle"),
  test_name = "Treatment_Main_Effect_rawP",
  apply_shrinkage = FALSE
)

# Sex effect
full_results$sex_main <- run_wald_test(
  counts, metadata,
  design_formula = ~ sex,
  contrast = c("sex", "male", "female"),
  test_name = "Sex_Main_Effect_rawP",
  apply_shrinkage = FALSE
)

# Strain effect
full_results$strain_main <- run_wald_test(
  counts, metadata,
  design_formula = ~ strain,
  contrast = c("strain", "DBA_2J", "C57BL_6"),
  test_name = "Strain_Main_Effect_rawP",
  apply_shrinkage = FALSE
)

# ===============================================================
# b. TWO-WAY INTERACTION ANALYSIS
# ===============================================================

message("\n==== TWO-WAY INTERACTIONS ====")

# Treatment × Sex
full_results$treatment_sex_interaction <- run_wald_test(
  counts, metadata,
  design_formula = ~ treatment * sex,
  contrast = NULL,
  test_name = "Treatment_x_Sex_Interaction_rawP",
  apply_shrinkage = FALSE
)

# Treatment × Strain
full_results$treatment_strain_interaction <- run_wald_test(
  counts, metadata,
  design_formula = ~ treatment * strain,
  contrast = NULL,
  test_name = "Treatment_x_Strain_Interaction_rawP",
  apply_shrinkage = FALSE
)

# Sex × Strain
full_results$sex_strain_interaction <- run_wald_test(
  counts, metadata,
  design_formula = ~ sex * strain,
  contrast = NULL,
  test_name = "Sex_x_Strain_Interaction_rawP",
  apply_shrinkage = FALSE
)

# ===============================================================
# c. THREE-WAY INTERACTION ANALYSIS
# ===============================================================

message("\n==== THREE-WAY INTERACTION ====")

full_results$three_way_interaction <- run_wald_test(
  counts, metadata,
  design_formula = ~ treatment * sex * strain,
  contrast = NULL,
  test_name = "Three_Way_Interaction_rawP",
  apply_shrinkage = FALSE
)

# ===============================================================
# d. STRATIFIED ANALYSES
# ===============================================================

message("\n==== STRATIFIED ANALYSES ====")
message("\n==== SEX-BASED ANALYSES ====")

# Treatment effect in females
female_idx <- metadata$sex == "female"
full_results$treatment_effect_in_all_females <- run_wald_test(
  counts[, female_idx],
  metadata[female_idx, ],
  design_formula = ~ treatment,
  contrast = c("treatment", "THC", "vehicle"),
  test_name = "Treatment_in_females_rawP",
  apply_shrinkage = FALSE
)

message("\n==== TREATMENT/STRAIN INTERACTION IN FEMALES ====")

# Treatment and Strain interaction in females
full_results$treatment_strain_interaction_in_all_females <- run_wald_test(
  counts[, female_idx],
  metadata[female_idx, ],
  design_formula = ~ treatment * strain,
  contrast = NULL,
  test_name = "Treatment_strain_interaction_in_females_rawP",
  apply_shrinkage = FALSE
)

# Treatment effect in males
male_idx <- metadata$sex == "male"
full_results$treatment_effect_in_all_males <- run_wald_test(
  counts[, male_idx],
  metadata[male_idx, ],
  design_formula = ~ treatment,
  contrast = c("treatment", "THC", "vehicle"),
  test_name = "Treatment_in_males_rawP",
  apply_shrinkage = FALSE
)

message("\n==== TREATMENT/STRAIN INTERACTION IN MALES ====")

# Treatment and Strain interaction in males
full_results$treatment_strain_interaction_in_all_males <- run_wald_test(
  counts[, male_idx],
  metadata[male_idx, ],
  design_formula = ~ treatment * strain,
  contrast = NULL,
  test_name = "Treatment_strain_interaction_in_males_rawP",
  apply_shrinkage = FALSE
)

message("\n==== STRAIN-BASED ANALYSES ====")

# Treatment effect in C57BL_6
c57_idx <- metadata$strain == "C57BL_6"
full_results$treatment_effect_in_all_C57BL6 <- run_wald_test(
  counts[, c57_idx],
  metadata[c57_idx, ],
  design_formula = ~ treatment,
  contrast = c("treatment", "THC", "vehicle"),
  test_name = "Treatment_in_C57BL6_rawP",
  apply_shrinkage = FALSE
)

# Treatment/sex interaction in C57BL_6
full_results$treatment_sex_interaction_in_all_C57BL6 <- run_wald_test(
  counts[, c57_idx],
  metadata[c57_idx, ],
  design_formula = ~ treatment * sex,
  contrast = NULL,
  test_name = "Treatment_sex_interaction_in_C57BL6_rawP",
  apply_shrinkage = FALSE
)

# Treatment effect in DBA_2J
dba_idx <- metadata$strain == "DBA_2J"
full_results$treatment_effect_in_all_DBA2J <- run_wald_test(
  counts[, dba_idx],
  metadata[dba_idx, ],
  design_formula = ~ treatment,
  contrast = c("treatment", "THC", "vehicle"),
  test_name = "Treatment_in_DBA2J_rawP",
  apply_shrinkage = FALSE
)

# Treatment/sex interaction in DBA_2J
full_results$treatment_sex_interaction_in_all_DBA2J <- run_wald_test(
  counts[, dba_idx],
  metadata[dba_idx, ],
  design_formula = ~ treatment * sex,
  contrast = NULL,
  test_name = "Treatment_sex_interaction_in_DBA2J_rawP",
  apply_shrinkage = FALSE
)

message("\n==== TREATMENT EFFECT BY STRAIN & SEX ====")

# Treatment effect in C57BL_6 Females
female_c57_idx <- metadata$strain == "C57BL_6" & metadata$sex == "female"
full_results$treatment_effect_in_female_C57BL6 <- run_wald_test(
  counts[, female_c57_idx],
  metadata[female_c57_idx, ],
  design_formula = ~ treatment,
  contrast = c("treatment", "THC", "vehicle"),
  test_name = "Treatment_in_female_C57BL6_rawP",
  apply_shrinkage = FALSE
)

# Treatment effect in DBA_2J Females
female_dba_idx <- metadata$strain == "DBA_2J" & metadata$sex == "female"
full_results$treatment_effect_in_female_DBA2J <- run_wald_test(
  counts[, female_dba_idx],
  metadata[female_dba_idx, ],
  design_formula = ~ treatment,
  contrast = c("treatment", "THC", "vehicle"),
  test_name = "Treatment_in_female_DBA2J_rawP",
  apply_shrinkage = FALSE
)

# Treatment effect in C57BL_6 Males
male_c57_idx <- metadata$strain == "C57BL_6" & metadata$sex == "male"
full_results$treatment_effect_in_male_C57BL6 <- run_wald_test(
  counts[, male_c57_idx],
  metadata[male_c57_idx, ],
  design_formula = ~ treatment,
  contrast = c("treatment", "THC", "vehicle"),
  test_name = "Treatment_in_male_C57BL6_rawP",
  apply_shrinkage = FALSE
)

# Treatment effect in DBA_2J Males
male_dba_idx <- metadata$strain == "DBA_2J" & metadata$sex == "male"
full_results$treatment_effect_in_male_DBA2J <- run_wald_test(
  counts[, male_dba_idx],
  metadata[male_dba_idx, ],
  design_formula = ~ treatment,
  contrast = c("treatment", "THC", "vehicle"),
  test_name = "Treatment_in_male_DBA2J_rawP",
  apply_shrinkage = FALSE
)

#-----------------------------------------------#
# 5. PCA Visualization (USING RAW P-VALUES)
#-----------------------------------------------#

create_interaction_pca <- function(dds, interaction_name, highlight_genes = NULL) {
  vsd <- vst(dds, blind = FALSE)
  
  # Feature selection BASED ON RAW P-VALUES
  if (!is.null(highlight_genes)) {
    gene_subset <- intersect(highlight_genes, rownames(vsd))
    if (length(gene_subset) < 50) {
      rv <- rowVars(assay(vsd))
      additional <- setdiff(names(sort(rv, decreasing = TRUE)[1:500]), gene_subset)
      gene_subset <- c(gene_subset, additional[1:min(450, length(additional))])
    }
  } else {
    rv <- rowVars(assay(vsd))
    gene_subset <- names(sort(rv, decreasing = TRUE)[1:500])
  }
  
  # PCA
  pca <- prcomp(t(assay(vsd)[gene_subset, ]))
  
  # Assemble data
  pca_data <- as.data.frame(pca$x[, 1:4])
  pca_data$sample_name <- rownames(pca_data)
  
  # Get metadata
  metadata_cols <- as.data.frame(colData(vsd))
  pca_data <- cbind(pca_data, metadata_cols)
  
  var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  
  list(
    pca_data = pca_data,
    var_explained = var_explained,
    pca_object = pca,
    gene_subset = gene_subset
  )
}

# PCA for each interaction model - USING RAW P-VALUES
pca_results <- list()
for (name in names(full_results)) {
  sig_genes <- full_results[[name]]$results %>%
    filter(!is.na(pvalue) & pvalue <= 0.05) %>%
    pull(ensembl_id)
  
  if (length(sig_genes) >= 1) {
    pca_results[[name]] <- create_interaction_pca(
      full_results[[name]]$dds,
      name,
      highlight_genes = sig_genes[1:min(100, length(sig_genes))]
    )
  }
}

# ---------------------------------------------------------------
# 6. PCA PLOTS / VOLCANO PLOT / HEATMAP
# ---------------------------------------------------------------

# Function for individual full-size PCA plots
create_individual_pca_plot <- function(pca_data, var_explained, title,
                                       color_by, shape_by = NULL,
                                       show_labels = TRUE, label_size = 3) {
  p <- ggplot(pca_data, aes_string(x = "PC1", y = "PC2", color = color_by)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(level = 0.95, linetype = "dashed", alpha = 0.5) +
    scale_color_brewer(palette = "Set1") +
    xlab(paste0("PC1: ", var_explained[1], "%")) +
    ylab(paste0("PC2: ", var_explained[2], "%")) +
    ggtitle(title) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
  
  if (!is.null(shape_by) && shape_by %in% names(pca_data)) {
    p <- p + aes_string(shape = shape_by)
  }
  
  if (show_labels) {
    p <- p + geom_text_repel(
      aes(label = sample_name),
      size = label_size,
      max.overlaps = 20,
      box.padding = 0.35,
      point.padding = 0.3,
      segment.color = 'grey50',
      segment.size = 0.2,
      alpha = 0.8
    )
  }
  p
}

# Function for combined multi-panel PCA plots
create_combined_pca_panel <- function(pca_data, var_explained, title,
                                      color_by, shape_by = NULL,
                                      show_labels = TRUE, label_size = 2.5) {
  p <- ggplot(pca_data, aes_string(x = "PC1", y = "PC2", color = color_by)) +
    geom_point(size = 3.5, alpha = 0.8) +
    stat_ellipse(level = 0.95, linetype = "dashed", alpha = 0.5) +
    scale_color_brewer(palette = "Set1") +
    xlab(paste0("PC1: ", var_explained[1], "%")) +
    ylab(paste0("PC2: ", var_explained[2], "%")) +
    ggtitle(title) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9)
    )
  
  if (!is.null(shape_by) && shape_by %in% names(pca_data)) {
    p <- p + aes_string(shape = shape_by)
  }
  
  if (show_labels) {
    p <- p + geom_text_repel(
      aes(label = sample_name),
      size = label_size,
      max.overlaps = 20,
      box.padding = 0.3,
      point.padding = 0.25,
      segment.color = 'grey50',
      segment.size = 0.2,
      alpha = 0.8
    )
  }
  p
}

# Directory for individual plots
individual_plots_dir <- file.path(output_dir, "individual_pca_plots")
dir.create(individual_plots_dir, showWarnings = FALSE)

# Create lists to store plots
individual_pca_plots <- list()
combined_pca_plots <- list()

# Generate individual and combined plots
for (model_name in names(pca_results)) {
  pca_data <- pca_results[[model_name]]$pca_data
  var_explained <- pca_results[[model_name]]$var_explained
  
  # Determine color and shape variables based on model
  if (grepl("treatment.*sex", model_name)) {
    color_var <- "treatment"; shape_var <- "sex"
  } else if (grepl("treatment.*strain", model_name)) {
    color_var <- "treatment"; shape_var <- "strain"
  } else if (grepl("sex.*strain", model_name)) {
    color_var <- "sex"; shape_var <- "strain"
  } else if (grepl("treatment", model_name)) {
    color_var <- "treatment"; shape_var <- NULL
  } else if (grepl("sex", model_name)) {
    color_var <- "sex"; shape_var <- NULL
  } else if (grepl("strain", model_name)) {
    color_var <- "strain"; shape_var <- NULL
  } else {
    color_var <- "treatment"; shape_var <- "sex"
  }
  
  # Individual full-size plot
  individual_plot <- create_individual_pca_plot(
    pca_data, var_explained,
    paste("PCA:", gsub("_rawP", " (raw p-value)", gsub("_", " ", model_name))),
    color_var, shape_var,
    show_labels = TRUE,
    label_size = 3
  )
  
  # Save individual plot
  ggsave(file.path(individual_plots_dir, paste0(model_name, "_pca_with_labels.png")),
         individual_plot,
         width = 10,
         height = 8,
         dpi = 300)
  
  individual_pca_plots[[model_name]] <- individual_plot
  
  # Plot for combined panel
  combined_plot <- create_combined_pca_panel(
    pca_data, var_explained,
    paste(gsub("_rawP", " (raw p)", gsub("_", " ", model_name))),
    color_var, shape_var,
    show_labels = TRUE,
    label_size = 2.2
  )
  
  combined_pca_plots[[model_name]] <- combined_plot
}

# Create combined multi-panel figure
if (length(combined_pca_plots) > 0) {
  n_plots <- length(combined_pca_plots)
  ncol_layout <- if(n_plots <= 2) 2 else if(n_plots <= 4) 2 else if(n_plots <= 6) 3 else 4
  nrow_layout <- ceiling(n_plots / ncol_layout)
  
  combined_pca <- wrap_plots(combined_pca_plots, ncol = ncol_layout, nrow = nrow_layout) +
    plot_annotation(
      title = "Gene Expression PCA Analysis: All Models (Raw P-values)",
      subtitle = "Based on significant genes (pvalue ≤ 0.05) from each model | No PC2 covariate",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5)
      )
    )
  
  output_filename <- file.path(output_dir, "combined_pca_all_models_rawPvalue.pdf")
  
  plot_width <- min(30, ncol_layout * 7)
  plot_height <- min(30, nrow_layout * 6)
  
  ggsave(output_filename,
         combined_pca,
         width = plot_width,
         height = plot_height,
         dpi = 300,
         limitsize = FALSE)
}

# ---------------------------------------------------------------
# 6.5 VOLCANO PLOTS (USING RAW P-VALUES)
# ---------------------------------------------------------------

create_volcano_plot <- function(res_df, title, my_top = 10, 
                                lfc_threshold = 1, p_threshold = 0.05) {
  
  total_genes <- nrow(res_df)
  
  volcano_data <- res_df %>%
    filter(!is.na(pvalue) & !is.na(log2FoldChange)) %>%
    mutate(
      sig_status = case_when(
        pvalue < p_threshold & log2FoldChange > lfc_threshold ~ "Up",
        pvalue < p_threshold & log2FoldChange < -lfc_threshold ~ "Down",
        pvalue < p_threshold ~ "Sig_no_LFC",
        TRUE ~ "NS"
      ),
      neg_log10_pvalue = pmin(-log10(pvalue + 1e-300), 300)
    )
  
  plotted_genes <- nrow(volcano_data)
  
  top_genes <- volcano_data %>%
    filter(sig_status %in% c("Up", "Down")) %>%
    arrange(pvalue) %>%
    slice_head(n = my_top) %>%
    pull(gene_symbol)
  
  p <- ggplot(volcano_data, aes(x = log2FoldChange, y = neg_log10_pvalue)) +
    geom_point(aes(color = sig_status), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c(
      "Up" = "red",
      "Down" = "blue",
      "Sig_no_LFC" = "orange",
      "NS" = "grey60"
    ),
    labels = c(
      "Up" = "Up-regulated",
      "Down" = "Down-regulated",
      "Sig_no_LFC" = "Sig (|LFC| < 1)",
      "NS" = "Not significant"
    )) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), 
               linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(p_threshold),
               linetype = "dashed", alpha = 0.5) +
    geom_text_repel(
      data = filter(volcano_data, gene_symbol %in% top_genes),
      aes(label = gene_symbol),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = 'grey50',
      segment.size = 0.2
    ) +
    labs(
      x = "log2 Fold Change",
      y = "-log10(p-value)",
      title = title,
      subtitle = sprintf("Sig genes: %d up, %d down | %d/%d genes plotted (raw p < 0.05)",
                         sum(volcano_data$sig_status == "Up"),
                         sum(volcano_data$sig_status == "Down"),
                         plotted_genes, total_genes),
      color = "Status"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "right"
    )
  
  return(p)
}

# Create volcano plots directory
volcano_plots_dir <- file.path(output_dir, "volcano_plots")
dir.create(volcano_plots_dir, showWarnings = FALSE)

# Generate volcano plots - USING RAW P-VALUES
volcano_plots <- list()
for (model_name in names(full_results)) {
  res_df <- full_results[[model_name]]$results
  n_sig <- sum(res_df$pvalue < 0.05, na.rm = TRUE)
  
  if (n_sig >= 5) {
    volcano_plot <- create_volcano_plot(
      res_df,
      title = paste("Volcano:", gsub("_rawP", " (raw p)", gsub("_", " ", model_name))),
      my_top = 10,
      lfc_threshold = 1,
      p_threshold = 0.05
    )
    
    ggsave(
      file.path(volcano_plots_dir, paste0(model_name, "_volcano.png")),
      volcano_plot,
      width = 8,
      height = 6,
      dpi = 300
    )
    
    volcano_plots[[model_name]] <- volcano_plot
  }
}

# Combined volcano plot
if (length(volcano_plots) > 1) {
  n_plots <- length(volcano_plots)
  ncol_volcano <- if(n_plots <= 2) 2 else if(n_plots <= 4) 2 else 3
  nrow_volcano <- ceiling(n_plots / ncol_volcano)
  
  combined_volcano <- wrap_plots(volcano_plots, ncol = ncol_volcano, nrow = nrow_volcano) +
    plot_annotation(
      title = "Volcano Plots: Differential Expression Analysis (Raw P-values)",
      subtitle = "Top 10 DEGs labeled | Based on raw p-values (p < 0.05) | No PC2 covariate",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )
  
  ggsave(
    file.path(output_dir, "combined_volcano_plots_rawPvalue.pdf"),
    combined_volcano,
    width = min(24, ncol_volcano * 8),
    height = min(24, nrow_volcano * 6),
    dpi = 300,
    limitsize = FALSE
  )
  
  message("\nVolcano plots generated for ", length(volcano_plots), " models (using raw p-values)")
}

# ---------------------------------------------------------------
# 6.6 HEATMAPS (USING RAW P-VALUES)
# ---------------------------------------------------------------

if (!require(pheatmap)) {
  install.packages("pheatmap")
  library(pheatmap)
}

create_deg_heatmap <- function(dds, res_df, title, my_top = 20) {
  
  # Rank by pvalue
  top_genes <- res_df %>%
    filter(!is.na(pvalue)) %>%
    arrange(pvalue) %>%
    slice_head(n = my_top) %>%
    pull(ensembl_id)
  
  if (length(top_genes) < 2) {
    message("  Not enough genes for heatmap")
    return(NULL)
  }
  
  vsd <- vst(dds, blind = FALSE)
  mat <- assay(vsd)[top_genes, ]
  mat_scaled <- t(scale(t(mat)))
  
  gene_symbols <- res_df %>%
    filter(ensembl_id %in% top_genes) %>%
    arrange(match(ensembl_id, top_genes)) %>%
    pull(gene_symbol)
  
  rownames(mat_scaled) <- gene_symbols
  
  sample_annotation <- as.data.frame(colData(dds)[, c("treatment", "sex", "strain")])
  
  ann_colors <- list(
    treatment = c(vehicle = "blue", THC = "red"),
    sex = c(female = "pink", male = "lightblue"),
    strain = c(C57BL_6 = "green", DBA_2J = "orange")
  )
  
  return(list(
    mat = mat_scaled,
    annotation_col = sample_annotation,
    annotation_colors = ann_colors,
    title = title
  ))
}

# Create heatmaps directory
heatmaps_dir <- file.path(output_dir, "heatmaps")
dir.create(heatmaps_dir, showWarnings = FALSE)

# Generate heatmaps - USING RAW P-VALUES
for (model_name in names(full_results)) {
  res_df <- full_results[[model_name]]$results
  n_sig <- sum(res_df$pvalue < 0.05, na.rm = TRUE)
  
  if (n_sig >= 1) {
    message(sprintf("\nCreating heatmap for %s (ranked by raw p-value)...", model_name))
    
    heatmap_data <- tryCatch({
      create_deg_heatmap(
        full_results[[model_name]]$dds,
        res_df,
        title = paste("Top 20 DEGs:", gsub("_rawP", " (raw p)", gsub("_", " ", model_name))),
        my_top = 20
      )
    }, error = function(e) {
      message(sprintf("  Error: %s", e$message))
      NULL
    })
    
    if (!is.null(heatmap_data)) {
      pdf(file.path(heatmaps_dir, paste0(model_name, "_heatmap.pdf")), 
          width = 10, height = 8)
      
      pheatmap(
        heatmap_data$mat,
        main = heatmap_data$title,
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        clustering_method = "complete",
        annotation_col = heatmap_data$annotation_col,
        annotation_colors = heatmap_data$annotation_colors,
        show_colnames = TRUE,
        show_rownames = TRUE,
        fontsize_row = 8,
        fontsize_col = 8,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        breaks = seq(-2, 2, length.out = 101),
        border_color = NA
      )
      
      dev.off()
      
      png(file.path(heatmaps_dir, paste0(model_name, "_heatmap.png")), 
          width = 10, height = 8, units = "in", res = 300)
      
      pheatmap(
        heatmap_data$mat,
        main = heatmap_data$title,
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        clustering_method = "complete",
        annotation_col = heatmap_data$annotation_col,
        annotation_colors = heatmap_data$annotation_colors,
        show_colnames = TRUE,
        show_rownames = TRUE,
        fontsize_row = 8,
        fontsize_col = 8,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        breaks = seq(-2, 2, length.out = 101),
        border_color = NA
      )
      
      dev.off()
      
      message(sprintf("  Saved %s heatmaps", model_name))
    }
  }
}

# ---------------------------------------------------------------
# 7. SUMMARY PLOTS AND REPORTS (USING RAW P-VALUES)
# ---------------------------------------------------------------

# DEG summary bar plot - USING RAW P-VALUES
deg_summary <- data.frame(
  Model = names(full_results),
  Up = sapply(full_results, function(x) {
    sum(x$results$pvalue < 0.05 & x$results$log2FoldChange > 0, na.rm = TRUE)
  }),
  Down = sapply(full_results, function(x) {
    sum(x$results$pvalue < 0.05 & x$results$log2FoldChange < 0, na.rm = TRUE)
  })
) %>%
  pivot_longer(cols = c(Up, Down), names_to = "Direction", values_to = "Count") %>%
  mutate(Model = factor(Model, levels = unique(Model)))

deg_barplot <- ggplot(deg_summary, aes(x = Model, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
  labs(
    title = "Differentially Expressed Genes Across All Models (Raw P-values)",
    subtitle = "Based on raw p-value < 0.05 | No PC2 covariate",
    x = "Model",
    y = "Number of DEGs (p < 0.05)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "top"
  ) +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3)

ggsave(file.path(output_dir, "DEG_summary_barplot_rawPvalue.pdf"),
       deg_barplot, width = 14, height = 6, dpi = 300)

ggsave(file.path(output_dir, "DEG_summary_barplot_rawPvalue.png"),
       deg_barplot, width = 14, height = 6, dpi = 300)

# ---------------------------------------------------------------
# 8. FINAL SUMMARY REPORT (WITH BOTH RAW P-VALUE AND PADJ COUNTS)
# ---------------------------------------------------------------

summary_data <- data.frame(
  Model = names(full_results),
  Total_Genes_Tested = sapply(full_results, function(x) sum(!is.na(x$results$pvalue))),
  # Raw p-value counts
  Significant_pvalue_0.05 = sapply(full_results, function(x) sum(x$results$pvalue < 0.05, na.rm = TRUE)),
  Significant_pvalue_0.01 = sapply(full_results, function(x) sum(x$results$pvalue < 0.01, na.rm = TRUE)),
  Significant_pvalue_0.001 = sapply(full_results, function(x) sum(x$results$pvalue < 0.001, na.rm = TRUE)),
  # Adjusted p-value counts for comparison
  Significant_padj_0.05 = sapply(full_results, function(x) sum(x$results$padj < 0.05, na.rm = TRUE)),
  Significant_padj_0.01 = sapply(full_results, function(x) sum(x$results$padj < 0.01, na.rm = TRUE)),
  Shrinkage_Applied = sapply(full_results, function(x) x$shrinkage_applied),
  PC2_Covariate = "No",
  Analysis_Type = "Raw_P_values"
)

sample_summary <- metadata %>%
  group_by(treatment, sex, strain) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  arrange(treatment, sex, strain)

write.csv(summary_data,
          file.path(output_dir, "analysis_summary_rawPvalue.csv"),
          row.names = FALSE)

write.csv(sample_summary,
          file.path(output_dir, "sample_group_summary.csv"),
          row.names = FALSE)

# Print summaries
message("\n")
message("========================================================")
message("  ANALYSIS COMPLETE - RAW P-VALUES")
message("========================================================")
print(summary_data)
message("\n--- Sample Summary ---")
print(sample_summary)

message("\n========================================")
message("Results saved in: ", output_dir)
message("========================================")
message("\nKey outputs:")
message(" - Analysis uses RAW P-VALUES (p < 0.05) for all filtering")
message(" - NO PC2 covariate in any models")
message(" - All result CSVs include BOTH pvalue and padj columns")
message(" - All plots based on raw p-values")
message(" - Combined PCA: combined_pca_all_models_rawPvalue.pdf")
message(" - Individual PCAs in: ", individual_plots_dir)
message(" - Volcano plots in: ", volcano_plots_dir)
message(" - Heatmaps in: ", heatmaps_dir)
message(" - DEG summary: DEG_summary_barplot_rawPvalue.pdf")