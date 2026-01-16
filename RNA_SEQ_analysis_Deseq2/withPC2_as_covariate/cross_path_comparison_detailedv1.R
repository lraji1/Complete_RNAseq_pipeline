#------------------------------------------------------------------------------
# DETAILED CROSS-MODEL GENE COMPARISON
#------------------------------------------------------------------------------

message("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
message("CREATING CROSS-MODEL GENE COMPARISONS")
message(paste(rep("=", 60), collapse = ""), "\n", sep = "")

# Create directory for cross-model gene analysis
cross_model_genes_dir <- file.path(fgsea_dir, "cross_model_genes")
dir.create(cross_model_genes_dir, showWarnings = FALSE, recursive = TRUE)

# Collect all pathways and their leading edge genes across all models and databases
message("\nCollecting pathways and leading edge genes across all models and databases...")

# Initialize data structures
all_pathway_genes <- list()
pathway_model_nes <- list()
pathway_model_padj <- list()

# Iterate through all models and collections
for (model_name in names(full_fgsea_results)) {
  message(" Processing model: ", model_name)
  
  for (collection_name in names(full_fgsea_results[[model_name]])) {
    if (collection_name %in% c("celltype", "micro_RNA" , "pertubations" , "immune_sig")) next
    
    fgsea_res <- full_fgsea_results[[model_name]][[collection_name]]
    
    if (!is.null(fgsea_res) && nrow(fgsea_res) > 0) {
      # Filter for significant pathways by abs_NES > 1.4 & padj <= 0.25
      sig_pathways <- fgsea_res %>%
        filter(abs_NES > 1.4 & padj <= 0.25)
      
      if (nrow(sig_pathways) > 0) {
        for (i in seq_len(nrow(sig_pathways))) {
          pathway_name <- sig_pathways$pathway[i]
          pathway_key <- paste0(collection_name, "::", pathway_name)
          
          # Store leading edge genes
          le_genes <- unlist(sig_pathways$leadingEdge[i], use.names = FALSE)
          
          if (length(le_genes) > 0) {
            if (!pathway_key %in% names(all_pathway_genes)) {
              all_pathway_genes[[pathway_key]] <- list()
              pathway_model_nes[[pathway_key]] <- list()
              pathway_model_padj[[pathway_key]] <- list()
            }
            
            all_pathway_genes[[pathway_key]][[model_name]] <- le_genes
            
            # Store NES and padj values
            pathway_model_nes[[pathway_key]][[model_name]] <- sig_pathways$NES[i]
            pathway_model_padj[[pathway_key]][[model_name]] <- sig_pathways$padj[i]
          }
        }
      }
    }
  }
}

message(" Total pathways collected: ", length(all_pathway_genes))

#------------------------------------------------------------------------------
#  Filter genes that appear in ≥2 pathways (BEFORE filtering pathways)
#------------------------------------------------------------------------------

message("\n[STEP 2] Identifying genes that appear in multiple pathways...")

# Create a gene-to-pathways mapping from ALL pathways (not filtered yet)
gene_pathway_map <- list()
for (pathway_key in names(all_pathway_genes)) {
  for (model_name in names(all_pathway_genes[[pathway_key]])) {
    genes <- all_pathway_genes[[pathway_key]][[model_name]]
    for (gene in genes) {
      if (!gene %in% names(gene_pathway_map)) {
        gene_pathway_map[[gene]] <- list()
      }
      # Store pathway and model information
      gene_pathway_map[[gene]][[length(gene_pathway_map[[gene]]) + 1]] <-
        list(pathway = pathway_key, model = model_name)
    }
  }
}

# Filter genes that appear in at least 2 pathways
multi_pathway_genes <- list()
for (gene in names(gene_pathway_map)) {
  unique_pathways <- unique(sapply(gene_pathway_map[[gene]], function(x) x$pathway))
  if (length(unique_pathways) >= 2) {
    multi_pathway_genes[[gene]] <- gene_pathway_map[[gene]]
  }
}

message(" Genes appearing in 2+ pathways: ", length(multi_pathway_genes))

# Create gene frequency summary
gene_pathway_counts <- sapply(names(multi_pathway_genes), function(gene) {
  length(unique(sapply(multi_pathway_genes[[gene]], function(x) x$pathway)))
})

# Sort genes by frequency
gene_pathway_counts <- sort(gene_pathway_counts, decreasing = TRUE)

# Save gene frequency table
gene_freq_df <- data.frame(
  gene = names(gene_pathway_counts),
  pathway_count = gene_pathway_counts,
  stringsAsFactors = FALSE
)
write.csv(gene_freq_df,
          file.path(cross_model_genes_dir, "gene_pathway_frequency.csv"),
          row.names = FALSE)
message(" Saved gene frequency table")

#------------------------------------------------------------------------------
#  Filter pathways that appear in ≥2 models (AFTER filtering genes)
#------------------------------------------------------------------------------

message("\n[STEP 1] Filtering pathways enriched in ≥2 models...")

# Get all pathways that contain at least one multi-pathway gene
relevant_pathways <- unique(unlist(lapply(names(multi_pathway_genes), function(gene) {
  unique(sapply(multi_pathway_genes[[gene]], function(x) x$pathway))
})))

# Filter these pathways to those enriched in ≥2 models
recurrent_pathways <- relevant_pathways[sapply(relevant_pathways, function(pathway) {
  length(names(all_pathway_genes[[pathway]])) >= 2
})]

message(" Pathways (with multi-pathway genes) enriched in >=2 models: ", length(recurrent_pathways))

#------------------------------------------------------------------------------
# CREATE VISUALIZATIONS
#------------------------------------------------------------------------------

if (length(multi_pathway_genes) > 0 && length(recurrent_pathways) > 0) {
  message("\nCreating heatmaps for ", length(multi_pathway_genes), " multi-pathway genes...")
  message(" Using ", length(recurrent_pathways), " pathways enriched in >=2 models")
  
  # Create NES matrix for pathways across models
  models <- names(full_fgsea_results)
  
  # Initialize matrices
  nes_matrix <- matrix(0,
                       nrow = length(recurrent_pathways),
                       ncol = length(models),
                       dimnames = list(recurrent_pathways, models))
  
  padj_matrix <- matrix(1,
                        nrow = length(recurrent_pathways),
                        ncol = length(models),
                        dimnames = list(recurrent_pathways, models))
  
  # Fill matrices
  for (pathway_key in recurrent_pathways) {
    for (model_name in models) {
      if (model_name %in% names(pathway_model_nes[[pathway_key]])) {
        nes_matrix[pathway_key, model_name] <- pathway_model_nes[[pathway_key]][[model_name]]
        padj_matrix[pathway_key, model_name] <- pathway_model_padj[[pathway_key]][[model_name]]
      }
    }
  }
  
  # Create display matrix with significance markers
  display_matrix <- matrix("",
                           nrow = nrow(nes_matrix),
                           ncol = ncol(nes_matrix),
                           dimnames = dimnames(nes_matrix))
  
  for (i in 1:nrow(nes_matrix)) {
    for (j in 1:ncol(nes_matrix)) {
      nes_val <- round(nes_matrix[i, j], 2)
      padj_val <- padj_matrix[i, j]
      
      if (padj_val <= 0.001) {
        display_matrix[i, j] <- paste0(nes_val, "***")
      } else if (padj_val <= 0.01) {
        display_matrix[i, j] <- paste0(nes_val, "**")
      } else if (padj_val <= 0.05) {
        display_matrix[i, j] <- paste0(nes_val, "*")
      } else {
        display_matrix[i, j] <- as.character(nes_val)
      }
    }
  }
  
  # Clean pathway names for display
  pathway_labels <- rownames(nes_matrix)
  pathway_labels <- gsub("^[^:]+::", "", pathway_labels)
  pathway_labels <- ifelse(nchar(pathway_labels) > 60,
                           paste0(substr(pathway_labels, 1, 57), "..."),
                           pathway_labels)
  rownames(nes_matrix) <- make.unique(pathway_labels, sep = "_")
  rownames(display_matrix) <- rownames(nes_matrix)
  
  # Create main heatmap: Pathways × Models (ALL PATHWAYS)
  if (sum(nes_matrix != 0) > 0) {
    out_png <- file.path(cross_model_genes_dir, "pathways_cross_model_NES_heatmap.png")
    
    # Calculate dimensions in inches
    width_in <- max(8, ncol(nes_matrix) * 1.5 + 6)
    height_in <- max(10, nrow(nes_matrix) * 0.25 + 3)
    
    # Cap maximum dimensions
    width_in <- min(width_in, 40)
    height_in <- min(height_in, 100)
    
    message(" Plot size: ", round(width_in, 1), " x ", round(height_in, 1), " inches")
    
    # Adjust font sizes based on number of items
    row_fontsize <- max(4, min(10, 250/nrow(nes_matrix)))
    col_fontsize <- max(8, min(12, 100/ncol(nes_matrix)))
    num_fontsize <- max(3, min(7, 150/nrow(nes_matrix)))
    
    tryCatch({
      png(out_png, width = width_in, height = height_in, units = "in", res = 300)
      
      pheatmap(nes_matrix,
               color = colorRampPalette(c("blue", "white", "red"))(100),
               breaks = seq(-3, 3, length.out = 101),
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               main = "Pathways (≥2 Models) × Models\nContaining Genes in ≥2 Pathways\n* p≤0.05, ** p≤0.01, *** p≤0.001",
               fontsize_row = row_fontsize,
               fontsize_col = col_fontsize,
               fontsize = 11,
               display_numbers = display_matrix,
               number_format = "%s",
               number_color = "black",
               fontsize_number = num_fontsize,
               border_color = NA,
               na_col = "grey95",
               legend_title = "NES")
      
      dev.off()
      message(" Created pathway NES heatmap (all pathways)")
      
    }, error = function(e) {
      message(" PNG device failed: ", e$message)
      message(" Attempting PDF alternative...")
      
      # Try PDF as fallback
      out_pdf <- file.path(cross_model_genes_dir, "pathways_cross_model_NES_heatmap.pdf")
      
      tryCatch({
        pdf(out_pdf, width = width_in, height = height_in)
        
        pheatmap(nes_matrix,
                 color = colorRampPalette(c("blue", "white", "red"))(100),
                 breaks = seq(-3, 3, length.out = 101),
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 main = "All Pathways Common in ≥2 Models with Models\nContaining Genes in ≥2 Pathways\n* p≤0.05, ** p≤0.01, *** p≤0.001",
                 fontsize_row = row_fontsize,
                 fontsize_col = col_fontsize,
                 fontsize = 11,
                 display_numbers = display_matrix,
                 number_format = "%s",
                 number_color = "black",
                 fontsize_number = num_fontsize,
                 border_color = NA,
                 na_col = "grey95",
                 legend_title = "NES")
        
        dev.off()
        message(" Created pathway heatmap as PDF: ", out_pdf)
        
      }, error = function(e2) {
        message(" PDF also failed: ", e2$message)
        message(" Consider subsetting the data or increasing system memory")
      })
    })
  }
  
  #------------------------------------------------------------------------------
  #  Create TOP 30 PATHWAYS heatmap
  #------------------------------------------------------------------------------
  
  message("\nCreating TOP 30 pathways heatmap...")
  
  # Calculate a ranking metric for pathways (mean absolute NES across models)
  pathway_scores <- apply(nes_matrix, 1, function(row) {
    mean(abs(row[row != 0]))  # mean of non-zero absolute NES values
  })
  
  # Sort and select top 30 pathways
  top_pathways_idx <- order(pathway_scores, decreasing = TRUE)[1:min(30, length(pathway_scores))]
  
  # Subset matrices for top 30 pathways
  nes_matrix_top30 <- nes_matrix[top_pathways_idx, , drop = FALSE]
  display_matrix_top30 <- display_matrix[top_pathways_idx, , drop = FALSE]
  
  message(" Selected top ", nrow(nes_matrix_top30), " pathways")
  
  if (sum(nes_matrix_top30 != 0) > 0) {
    out_png_top30 <- file.path(cross_model_genes_dir, "pathways_cross_model_NES_heatmap_top30.png")
    
    # Calculate dimensions for top 30
    width_in_top30 <- max(8, ncol(nes_matrix_top30) * 1.5 + 6)
    height_in_top30 <- max(10, nrow(nes_matrix_top30) * 0.4 + 3)
    
    # Cap dimensions
    width_in_top30 <- min(width_in_top30, 40)
    height_in_top30 <- min(height_in_top30, 50)
    
    message(" Top 30 plot size: ", round(width_in_top30, 1), " x ", round(height_in_top30, 1), " inches")
    
    # Adjust font sizes for top 30
    row_fontsize_top30 <- max(6, min(10, 300/nrow(nes_matrix_top30)))
    col_fontsize_top30 <- max(8, min(12, 100/ncol(nes_matrix_top30)))
    num_fontsize_top30 <- max(5, min(8, 200/nrow(nes_matrix_top30)))
    
    tryCatch({
      png(out_png_top30, width = width_in_top30, height = height_in_top30, units = "in", res = 300)
      
      pheatmap(nes_matrix_top30,
               color = colorRampPalette(c("blue", "white", "red"))(100),
               breaks = seq(-3, 3, length.out = 101),
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               main = "Top 30 Pathways Shared Between ≥2 Models \nContaining Genes Enriched in ≥2 Pathways \np≤0.05, ** p≤0.01, *** p≤0.001" ,
               fontsize_row = row_fontsize_top30,
               fontsize_col = col_fontsize_top30,
               fontsize = 11,
               display_numbers = display_matrix_top30,
               number_format = "%s",
               number_color = "black",
               fontsize_number = num_fontsize_top30,
               border_color = NA,
               na_col = "grey95",
               legend_title = "NES")
      
      dev.off()
      message(" Created TOP 30 pathway NES heatmap")
      
    }, error = function(e) {
      message(" PNG device failed for top 30: ", e$message)
      message(" Attempting PDF alternative...")
      
      # Try PDF as fallback
      out_pdf_top30 <- file.path(cross_model_genes_dir, "pathways_cross_model_NES_heatmap_top30.pdf")
      
      tryCatch({
        pdf(out_pdf_top30, width = width_in_top30, height = height_in_top30)
        
        pheatmap(nes_matrix_top30,
                 color = colorRampPalette(c("blue", "white", "red"))(100),
                 breaks = seq(-3, 3, length.out = 101),
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 main = "Top 30 Pathways Shared Between ≥2 Models \nContaining Genes Enriched in ≥2 Pathways \np≤0.05, ** p≤0.01, *** p≤0.001",
                 fontsize_row = row_fontsize_top30,
                 fontsize_col = col_fontsize_top30,
                 fontsize = 11,
                 display_numbers = display_matrix_top30,
                 number_format = "%s",
                 number_color = "black",
                 fontsize_number = num_fontsize_top30,
                 border_color = NA,
                 na_col = "grey95",
                 legend_title = "NES")
        
        dev.off()
        message(" Created TOP 30 pathway heatmap as PDF: ", out_pdf_top30)
        
      }, error = function(e2) {
        message(" PDF also failed for top 30: ", e2$message)
      })
    })
  }
  
  #------------------------------------------------------------------------------
  # Create GENE-LEVEL heatmap for top 30 pathways
  #------------------------------------------------------------------------------
  
  message("\nCreating gene-level heatmap for top 30 pathways...")
  
  # Get the pathway keys for the top 30 pathways
  top30_pathway_keys <- recurrent_pathways[top_pathways_idx]
  top30_pathway_labels <- rownames(nes_matrix_top30)
  
  # Create a matrix to display gene names
  gene_display_matrix <- matrix("",
                                nrow = nrow(nes_matrix_top30),
                                ncol = ncol(nes_matrix_top30),
                                dimnames = list(top30_pathway_labels, colnames(nes_matrix_top30)))
  
  # Create a matrix for NES values (for coloring)
  gene_nes_matrix <- matrix(0,
                            nrow = nrow(nes_matrix_top30),
                            ncol = ncol(nes_matrix_top30),
                            dimnames = list(top30_pathway_labels, colnames(nes_matrix_top30)))
  
  # For each pathway and model, find the most representative multi-pathway gene
  for (i in seq_along(top30_pathway_keys)) {
    pathway_key <- top30_pathway_keys[i]
    pathway_label <- top30_pathway_labels[i]
    
    for (model_name in colnames(nes_matrix_top30)) {
      # Check if this pathway is enriched in this model
      if (model_name %in% names(all_pathway_genes[[pathway_key]])) {
        # Get leading edge genes for this pathway-model combination
        le_genes <- all_pathway_genes[[pathway_key]][[model_name]]
        
        # Find which of these genes are multi-pathway genes
        multi_genes_in_pathway <- intersect(le_genes, names(multi_pathway_genes))
        
        if (length(multi_genes_in_pathway) > 0) {
          # Rank genes by how many pathways they appear in
          gene_frequencies <- gene_pathway_counts[multi_genes_in_pathway]
          
          # Select the gene that appears in the most pathways (most "hub-like")
          top_gene <- names(gene_frequencies)[which.max(gene_frequencies)]
          
          # Get NES and significance for display
          nes_val <- nes_matrix_top30[pathway_label, model_name]
          padj_val <- padj_matrix[pathway_key, model_name]
          
          # Store NES for coloring
          gene_nes_matrix[pathway_label, model_name] <- nes_val
          
          # Create display text with significance markers
          if (padj_val <= 0.001) {
            gene_display_matrix[pathway_label, model_name] <- paste0(top_gene, "***")
          } else if (padj_val <= 0.01) {
            gene_display_matrix[pathway_label, model_name] <- paste0(top_gene, "**")
          } else if (padj_val <= 0.05) {
            gene_display_matrix[pathway_label, model_name] <- paste0(top_gene, "*")
          } else {
            gene_display_matrix[pathway_label, model_name] <- top_gene
          }
        }
      }
    }
  }
  
  # Create gene-level heatmap
  out_png_genes <- file.path(cross_model_genes_dir, "pathways_top30_gene_level.png")
  
  # Calculate dimensions
  width_in_genes <- max(10, ncol(gene_display_matrix) * 2 + 6)
  height_in_genes <- max(12, nrow(gene_display_matrix) * 0.5 + 3)
  
  # Cap dimensions
  width_in_genes <- min(width_in_genes, 50)
  height_in_genes <- min(height_in_genes, 60)
  
  message(" Gene-level plot size: ", round(width_in_genes, 1), " x ", round(height_in_genes, 1), " inches")
  
  # Adjust font sizes
  row_fontsize_genes <- max(6, min(10, 300/nrow(gene_display_matrix)))
  col_fontsize_genes <- max(8, min(12, 100/ncol(gene_display_matrix)))
  cell_fontsize_genes <- max(4, min(7, 200/nrow(gene_display_matrix)))
  
  tryCatch({
    png(out_png_genes, width = width_in_genes, height = height_in_genes, units = "in", res = 300)
    
    pheatmap(gene_nes_matrix,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             breaks = seq(-3, 3, length.out = 101),
             cluster_rows = TRUE,  # Keep ranking order
             cluster_cols = TRUE,
             main = "Top 30 Pathways Shared Between ≥2 Models Displaying Genes Enriched in ≥2 Pathways \nGenes are Colored by NES direction\n* p≤0.05, ** p≤0.01, *** p≤0.001",
             fontsize_row = row_fontsize_genes,
             fontsize_col = col_fontsize_genes,
             fontsize = 11,
             display_numbers = gene_display_matrix,
             number_format = "%s",
             number_color = "black",
             fontsize_number = cell_fontsize_genes,
             border_color = NA,
             na_col = "grey95",
             legend_title = "NES")
    
    dev.off()
    message(" Created gene-level heatmap for top 30 pathways")
    
  }, error = function(e) {
    message(" PNG device failed for gene-level heatmap: ", e$message)
    
    # Try PDF as fallback
    out_pdf_genes <- file.path(cross_model_genes_dir, "pathways_top30_gene_level.pdf")
    
    tryCatch({
      pdf(out_pdf_genes, width = width_in_genes, height = height_in_genes)
      
      pheatmap(gene_nes_matrix,
               color = colorRampPalette(c("blue", "white", "red"))(100),
               breaks = seq(-3, 3, length.out = 101),
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               main = "Top 30 Pathways Shared Between ≥2 Models Displaying Genes Enriched in ≥2 Pathways \nGenes are Colored by NES direction\n* p≤0.05, ** p≤0.01, *** p≤0.001" ,
               fontsize_row = row_fontsize_genes,
               fontsize_col = col_fontsize_genes,
               fontsize = 11,
               display_numbers = gene_display_matrix,
               number_format = "%s",
               number_color = "black",
               fontsize_number = cell_fontsize_genes,
               border_color = NA,
               na_col = "grey95",
               legend_title = "NES")
      
      dev.off()
      message(" Created gene-level heatmap as PDF: ", out_pdf_genes)
      
    }, error = function(e2) {
      message(" PDF also failed for gene-level heatmap: ", e2$message)
    })
  })
  
  # Save gene-pathway-model mapping
  message("\nSaving detailed gene-pathway-model mappings...")
  
  detailed_mapping <- data.frame()
  for (gene in names(multi_pathway_genes)) {
    for (item in multi_pathway_genes[[gene]]) {
      # Extract collection and pathway name
      pathway_parts <- strsplit(item$pathway, "::")[[1]]
      collection <- pathway_parts[1]
      pathway_name <- pathway_parts[2]
      
      # Get NES and padj if available
      nes_val <- NA
      padj_val <- NA
      if (item$pathway %in% names(pathway_model_nes) &&
          item$model %in% names(pathway_model_nes[[item$pathway]])) {
        nes_val <- pathway_model_nes[[item$pathway]][[item$model]]
        padj_val <- pathway_model_padj[[item$pathway]][[item$model]]
      }
      
      detailed_mapping <- rbind(detailed_mapping, data.frame(
        gene = gene,
        collection = collection,
        pathway = pathway_name,
        model = item$model,
        NES = nes_val,
        padj = padj_val,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Sort by gene frequency and padj
  detailed_mapping <- detailed_mapping %>%
    arrange(gene, padj)
  
  write.csv(detailed_mapping,
            file.path(cross_model_genes_dir, "gene_pathway_model_detailed_mapping.csv"),
            row.names = FALSE)
  
  message(" Saved detailed mapping CSV")
}

# Create summary statistics
summary_stats <- data.frame(
  Metric = c("Total pathways (before filtering)",
             "Pathways in ≥2 models (after gene filtering)",
             "Total unique genes",
             "Genes in ≥2 pathways",
             "Genes in ≥3 pathways",
             "Genes in ≥5 pathways",
             "Models analyzed"),
  Value = c(length(all_pathway_genes),
            length(recurrent_pathways),
            length(gene_pathway_map),
            length(multi_pathway_genes),
            sum(gene_pathway_counts >= 3),
            sum(gene_pathway_counts >= 5),
            length(names(full_fgsea_results)))
)

write.csv(summary_stats,
          file.path(cross_model_genes_dir, "analysis_summary.csv"),
          row.names = FALSE)

message("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
message("CROSS-MODEL GENE ANALYSIS COMPLETE!")
message(paste(rep("=", 60), collapse = ""), "\n", sep = "")
message("\nSummary:")
message(" - Total pathways collected: ", length(all_pathway_genes))
message(" - Pathways in ≥2 models: ", length(recurrent_pathways))
message(" - Genes in ≥2 pathways: ", length(multi_pathway_genes))
message(" - Output directory: ", cross_model_genes_dir)
message(paste(rep("=", 60), collapse = ""), "\n", sep = "")