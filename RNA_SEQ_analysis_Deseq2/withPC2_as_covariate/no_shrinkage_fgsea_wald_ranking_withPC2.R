start_time <- Sys.time()

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List of Bioconductor packages
bioc_packages <- c("fgsea")

# List of CRAN packages
cran_packages <- c("msigdbr", "data.table", "ggridges", "pheatmap", "dplyr",
                   "ggplot2", "tidyr", "tibble", "RColorBrewer", "ggrepel", "Cairo")


# Install Bioconductor packages if not installed
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Install CRAN packages if not installed
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
}

# Load all libraries
library(fgsea)
library(msigdbr)
library(data.table)
library(ggridges)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(ggrepel)
library(Cairo)
library(grid)
library(parallel)
library(KEGGREST)
library(biomaRt)



safe_filename <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

# Create fgsea output directory structure
fgsea_dir <- file.path(output_dir, "fgsea_enrichment")
dir.create(fgsea_dir, showWarnings = FALSE, recursive = TRUE)
fgsea_subdirs <- c("positional" , "go_bp", "go_mf", "go_cc", "wikipathways" 
                   , "hallmark" ,  "biocarta" , "reactome" , "micro_RNA" ,"pertubations" , "celltype", "immune_sig" , "kegg" , "plots", "tables")
for (subdir in fgsea_subdirs) {
  dir.create(file.path(fgsea_dir, subdir), showWarnings = FALSE, recursive = TRUE)
}
#-----------------------------------------------------------------------------
#1. # FUNCTION: Prepare Gene Sets from MSigDB
#-----------------------------------------------------------------------------
prepare_gene_sets <- function(db_species = "MM" ,
                              species = "Mus musculus", 
                              include_kegg = TRUE,
                              kegg_cache_dir = NULL,
                              kegg_force_refresh = FALSE,
                              collections = NULL) {
  
  message("\nPreparing gene sets from MSigDB...")
  
  if (is.null(collections)) {
    collections <- list(
      positional = c("M1", NULL),
      go_bp = c("M5", "GO:BP"),
      go_mf = c("M5", "GO:MF"),
      go_cc = c("M5", "GO:CC"),
      wikipathways = c("M2", "CP:WIKIPATHWAYS"),
      hallmark = c("MH", NULL),
      biocarta = c("M2", "CP:BIOCARTA"),
      reactome = c("M2", "CP:REACTOME"),
      micro_RNA = c("M3", "MIRDB") ,
      pertubations = c("M2" , "CGP"),
      celltype = c("M8" , NULL) ,
      immune_sig = c("M7" , NULL)
    )
  }
  
  gene_sets <- list()
  
  # Load MSigDB collections
  for (name in names(collections)) {
    message("  Loading ", name, " gene sets...")
    category_val <- collections[[name]][1]
    subcategory_val <- collections[[name]][2]
    
    tryCatch({
      
      msig_data <- if (is.na(subcategory_val)) {
        msigdbr(db_species = db_species,
                species = species,
                collection = category_val)
      } else {
        msigdbr(db_species = db_species,
                species = species,
                collection = category_val,
                subcollection = subcategory_val)
      }
      
      gene_sets[[name]] <- split(msig_data$gene_symbol, msig_data$gs_name)
      message("  Loaded ", length(gene_sets[[name]]), " gene sets")
    }, error = function(e) {
      message("  Failed to load ", name, ": ", e$message)
    })
  }
  
  # Add KEGG pathways from API with caching
  if (include_kegg) {
    message("  Loading KEGG pathways for mouse...")
    kegg_sets <- get_mouse_kegg_pathways(
      cache_dir = kegg_cache_dir,
      force_refresh = kegg_force_refresh,
      cache_days = 7,  # A new KEGG cache is created every week
      verbose = TRUE
    )
    
    if (length(kegg_sets) > 0) {
      gene_sets[["kegg"]] <- kegg_sets
      message("  Loaded ", length(kegg_sets), " KEGG pathways")
    }
  }
  
  gene_sets
}
#  KEGG pathway retrieval
get_mouse_kegg_pathways <- function(cache_dir = NULL, force_refresh = FALSE, 
                                    cache_days = 7, verbose = TRUE) {
  require(KEGGREST)
  require(biomaRt)
  require(digest)
  
  # Set up caching directory
  if (is.null(cache_dir)) {
    cache_dir <- file.path(getwd(), ".kegg_cache")
  }
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Create cache filenames with timestamp
  cache_file <- file.path(cache_dir, "mouse_kegg_pathways.rds")
  cache_meta <- file.path(cache_dir, "cache_metadata.rds")
  
  # Check if cache exists and is fresh
  use_cache <- FALSE
  if (!force_refresh && file.exists(cache_file) && file.exists(cache_meta)) {
    metadata <- readRDS(cache_meta)
    cache_age <- difftime(Sys.Date(), metadata$date, units = "days")
    
    if (cache_age < cache_days) {
      use_cache <- TRUE
      if (verbose) {
        message(sprintf("  Using cached KEGG data (%.1f days old)", as.numeric(cache_age)))
      }
    }
  }
  
  if (use_cache) {
    kegg_pathways <- readRDS(cache_file)
    return(kegg_pathways)
  }
  
  # If not using cache, fetch fresh data
  kegg_pathways <- list()
  
  tryCatch({
    # Record fetch timestamp
    fetch_time <- Sys.time()
    
    # Get all mouse pathways
    if (verbose) message("  Fetching fresh mouse pathway list from KEGG...")
    pathways <- keggList("pathway", "mmu")
    pathway_ids <- names(pathways)
    
    # Clean pathway names
    pathway_names <- gsub(" - Mus musculus \\(house mouse\\)", "", pathways)
    
    if (verbose) {
      message(sprintf("  Found %d mouse pathways", length(pathway_ids)))
      message("  Retrieving genes for each pathway...")
    }
    
    # Collect all unique Entrez IDs first for batch conversion
    all_entrez_ids <- character()
    pathway_entrez_map <- list()
    
    max_retries <- 3
    retry_delay <- 2  
    
    if (verbose) pb <- txtProgressBar(min = 0, max = length(pathway_ids), style = 3)
    
    for (i in seq_along(pathway_ids)) {
      pathway_id <- pathway_ids[i]
      success <- FALSE
      retry_count <- 0
      
      while (!success && retry_count < max_retries) {
        tryCatch({
          pathway_info <- keggGet(pathway_id)
          
          if (!is.null(pathway_info[[1]]$GENE)) {
            genes <- pathway_info[[1]]$GENE
            entrez_ids <- genes[seq(1, length(genes), by = 2)]
            pathway_entrez_map[[pathway_id]] <- entrez_ids
            all_entrez_ids <- unique(c(all_entrez_ids, entrez_ids))
          }
          
          success <- TRUE
          Sys.sleep(0.1) 
          
        }, error = function(e) {
          retry_count <<- retry_count + 1
          if (retry_count < max_retries) {
            Sys.sleep(retry_delay * retry_count)
          }
          NULL
        })
      }
      
      if (verbose) setTxtProgressBar(pb, i)
    }
    if (verbose) close(pb)
    
    # Sort Entrez IDs for consistency
    all_entrez_ids <- sort(unique(all_entrez_ids))
    
    # Convert to gene symbols with consistent biomaRt settings
    if (verbose) message("\n  Converting Entrez IDs to gene symbols using biomaRt...")
    
    # Use specific Ensembl archive 
    
    mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl" ,
                    host = "https://sep2025.archive.ensembl.org")
    
    # Get both MGI symbols and external gene names
    chunk_size <- 500  # Smaller chunks for stability
    entrez_to_symbol <- character()
    
    for (i in seq(1, length(all_entrez_ids), chunk_size)) {
      chunk_end <- min(i + chunk_size - 1, length(all_entrez_ids))
      chunk_ids <- all_entrez_ids[i:chunk_end]
      
      retry_count <- 0
      success <- FALSE
      
      while (!success && retry_count < max_retries) {
        tryCatch({
          chunk_conversion <- getBM(
            attributes = c("entrezgene_id", "mgi_symbol", "external_gene_name"),
            filters = "entrezgene_id",
            values = chunk_ids,
            mart = mart
          )
          
          if (nrow(chunk_conversion) > 0) {
            # Prefer MGI symbols, fall back to external names
            chunk_conversion$gene_symbol <- ifelse(
              chunk_conversion$mgi_symbol != "",
              chunk_conversion$mgi_symbol,
              chunk_conversion$external_gene_name
            )
            
            chunk_map <- setNames(
              chunk_conversion$gene_symbol,
              as.character(chunk_conversion$entrezgene_id)
            )
            entrez_to_symbol <- c(entrez_to_symbol, chunk_map)
          }
          
          success <- TRUE
          
        }, error = function(e) {
          retry_count <<- retry_count + 1
          if (retry_count < max_retries) {
            if (verbose) message(sprintf("\n  Retry %d for chunk %d", retry_count, i))
            Sys.sleep(retry_delay * retry_count)
          } else {
            warning(sprintf("Failed to convert chunk starting at %d: %s", i, e$message))
          }
          NULL
        })
      }
    }
    
    # Build pathway gene sets 
    pathway_ids_sorted <- sort(names(pathway_entrez_map))
    
    for (pathway_id in pathway_ids_sorted) {
      pathway_name <- paste0("KEGG_", 
                             toupper(gsub("[^A-Z0-9]", "_", 
                                          toupper(pathway_names[pathway_id]))))
      
      entrez_ids <- pathway_entrez_map[[pathway_id]]
      gene_symbols <- unique(na.omit(entrez_to_symbol[entrez_ids]))
      gene_symbols <- sort(gene_symbols)  # Sort for consistency
      
      if (length(gene_symbols) >= 10) {
        kegg_pathways[[pathway_name]] <- gene_symbols
      }
    }
    
    # Sort pathway names for consistent ordering
    kegg_pathways <- kegg_pathways[sort(names(kegg_pathways))]
    
    # Save to cache with metadata
    cache_metadata <- list(
      date = Sys.Date(),
      fetch_time = fetch_time,
      n_pathways = length(kegg_pathways),
      n_genes = length(unique(unlist(kegg_pathways))),
      kegg_version = packageVersion("KEGGREST"),
      biomart_version = packageVersion("biomaRt"),
      checksum = digest::digest(kegg_pathways, algo = "md5")
    )
    
    saveRDS(kegg_pathways, cache_file)
    saveRDS(cache_metadata, cache_meta)
    
    if (verbose) {
      message(sprintf("  Successfully loaded %d KEGG pathways", length(kegg_pathways)))
      message(sprintf("  Cached to: %s", cache_file))
      message(sprintf("  Cache checksum: %s", cache_metadata$checksum))
    }
    
  }, error = function(e) {
    message(sprintf("\n  Error fetching KEGG pathways: %s", e$message))
    
    # Try to load old cache as fallback
    if (file.exists(cache_file)) {
      if (verbose) message("  Attempting to load old cache as fallback...")
      kegg_pathways <- readRDS(cache_file)
      warning("Using potentially outdated cache due to fetch error")
    } else {
      return(list())
    }
  })
  
  return(kegg_pathways)
}

verify_kegg_consistency <- function(cache_dir = NULL) {
  if (is.null(cache_dir)) {
    cache_dir <- file.path(getwd(), ".kegg_cache")
  }
  
  cache_meta <- file.path(cache_dir, "cache_metadata.rds")
  
  if (file.exists(cache_meta)) {
    metadata <- readRDS(cache_meta)
    cat("\nKEGG Cache Information:\n")
    cat(sprintf("  Last updated: %s\n", metadata$date))
    cat(sprintf("  Fetch time: %s\n", metadata$fetch_time))
    cat(sprintf("  Number of pathways: %d\n", metadata$n_pathways))
    cat(sprintf("  Number of unique genes: %d\n", metadata$n_genes))
    cat(sprintf("  MD5 checksum: %s\n", metadata$checksum))
    cat(sprintf("  Cache age: %.1f days\n", 
                as.numeric(difftime(Sys.Date(), metadata$date, units = "days"))))
    return(metadata)
  } else {
    cat("No cache found\n")
    return(NULL)
  }
}
#-----------------------------------------------------------------------------
#2. # FUNCTION: RUN fgsea analysis
#-----------------------------------------------------------------------------

run_fgsea_analysis <- function(res_df, gene_sets, test_name, output_dir,
                               min_size = 15, max_size = 500, nperm = 10000,
                               score_type = "std", eps = 1e-10, seed = 42) {
  
  message("\nRunning fgsea for: ", test_name)
  message(" Mode: running with Deterministic tie-breaking and nproc = 1 & wald ranking")
  
  set.seed(seed)
  
  # Convert to data frame
  res_df <- as.data.frame(res_df)
  
  # Filter and handle duplicates
  rank_data <- res_df[!is.na(res_df$stat) & 
                        !is.na(res_df$gene_symbol) & 
                        res_df$gene_symbol != "", ]
  
  rank_data <- rank_data[order(rank_data$gene_symbol, rank_data$pvalue), ]
  rank_data <- rank_data[!duplicated(rank_data$gene_symbol), ]
  rank_data <- rank_data[order(rank_data$gene_symbol), ]
  
  # Create ranks by wald
  ranks <- setNames(rank_data$stat, rank_data$gene_symbol)
  ranks <- ranks[is.finite(ranks)]
  
  # Deterministic tie-breaking
  if(any(duplicated(ranks))) {
    warning(paste("Tied ranks detected for", test_name, "- using deterministic tie-breaking"))
    gene_order <- rank(names(ranks), ties.method = "first")
    ranks <- ranks + (gene_order * 1e-12)
  }
  
  ranks <- sort(ranks, decreasing = TRUE)
  
  message(" Genes in ranking: ", length(ranks))
  
  # Save rank data to CSV
  rank_output <- data.frame(
    gene_symbol = names(ranks),
    rank_stat = ranks,
    row.names = NULL
  )
  
  rank_csv_file <- file.path(output_dir, paste0(test_name, "_gene_ranks.csv"))
  write.csv(rank_output, rank_csv_file, row.names = FALSE)
  message(" Rank data saved to: ", rank_csv_file)
  
  fgsea_results <- list()
  
  for (collection_name in names(gene_sets)) {
    message(" Running enrichment for ", collection_name, "...")
    
    if (collection_name %in% c("go_bp", "go_mf", "go_cc")) {
      coll_min_size <- 15
      coll_max_size <- 500
    } else if (collection_name == "hallmark") {
      coll_min_size <- 15
      coll_max_size <- 500
    } else {
      coll_min_size <- min_size
      coll_max_size <- max_size
    }
    
    gs_filtered <- gene_sets[[collection_name]]
    if (is.null(gs_filtered) || !length(gs_filtered)) {
      message(" No gene sets found for ", collection_name)
      next
    }
    
    gs_sizes <- vapply(gs_filtered, length, integer(1))
    gs_filtered <- gs_filtered[gs_sizes >= coll_min_size & gs_sizes <= coll_max_size]
    
    if (!length(gs_filtered)) {
      message(" No gene sets pass size filter")
      next
    }
    
    set.seed(seed)
    fgsea_res <- tryCatch({
      fgseaMultilevel(
        pathways = gs_filtered,
        stats = ranks,
        minSize = coll_min_size,
        maxSize = coll_max_size,
        scoreType = score_type,
        eps = 0,
        nPermSimple = nperm,
        nproc = 1  
      )
    }, error = function(e) {
      message(" Multilevel failed, trying simple fgsea: ", e$message)
      fgseaSimple(
        pathways = gs_filtered,
        stats = ranks,
        minSize = coll_min_size,
        maxSize = coll_max_size,
        nperm = nperm,
        scoreType = score_type,
        nproc = 1  
      )
    })
    
    if (nrow(fgsea_res) > 0) {
      collapsed <- collapsePathways(fgsea_res, pathways = gs_filtered, stats = ranks)
      fgsea_res <- fgsea_res[pathway %in% collapsed$mainPathways]
    }
    
    
    fgsea_res <- as.data.frame(fgsea_res)
    
    fgsea_res$collection <- collection_name
    fgsea_res$neg_log10_padj <- -log10(pmax(fgsea_res$padj, eps))
    fgsea_res$abs_NES <- abs(fgsea_res$NES)
    fgsea_res$direction <- ifelse(fgsea_res$NES > 0, "Up", "Down")
    fgsea_res$significance <- ifelse(fgsea_res$padj <= 0.001, "***",
                                     ifelse(fgsea_res$padj <= 0.01, "**",
                                            ifelse(fgsea_res$padj <= 0.05, "*",
                                                   ifelse(fgsea_res$padj <= 0.1, ".", "ns"))))
    
    fgsea_res <- fgsea_res[order(fgsea_res$padj), ]
    
    fgsea_results[[collection_name]] <- fgsea_res
    
    tryCatch({
      collection_dir <- file.path(output_dir, collection_name)
      if (!dir.exists(collection_dir)) {
        dir.create(collection_dir, showWarnings = FALSE, recursive = TRUE)
      }
      
      out_csv <- file.path(collection_dir,
                           paste0(test_name, "_", collection_name, "_fgsea.csv"))
      
      # Handle leadingEdge column
      save_res <- fgsea_res
      if ("leadingEdge" %in% names(save_res)) {
        save_res$leadingEdge <- vapply(save_res$leadingEdge, function(x) {
          if (length(x) > 0) {
            paste(head(x, 100), collapse = ";")
          } else {
            ""
          }
        }, character(1))
      }
      
      write.csv(save_res, out_csv, row.names = FALSE)
      
    }, error = function(e) {
      message(" Warning: Could not save CSV for ", collection_name, ": ", e$message)
    })
    
    message(" Found ", sum(fgsea_res$padj <= 0.05, na.rm = TRUE),
            " significant pathways (padj <= 0.05)")
  }
  
  fgsea_results
}
#------------------------------------------------------------------------------
#3. # FUNCTION: fgsea Visualizations
#------------------------------------------------------------------------------
create_fgsea_plots <- function(fgsea_results, gene_sets, ranks, test_name, output_dir, top_num = 20) {
  message("Creating fgsea visualizations for: ", test_name)
  
  # Create model-specific directory
  model_plot_dir <- file.path(output_dir, "plots", test_name)
  dir.create(model_plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Combine results
  full_results <- bind_rows(fgsea_results, .id = "collection")
  
  # a) Top pathways plots
  top_pathways <- full_results %>%
    group_by(collection) %>%
    filter(padj <= 0.05) %>%
    arrange(padj) %>%
    slice_head(n = top_num) %>%
    ungroup()
  
  if (nrow(top_pathways) > 0) {
    top_pathways <- top_pathways %>%
      mutate(pathway_short = ifelse(nchar(pathway) > 50,
                                    paste0(substr(pathway, 1, 47), "..."),
                                    pathway))
    
    # Combined dot plot
    p_dot_combined <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway_short, NES))) +
      geom_point(aes(size = -log10(padj), color = NES)) +
      facet_wrap(~ collection, scales = "free_y", ncol = 2) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      scale_size_continuous(range = c(3, 10)) +
      theme_bw() +
      labs(title = paste("Top", min(top_num, nrow(top_pathways)), "Enriched Pathways:", test_name),
           x = "Normalized Enrichment Score (NES)",
           y = "Pathway",
           size = "-log10(padj)",
           color = "NES") +
      theme(axis.text.y = element_text(size = 9),
            strip.text = element_text(face = "bold", size = 11),
            title = element_text(size = 14),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 9))
    
    n_pathways <- length(unique(top_pathways$pathway))
    plot_height <- max(10, min(20, n_pathways * 0.3))
    
    ggsave(file.path(model_plot_dir, paste0("combined_top", top_num, "_pathways_dotplot.png")),
           p_dot_combined, width = 18, height = plot_height, dpi = 300)
    
    # Individual plots per collection
    for (coll in unique(top_pathways$collection)) {
      collection_data <- top_pathways %>% filter(collection == coll)
      
      if (nrow(collection_data) > 0) {
        p_dot_individual <- ggplot(collection_data,
                                   aes(x = NES, y = reorder(pathway_short, NES))) +
          geom_point(aes(size = -log10(padj), color = NES)) +
          scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
          scale_size_continuous(range = c(3, 10)) +
          theme_bw() +
          labs(title = paste("Top", nrow(collection_data), coll, "Pathways:", test_name),
               x = "Normalized Enrichment Score (NES)",
               y = "Pathway",
               size = "-log10(padj)",
               color = "NES") +
          theme(axis.text.y = element_text(size = 10),
                title = element_text(size = 14),
                legend.title = element_text(size = 10),
                legend.text = element_text(size = 9))
        
        ind_height <- max(6, min(15, nrow(collection_data) * 0.4))
        
        ggsave(file.path(model_plot_dir,
                         paste0(coll, "_top", nrow(collection_data), "_pathways_dotplot.png")),
               p_dot_individual, width = 12, height = ind_height, dpi = 300)
      }
    }
  }
  
  # b) NES distribution ridgeline plot
  if (nrow(full_results) > 0) {
    p_dist <- ggplot(full_results, aes(x = NES, y = collection, fill = collection)) +
      ggridges::geom_density_ridges(alpha = 0.7, scale = 1.5) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
      scale_fill_brewer(palette = "Set2") +
      theme_bw() +
      labs(title = paste("NES Distribution Across Collections:", test_name),
           x = "Normalized Enrichment Score",
           y = "Collection") +
      theme(legend.position = "none",
            title = element_text(size = 14),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 11))
    
    ggsave(file.path(model_plot_dir, "NES_distribution.png"),
           p_dist, width = 12, height = 8, dpi = 300)
  }
  
  # c) Volcano plots
  sig_results <- full_results %>% filter(is.finite(padj), is.finite(NES))
  
  if (nrow(sig_results) > 0) {
    # Combined volcano
    top_to_label <- sig_results %>%
      filter(padj <= 0.05) %>%
      arrange(padj) %>%
      group_by(collection) %>%
      slice_head(n = 5) %>%
      ungroup() %>%
      mutate(pathway_label = ifelse(nchar(pathway) > 30,
                                    paste0(substr(pathway, 1, 27), "..."),
                                    pathway))
    
    p_volcano_combined <- ggplot(sig_results, aes(x = NES, y = neg_log10_padj)) +
      geom_point(aes(color = significance), alpha = 0.6, size = 2) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.5) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue", alpha = 0.5) +
      geom_text_repel(
        data = top_to_label,
        aes(label = pathway_label),
        size = 2.5,
        box.padding = 0.35,
        point.padding = 0.3,
        segment.color = 'grey50',
        max.overlaps = 20,
        force = 2,
        segment.size = 0.3
      ) +
      scale_color_manual(values = c("***" = "darkred", "**" = "red",
                                    "*" = "orange", "." = "yellow", "ns" = "grey70")) +
      facet_wrap(~ collection, scales = "free") +
      theme_bw() +
      labs(title = paste("Pathway Enrichment Volcano:", test_name),
           x = "Normalized Enrichment Score",
           y = "-log10(adjusted p-value)",
           color = "Significance") +
      theme(legend.position = "bottom",
            title = element_text(size = 14),
            strip.text = element_text(face = "bold", size = 11))
    
    ggsave(file.path(model_plot_dir, "combined_pathway_volcano.png"),
           p_volcano_combined, width = 14, height = 10, dpi = 300)
    
    # Individual volcano plots
    for (coll in unique(sig_results$collection)) {
      collection_volcano <- sig_results %>% filter(collection == coll)
      
      if (nrow(collection_volcano) > 0) {
        top_to_label_ind <- collection_volcano %>%
          filter(padj <= 0.05) %>%
          arrange(padj) %>%
          slice_head(n = 10) %>%
          mutate(pathway_label = ifelse(nchar(pathway) > 40,
                                        paste0(substr(pathway, 1, 37), "..."),
                                        pathway))
        
        p_volcano_ind <- ggplot(collection_volcano, aes(x = NES, y = neg_log10_padj)) +
          geom_point(aes(color = significance), alpha = 0.6, size = 2.5) +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.5) +
          geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue", alpha = 0.5) +
          scale_color_manual(values = c("***" = "darkred", "**" = "red",
                                        "*" = "orange", "." = "yellow", "ns" = "grey70"))
        
        if (nrow(top_to_label_ind) > 0) {
          p_volcano_ind <- p_volcano_ind +
            geom_text_repel(
              data = top_to_label_ind,
              aes(label = pathway_label),
              size = 3,
              box.padding = 0.5,
              point.padding = 0.3,
              segment.color = 'grey50',
              max.overlaps = 25,
              force = 2,
              segment.size = 0.3
            )
        }
        
        p_volcano_ind <- p_volcano_ind +
          theme_bw() +
          labs(title = paste(coll, "Pathway Volcano:", test_name),
               x = "Normalized Enrichment Score",
               y = "-log10(adjusted p-value)",
               color = "Significance") +
          theme(legend.position = "right",
                title = element_text(size = 14))
        
        ggsave(file.path(model_plot_dir, paste0(coll, "_pathway_volcano.png")),
               p_volcano_ind, width = 10, height = 8, dpi = 300)
      }
    }
  }
  
  
  
  # d) Per-pathway enrichment curves
  enrichment_dir <- file.path(model_plot_dir, "enrichment_curves")
  dir.create(enrichment_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (collection_name in names(fgsea_results)) {
    if (nrow(fgsea_results[[collection_name]]) == 0) next
    
    top_paths <- fgsea_results[[collection_name]] %>%
      filter(padj <= 0.05) %>%
      arrange(padj) %>%
      head(10)
    
    if (nrow(top_paths) > 0) {
      for (i in seq_len(nrow(top_paths))) {
        pathway_name <- top_paths$pathway[i]
        gs_collection <- gene_sets[[collection_name]]
        
        if (!is.null(gs_collection) && pathway_name %in% names(gs_collection)) {
          # Create the enrichment plot
          p_enr <- plotEnrichment(gs_collection[[pathway_name]], ranks) +
            labs(title = paste0(
              substr(pathway_name, 1, 60),
              if(nchar(pathway_name) > 60) "..." else "",
              "\nNES: ", round(top_paths$NES[i], 3),
              ", padj: ", formatC(top_paths$padj[i], format = "e", digits = 2)
            )) +
            theme_minimal() +
            theme(
              title = element_text(size = 11),
              plot.background = element_rect(fill = "white", color = NA),
              panel.background = element_rect(fill = "white", color = NA)
            )
          
          
          short_filename <- sprintf("%s_pathway_%02d.png", collection_name, i)
          out_file <- file.path(enrichment_dir, short_filename)
          
          
          ggsave(
            filename = out_file,
            plot = p_enr,
            device = "png",
            type = "cairo",
            width = 10,
            height = 6,
            dpi = 300,
            units = "in",
            bg = "white"
          )
          
          message(" Saved enrichment plot: ", basename(out_file))
        }
      }
    }
  }
  
  list(full_results = full_results, top_pathways = top_pathways)
}

#------------------------------------------------------------------------------
#4 FUNCTION: Create GSEA Tables
#------------------------------------------------------------------------------
create_gsea_tables <- function(fgsea_results, gene_sets, ranks, test_name, output_dir, top_num = 10) {
  message("Creating GSEA table plots for: ", test_name)
  
  gsea_table_dir <- file.path(output_dir, "plots", test_name, "gsea_tables")
  dir.create(gsea_table_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (collection_name in names(fgsea_results)) {
    fgsea_res <- fgsea_results[[collection_name]]
    
    if (is.null(fgsea_res) || nrow(fgsea_res) == 0) next
    
   
    fgsea_res <- as.data.frame(fgsea_res)
    
    # Get significant results using base R (avoid dplyr on list columns)
    sig_idx <- which(!is.na(fgsea_res$padj) & fgsea_res$padj < 0.01)
    
    if (length(sig_idx) == 0) {
      sig_idx <- which(!is.na(fgsea_res$padj) & fgsea_res$padj <= 0.05)
    }
    
    if (length(sig_idx) > 0) {
      sig_res <- fgsea_res[sig_idx, ]
      
      # Order by absolute NES
      sig_res <- sig_res[order(-abs(sig_res$NES)), ]
      
      # Select top pathways
      n_select <- min(top_num, nrow(sig_res))
      main_pathways <- sig_res$pathway[1:n_select]
      
      output_file <- file.path(gsea_table_dir,
                               paste0(collection_name, "_gsea_table.png"))
      
      n_pathways <- length(main_pathways)
      plot_height <- max(6, min(15, n_pathways * 0.8))
      
      tryCatch({
        p <- plotGseaTable(
          pathways = gene_sets[[collection_name]][main_pathways],
          stats = ranks,
          fgseaRes = fgsea_res,
          gseaParam = 0.5,
          colwidths = c(5, 3, 0.8, 1.2, 1.2),
          render = FALSE
        )
        
        ggsave(
          filename = output_file,
          plot = p,
          device = "png",
          width = 18,
          height = plot_height,
          dpi = 300,
          units = "in",
          bg = "white"
        )
        
        message(" Created GSEA table for ", collection_name)
        
      }, error = function(e) {
        message(" Failed to create GSEA table for ", collection_name, ": ", e$message)
        
        tryCatch({
          png(output_file,
              width = 12 * 300,
              height = plot_height * 300,
              res = 300,
              bg = "white")
          
          plotGseaTable(
            pathways = gene_sets[[collection_name]][main_pathways],
            stats = ranks,
            fgseaRes = fgsea_res,
            gseaParam = 0.5,
            colwidths = c(5, 3, 0.8, 1.2, 1.2)
          )
          
          dev.off()
          message(" Created GSEA table using alternative method for ", collection_name)
          
        }, error = function(e2) {
          message(" Both methods failed for ", collection_name, ": ", e2$message)
        })
      })
    }
  }
}

#------------------------------------------------------------------------------
#5. FUNCTION: Leading Edge Analysis
#------------------------------------------------------------------------------
analyze_leading_edge <- function(fgsea_results, res_df, test_name, output_dir, top_num = 15) {
  message("Analyzing leading edge genes for: ", test_name)
  dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
  
  leading_edge_summary <- list()
  
  for (collection_name in names(fgsea_results)) {
    fgsea_res <- fgsea_results[[collection_name]]
    if (is.null(fgsea_res) || !nrow(fgsea_res)) next
    
    top_pathways <- fgsea_res %>%
      dplyr::filter(padj <= 0.05) %>%
      dplyr::arrange(padj) %>%
      head(top_num)
    if (!nrow(top_pathways)) next
    
    le_genes <- setNames(vector("list", nrow(top_pathways)), top_pathways$pathway)
    for (i in seq_len(nrow(top_pathways))) {
      le_genes[[i]] <- unlist(top_pathways$leadingEdge[i], use.names = FALSE)
    }
    names(le_genes) <- top_pathways$pathway
    
    all_le_genes <- unlist(le_genes, use.names = FALSE)
    gene_freq <- table(all_le_genes)
    core_genes <- names(gene_freq)[gene_freq >= 2]
    
    if (length(core_genes) > 0) {
      core_gene_data <- res_df %>%
        dplyr::filter(gene_symbol %in% core_genes) %>%
        dplyr::select(gene_symbol, log2FoldChange, padj) %>%
        dplyr::arrange(padj)
      
      csv_file <- file.path(output_dir, "tables",
                            paste0(test_name, "_", collection_name, "_core_genes.csv"))
      write.csv(core_gene_data, csv_file, row.names = FALSE)
      
      if (length(le_genes) > 1 && length(core_genes) > 1) {
        cols <- core_genes[1:min(50, length(core_genes))]
        
        pathway_names <- names(le_genes)
        pathway_names_temp <- ifelse(nchar(pathway_names) > 40,
                                     paste0(substr(pathway_names, 1, 37), "..."),
                                     pathway_names)
        pathway_names_short <- make.unique(pathway_names_temp, sep = "_")
        
        pathway_gene_matrix <- matrix(
          0L, nrow = length(le_genes), ncol = length(cols),
          dimnames = list(pathway_names_short, cols)
        )
        
        for (i in seq_along(le_genes)) {
          present <- intersect(le_genes[[i]], cols)
          # Use the unique short name to assign the value
          if (length(present)) pathway_gene_matrix[pathway_names_short[i], present] <- 1L
        }
        
        png_file <- file.path(output_dir, "plots", test_name,
                              paste0(collection_name, "_leading_edge_heatmap.png"))
        
        w_px <- min(7500, max(3000, ncol(pathway_gene_matrix) * 120 + 900))
        h_px <- min(6000, max(2400, nrow(pathway_gene_matrix) * 120 + 600))
        
        tryCatch({
          png(png_file, width = w_px, height = h_px, res = 300)
          
          pheatmap(
            pathway_gene_matrix,
            color = c("#FFE4B5", "#4B0082"),
            breaks = c(-0.5, 0.5, 1.5),
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            main = paste("Leading Edge Genes:", collection_name, "-", test_name),
            fontsize_row = max(6, min(10, 100/nrow(pathway_gene_matrix))),
            fontsize_col = max(6, min(10, 200/ncol(pathway_gene_matrix))),
            fontsize = 12,
            border_color = "grey80",
            legend_breaks = c(0, 1),
            legend_labels = c("Absent", "Present")
          )
          
          dev.off()
        }, error = function(e) {
          message(" PNG device failed, using alternative method: ", e$message)
          p <- pheatmap(
            pathway_gene_matrix,
            color = c("#FFE4B5", "#4B0082"),
            breaks = c(-0.5, 0.5, 1.5),
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            main = paste("Leading Edge Genes:", collection_name, "-", test_name),
            fontsize_row = max(6, min(10, 100/nrow(pathway_gene_matrix))),
            fontsize_col = max(6, min(10, 200/ncol(pathway_gene_matrix))),
            fontsize = 12,
            border_color = "grey80",
            legend_breaks = c(0, 1),
            legend_labels = c("Absent", "Present"),
            silent = TRUE
          )
          w_in <- min(25, max(10, ncol(pathway_gene_matrix) * 0.4 + 3))
          h_in <- min(20, max(8, nrow(pathway_gene_matrix) * 0.4 + 2))
          ggsave(png_file, p, width = w_in, height = h_in, dpi = 300)
        })
      }
    }
    
    leading_edge_summary[[collection_name]] <- list(
      top_pathways = top_pathways,
      core_genes = core_genes,
      gene_frequencies = gene_freq
    )
  }
  
  leading_edge_summary
}
#------------------------------------------------------------------------------
# RUN FGSEA FOR SELECTED MODELS ONLY
#------------------------------------------------------------------------------
message("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
message("STARTING GENE SET ENRICHMENT ANALYSIS")
message(paste(rep("=", 60), collapse = ""), "\n", sep = "")

# Define which models to analyze
selected_models <- c(
  "treatment_main",
  "sex_main",
  "strain_main" ,
  "treatment_effect_in_male_C57BL6",
  "treatment_effect_in_male_DBA2J",
  "treatment_effect_in_female_DBA2J",
  "treatment_effect_in_female_C57BL6"
)

# Verify all selected models exist in full_results
available_models <- intersect(selected_models, names(full_results))
missing_models <- setdiff(selected_models, names(full_results))

if (length(missing_models) > 0) {
  warning("The following models were not found in full_results: ", 
          paste(missing_models, collapse = ", "))
}

message("Processing ", length(available_models), " selected models:")
message(paste(available_models, collapse = "\n"))
message("\n")

gene_sets <- prepare_gene_sets(species = "Mus musculus", include_kegg = TRUE)
full_fgsea_results <- list()
all_leading_edge <- list()

# Loop through ONLY the selected models
for (model_name in available_models) {
  message("\n", paste(rep("-", 40), collapse = ""), sep = "")
  message("Processing model: ", model_name)
  message(paste(rep("-", 40), collapse = ""), sep = "")
  
  
  fgsea_res <- run_fgsea_analysis(
    res_df = full_results[[model_name]]$results,
    gene_sets = gene_sets,
    test_name = model_name,
    output_dir = fgsea_dir,
    min_size = 15,
    max_size = 500,
    nperm = 10000
  )
  full_fgsea_results[[model_name]] <- fgsea_res
  
  # Count significant pathways
  total_significant <- 0
  for (collection in names(fgsea_res)) {
    if (!is.null(fgsea_res[[collection]])) {
      total_significant <- total_significant + sum(fgsea_res[[collection]]$padj <= 0.05, na.rm = TRUE)
    }
  }
  
  message(" Total significant pathways for ", model_name, ": ", total_significant)
  
  # Create plots if enough significant pathways
  if (total_significant >= 5) {
    message(" Creating visualizations (", total_significant, " significant pathways found)...")
    # Convert to data frame
    res_df <- as.data.frame(res_df)
    
    # Filter and handle duplicates
    rank_data <- res_df[!is.na(res_df$stat) & 
                          !is.na(res_df$gene_symbol) & 
                          res_df$gene_symbol != "", ]
    
    rank_data <- rank_data[order(rank_data$gene_symbol, rank_data$pvalue), ]
    rank_data <- rank_data[!duplicated(rank_data$gene_symbol), ]
    rank_data <- rank_data[order(rank_data$gene_symbol), ]
    
    # Create ranks by wald
    ranks <- setNames(rank_data$stat, rank_data$gene_symbol)
    ranks <- ranks[is.finite(ranks)]
    
    # Deterministic tie-breaking
    if(any(duplicated(ranks))) {
      warning(paste("Tied ranks detected for", test_name, "- using deterministic tie-breaking"))
      gene_order <- rank(names(ranks), ties.method = "first")
      ranks <- ranks + (gene_order * 1e-12)
    }
    
    # Sort by rank value for fgsea
    ranks <- sort(ranks, decreasing = TRUE)
    
    # Create plots (PNG at 300 DPI)
    create_fgsea_plots(
      fgsea_results = fgsea_res,
      gene_sets = gene_sets,
      ranks = ranks,
      test_name = model_name,
      output_dir = fgsea_dir,
      top_num = 20
    )
    
    # Create GSEA table plots
    create_gsea_tables(
      fgsea_results = fgsea_res,
      gene_sets = gene_sets,
      ranks = ranks,
      test_name = model_name,
      output_dir = fgsea_dir,
      top_num = 15
    )
    
    # Leading edge analysis
    le_analysis <- analyze_leading_edge(
      fgsea_results = fgsea_res,
      res_df = full_results[[model_name]]$results,
      test_name = model_name,
      output_dir = fgsea_dir
    )
    all_leading_edge[[model_name]] <- le_analysis
    
  } else {
    message(" Skipping visualizations (only ", total_significant, " significant pathways, need greater than or equal to 5)")
    
    # Save leading edge CSV data without plots
    le_analysis <- list()
    for (collection_name in names(fgsea_res)) {
      fgsea_res_col <- fgsea_res[[collection_name]]
      if (!is.null(fgsea_res_col) && nrow(fgsea_res_col) > 0) {
        top_pathways <- fgsea_res_col %>%
          dplyr::filter(padj <= 0.05) %>%
          dplyr::arrange(padj) %>%
          head(10)
        
        if (nrow(top_pathways) > 0) {
          le_genes <- setNames(vector("list", nrow(top_pathways)), top_pathways$pathway)
          for (i in seq_len(nrow(top_pathways))) {
            le_genes[[i]] <- unlist(top_pathways$leadingEdge[i], use.names = FALSE)
          }
          
          all_le_genes <- unlist(le_genes, use.names = FALSE)
          gene_freq <- table(all_le_genes)
          core_genes <- names(gene_freq)[gene_freq >= 2]
          
          if (length(core_genes) > 0) {
            core_gene_data <- full_results[[model_name]]$results %>%
              dplyr::filter(gene_symbol %in% core_genes) %>%
              dplyr::select(gene_symbol, log2FoldChange, padj) %>%
              dplyr::arrange(padj)
            
            csv_file <- file.path(fgsea_dir, "tables",
                                  paste0(model_name, "_", collection_name, "_core_genes.csv"))
            write.csv(core_gene_data, csv_file, row.names = FALSE)
          }
        }
      }
    }
    all_leading_edge[[model_name]] <- le_analysis
  }
}

#------------------------------------------------------------------------------
# CROSS-MODEL PATHWAY COMPARISON
#------------------------------------------------------------------------------
message("\nCreating cross-model pathway comparisons...")
combined_fgsea <- list()
for (collection_name in names(gene_sets)) {
  collection_results <- list()
  for (model_name in names(full_fgsea_results)) {
    if (collection_name %in% names(full_fgsea_results[[model_name]])) {
      model_res <- full_fgsea_results[[model_name]][[collection_name]] %>%
        dplyr::mutate(model = model_name) %>%
        dplyr::select(pathway, NES, padj, model)
      collection_results[[model_name]] <- model_res
    }
  }
  if (length(collection_results)) {
    combined_fgsea[[collection_name]] <- bind_rows(collection_results)
  }
}
for (collection_name in names(combined_fgsea)) {
  combined_data <- combined_fgsea[[collection_name]]
  
  sig_pathways <- combined_data %>%
    filter(padj <= 0.05) %>%
    group_by(pathway) %>%
    summarise(min_padj = min(padj), .groups = "drop") %>%
    arrange(min_padj) %>%
    head(30) %>%
    pull(pathway)
  
  if (length(sig_pathways) > 2) {
    pathway_df <- data.frame(pathway = sig_pathways) %>%
      mutate(pathway_short = ifelse(nchar(pathway) > 50,
                                    paste0(substr(pathway, 1, 47), "..."),
                                    pathway))
    pathway_df$label <- make.unique(pathway_df$pathway_short, sep = "_")
    
    combined_data_plot <- combined_data %>%
      dplyr::filter(pathway %in% sig_pathways) %>%
      mutate(
        NES = as.numeric(NES),
        padj = as.numeric(padj),
        NES_display = case_when(
          padj <= 0.05 & abs(NES) >= 1 ~ paste0(round(NES, 2), "*"),
          TRUE ~ as.character(round(NES, 2))
        )
      ) %>%
      dplyr::filter(is.finite(NES)) %>%
      left_join(pathway_df %>% dplyr::select(pathway, label), by = "pathway")
    
    wide_color <- combined_data_plot %>%
      dplyr::select(label, model, NES) %>%
      tidyr::pivot_wider(
        names_from = model,
        values_from = NES,
        values_fill = 0
      ) %>%
      as.data.frame()
    
    wide_display <- combined_data_plot %>%
      dplyr::select(label, model, NES_display) %>%
      tidyr::pivot_wider(
        names_from = model,
        values_from = NES_display,
        values_fill = "0"
      ) %>%
      as.data.frame()
    
    row.names(wide_color) <- wide_color$label
    row.names(wide_display) <- wide_display$label
    
    nes_matrix <- as.matrix(wide_color[, setdiff(names(wide_color), "label"), drop = FALSE])
    display_matrix <- as.matrix(wide_display[, setdiff(names(wide_display), "label"), drop = FALSE])
    
    if (nrow(nes_matrix) > 1 && ncol(nes_matrix) > 1) {
      out_png <- file.path(fgsea_dir, "plots",
                           paste0("cross_model_", collection_name, "_NES_heatmap.png"))
      
      w_px <- max(3000, ncol(nes_matrix) * 600 + 1200)
      h_px <- max(3000, nrow(nes_matrix) * 105 + 600)
      
      tryCatch({
        png(out_png, width = w_px, height = h_px, res = 300)
        
        pheatmap(nes_matrix,
                 color = colorRampPalette(c("blue", "white", "red"))(100),
                 breaks = seq(-3, 3, length.out = 101),
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 main = paste("Cross-Model NES Comparison for :", collection_name , "          padj <= 0.05 is denoted as * " ),
                 fontsize_row = max(7, min(10, 300/nrow(nes_matrix))),
                 fontsize_col = 11,
                 fontsize = 13,
                 display_numbers = display_matrix,
                 number_format = "%s",
                 number_color = "black",
                 fontsize_number = max(6, min(9, 200/nrow(nes_matrix))),
                 border_color = NA,
                 na_col = "grey90")
        
        dev.off()
      }, error = function(e) {
        message(" PNG device failed for cross-model heatmap, using alternative: ", e$message)
        p <- pheatmap(nes_matrix,
                      color = colorRampPalette(c("blue", "white", "red"))(100),
                      breaks = seq(-3, 3, length.out = 101),
                      cluster_rows = TRUE,
                      cluster_cols = TRUE,
                      main = paste("Cross-Model NES Comparison:", collection_name),
                      fontsize_row = max(7, min(10, 300/nrow(nes_matrix))),
                      fontsize_col = 11,
                      fontsize = 13,
                      display_numbers = display_matrix,
                      number_format = "%s",
                      number_color = "black",
                      fontsize_number = max(6, min(9, 200/nrow(nes_matrix))),
                      border_color = NA,
                      na_col = "grey90",
                      silent = TRUE)
        
        w_in <- max(10, ncol(nes_matrix) * 2 + 4)
        h_in <- max(10, nrow(nes_matrix) * 0.35 + 2)
        ggsave(out_png, p, width = w_in, height = h_in, dpi = 300)
      })
    }
  }
}
#------------------------------------------------------------------------------
# GLOBAL FDR CORRECTION ACROSS ALL COLLECTIONS
#------------------------------------------------------------------------------
message("\nApplying global FDR correction across all collections...")
for (model_name in names(full_fgsea_results)) {
  all_pvals <- c()
  pathway_info <- list()
  
  for (collection_name in names(full_fgsea_results[[model_name]])) {
    res <- full_fgsea_results[[model_name]][[collection_name]]
    if (nrow(res) > 0) {
      all_pvals <- c(all_pvals, res$pval)
      pathway_info[[length(pathway_info) + 1]] <- data.frame(
        collection = collection_name,
        pathway = res$pathway,
        original_padj = res$padj
      )
    }
  }
  
  if (length(all_pvals) > 0) {
    global_padj <- p.adjust(all_pvals, method = "BH")
    
    idx <- 1
    for (collection_name in names(full_fgsea_results[[model_name]])) {
      res <- full_fgsea_results[[model_name]][[collection_name]]
      if (nrow(res) > 0) {
        n_pathways <- nrow(res)
        res$global_padj <- global_padj[idx:(idx + n_pathways - 1)]
        full_fgsea_results[[model_name]][[collection_name]] <- res
        idx <- idx + n_pathways
      }
    }
  }
}

verify_kegg_consistency()
#------------------------------------------------------------------------------
# SUMMARY STATISTICS
#------------------------------------------------------------------------------
fgsea_summary <- data.frame()
for (model_name in names(full_fgsea_results)) {
  for (collection_name in names(full_fgsea_results[[model_name]])) {
    res <- full_fgsea_results[[model_name]][[collection_name]]
    summary_row <- data.frame(
      Model = model_name,
      Collection = collection_name,
      Total_Pathways = nrow(res),
      Sig_0.05 = sum(res$padj <= 0.05, na.rm = TRUE),
      Sig_0.01 = sum(res$padj <= 0.01, na.rm = TRUE),
      global_padj_0.05 = sum(res$global_padj <= 0.05, na.rm = TRUE) ,
      global_padj_0.01 = sum(res$global_padj <= 0.01, na.rm = TRUE) ,
      Upregulated = sum(res$padj <= 0.05 & res$NES > 0, na.rm = TRUE),
      Downregulated = sum(res$padj <= 0.05 & res$NES < 0, na.rm = TRUE),
      Mean_abs_NES = mean(abs(res$NES[res$padj <= 0.05]), na.rm = TRUE),
      Max_NES = ifelse(any(res$padj <= 0.05, na.rm = TRUE),
                       max(res$NES[res$padj <= 0.05], na.rm = TRUE), NA),
      Min_NES = ifelse(any(res$padj <= 0.05, na.rm = TRUE),
                       min(res$NES[res$padj <= 0.05], na.rm = TRUE), NA)
    )
    fgsea_summary <- rbind(fgsea_summary, summary_row)
  }
}
write.csv(fgsea_summary,
          file.path(fgsea_dir, "tables", "fgsea_summary_all_models.csv"),
          row.names = FALSE)
cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("FGSEA ANALYSIS COMPLETE!\n")
cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("\nFGSEA Summary Statistics:\n")
cat(paste(rep("-", 60), collapse = ""), "\n", sep = "")
print(fgsea_summary)
cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")
sessionInfo()

end_time <- Sys.time()
total_time <- end_time - start_time
message(" This script ran for :" , total_time )