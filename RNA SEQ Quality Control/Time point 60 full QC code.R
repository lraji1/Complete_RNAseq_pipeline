##-----------------------------------------------------------------------
# RNA-seq QC Pipeline (Timepoint 60 Only, Outlier Excluded) + PC2 Correlation Analysis
# Outputs: QC plots, correlation analysis, PC2 associations with metadata
# REVISED: Individual PCA plots + Outlier density comparison
##-----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(pheatmap)
  library(patchwork)
  library(scales)
  library(ggrepel)
  library(matrixStats)
})

##-----------------------------------------------------------------------
# 1. LOAD DATA
##-----------------------------------------------------------------------

message("=== Loading Data ===")

counts <- read.csv("count.csv", row.names = 1, check.names = FALSE)
met <- read.csv("metv2.csv", row.names = 1, check.names = FALSE)

# Keep only numeric count columns; make integer matrix
count_cols <- sapply(counts, is.numeric)
counts <- counts[, count_cols, drop = FALSE]
counts <- round(as.matrix(counts))
mode(counts) <- "integer"

# Match samples between counts and met
common_samples <- intersect(colnames(counts), rownames(met))
counts <- counts[, common_samples, drop = FALSE]
met <- met[common_samples, , drop = FALSE]

##-----------------------------------------------------------------------
# 1.5 FILTER FOR TIMEPOINT 60 (WITH OUTLIER FIRST)
##-----------------------------------------------------------------------

message("\n=== Filtering for Timepoint 60 ===")

if (!"Time.Point" %in% colnames(met)) {
  stop("Time.Point column not found in metadata!")
}

outlier_sample <- "69.B6.F.60.THC"

# Get ALL timepoint 60 samples (including outlier)
tp60_samples_with_outlier <- rownames(met)[met$Time.Point == "60"]

if (length(tp60_samples_with_outlier) == 0) {
  stop("No samples found at timepoint 60!")
}

message(sprintf("Total samples at timepoint 60 (with outlier): %d", length(tp60_samples_with_outlier)))

# Create dataset WITH outlier for initial density plot
counts_with_outlier <- counts[, tp60_samples_with_outlier, drop = FALSE]
met_with_outlier <- met[tp60_samples_with_outlier, , drop = FALSE]

##-----------------------------------------------------------------------
# 1.6 SET FACTOR LEVELS (FOR WITH-OUTLIER DATA)
##-----------------------------------------------------------------------

# Factors
if ("RNA.Batch" %in% colnames(met_with_outlier)) {
  met_with_outlier$RNA.Batch <- factor(met_with_outlier$RNA.Batch)
}
if ("Library.Batch" %in% colnames(met_with_outlier)) {
  met_with_outlier$Library.Batch <- factor(met_with_outlier$Library.Batch)
}
if ("Treatment" %in% colnames(met_with_outlier)) {
  met_with_outlier$Treatment <- factor(met_with_outlier$Treatment, levels = c("None", "vehicle", "THC"))
}
if ("Sex" %in% colnames(met_with_outlier)) {
  met_with_outlier$Sex <- factor(met_with_outlier$Sex, levels = c("female", "male"))
}
if ("Strain" %in% colnames(met_with_outlier)) {
  met_with_outlier$Strain <- factor(met_with_outlier$Strain)
}

# Determine primary batch variable
if ("RNA.Batch" %in% colnames(met_with_outlier)) {
  met_with_outlier$batch <- met_with_outlier$RNA.Batch
} else if ("Library.Batch" %in% colnames(met_with_outlier)) {
  met_with_outlier$batch <- met_with_outlier$Library.Batch
} else {
  stop("No batch column found.")
}

n_batches <- length(unique(met_with_outlier$batch))
batch_colors <- colorRampPalette(brewer.pal(9, "Set1"))(n_batches)
names(batch_colors) <- levels(met_with_outlier$batch)

##-----------------------------------------------------------------------
# 2.0 EXPRESSION DENSITY PLOT WITH OUTLIER
##-----------------------------------------------------------------------

message("\n=== Generating Expression Density Plot WITH Outlier ===")

qc_output_dir <- "QC_Results_Timepoint60_NoOutlier"
dir.create(qc_output_dir, showWarnings = FALSE, recursive = TRUE)

gene_means_with_outlier <- rowMeans(counts_with_outlier)
top_genes_with_outlier <- names(sort(gene_means_with_outlier, decreasing = TRUE)[1:min(23040, length(gene_means_with_outlier))])

log2_counts_with_outlier <- log2(counts_with_outlier[top_genes_with_outlier, , drop = FALSE] + 0.5)

set.seed(42)
genes_for_plot_with_outlier <- sample(top_genes_with_outlier, min(23040, length(top_genes_with_outlier)))

log2_long_with_outlier <- as.data.frame(log2_counts_with_outlier[genes_for_plot_with_outlier, , drop = FALSE]) %>%
  mutate(gene = rownames(.)) %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "log2_expression")

log2_long_with_outlier <- merge(
  log2_long_with_outlier,
  data.frame(sample = rownames(met_with_outlier), met_with_outlier),
  by = "sample"
)

# Calculate density for outlier sample to find peak location for labeling
outlier_data <- log2_long_with_outlier %>%
  filter(sample == outlier_sample)

if (nrow(outlier_data) > 0) {
  # Calculate density to find peak
  dens <- density(outlier_data$log2_expression)
  peak_x <- dens$x[which.max(dens$y)]
  peak_y <- max(dens$y)
} else {
  peak_x <- 0
  peak_y <- 0
}

# Density plot WITH outlier
p_density_with_outlier <- ggplot(log2_long_with_outlier, aes(x = log2_expression, fill = batch)) +
  geom_density(alpha = 0.5, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = batch_colors) +
  labs(
    title = "Gene Expression Distribution by Batch (Timepoint 60, WITH Outlier)",
    subtitle = "log2(counts + 0.5), sampled genes",
    x = "log2(counts + 0.5)",
    y = "Density",
    fill = "Batch"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2))

# Add annotation for outlier sample
if (nrow(outlier_data) > 0) {
  p_density_with_outlier <- p_density_with_outlier +
    annotate("text", x = peak_x, y = peak_y * 1.1, 
             label = outlier_sample,
             color = "red", fontface = "bold", size = 3.5,
             hjust = 0.5) +
    annotate("segment", x = peak_x, xend = peak_x,
             y = peak_y, yend = peak_y * 1.08,
             arrow = arrow(length = unit(0.2, "cm")),
             color = "red", linewidth = 0.8)
}

ggsave(file.path(qc_output_dir, "01a_expression_density_WITH_outlier.png"),
       p_density_with_outlier, width = 10, height = 7, dpi = 300)

message(sprintf("Expression density WITH outlier saved (outlier: %s)", outlier_sample))

##-----------------------------------------------------------------------
# 1.7 NOW FILTER TO EXCLUDE OUTLIER
##-----------------------------------------------------------------------

message("\n=== Filtering to EXCLUDE Outlier ===")

# Get timepoint 60 samples, excluding the outlier
tp60_samples <- rownames(met)[met$Time.Point == "60" & met$Name != outlier_sample]

if (length(tp60_samples) > 0) {
  message(sprintf("Found %d samples at timepoint 60 (after excluding outlier)", length(tp60_samples)))
  message(sprintf("Excluded outlier sample: %s", outlier_sample))
}

message(sprintf("Total samples before filtering: %d", ncol(counts)))
message(sprintf("Samples at timepoint 60: %d", length(tp60_samples)))

# Filter both counts and metadata to timepoint 60 only (no outlier)
counts <- counts[, tp60_samples, drop = FALSE]
met <- met[tp60_samples, , drop = FALSE]

message(sprintf("After timepoint 60 filtering: %d samples", ncol(counts)))

##-----------------------------------------------------------------------
# 1.8 SET FACTOR LEVELS (FOR NO-OUTLIER DATA)
##-----------------------------------------------------------------------

# Factors
if ("RNA.Batch" %in% colnames(met)) {
  met$RNA.Batch <- factor(met$RNA.Batch)
}
if ("Library.Batch" %in% colnames(met)) {
  met$Library.Batch <- factor(met$Library.Batch)
}
if ("Treatment" %in% colnames(met)) {
  met$Treatment <- factor(met$Treatment, levels = c("None", "vehicle", "THC"))
}
if ("Sex" %in% colnames(met)) {
  met$Sex <- factor(met$Sex, levels = c("female", "male"))
}
if ("Strain" %in% colnames(met)) {
  met$Strain <- factor(met$Strain)
}

# Add group1-group7 as factors (automatically detect levels)
for (i in 1:7) {
  group_col <- paste0("Group", i)
  if (group_col %in% colnames(met)) {
    met[[group_col]] <- factor(met[[group_col]])
    message(sprintf("%s levels: %s", group_col, paste(levels(met[[group_col]]), collapse = ", ")))
  }
}

# Create separate folder for individual PCA plots
pca_by_variable_dir <- file.path(qc_output_dir, "PCA_by_Variable")
dir.create(pca_by_variable_dir, showWarnings = FALSE, recursive = TRUE)

# Create separate folder for group PCA plots
group_pca_dir <- file.path(qc_output_dir, "Group_PCA_Plots")
dir.create(group_pca_dir, showWarnings = FALSE, recursive = TRUE)

message(sprintf("Samples: %d", ncol(counts)))
message(sprintf("Genes: %d", nrow(counts)))

# Determine primary batch variable for coloring
if ("RNA.Batch" %in% colnames(met)) {
  met$batch <- met$RNA.Batch
} else if ("Library.Batch" %in% colnames(met)) {
  met$batch <- met$Library.Batch
} else {
  stop("No batch column found.")
}

message(sprintf("RNA Batches: %s", paste(levels(met$RNA.Batch), collapse = ", ")))
if ("Library.Batch" %in% colnames(met)) {
  message(sprintf("Library Batches: %s", paste(levels(met$Library.Batch), collapse = ", ")))
}

##-----------------------------------------------------------------------
# 2. EXPRESSION DENSITY PLOTS (WITHOUT OUTLIER)
##-----------------------------------------------------------------------

message("\n=== Generating Expression Density Plots (WITHOUT Outlier) ===")

gene_means <- rowMeans(counts)
top_genes <- names(sort(gene_means, decreasing = TRUE)[1:min(23040, length(gene_means))])

log2_counts <- log2(counts[top_genes, , drop = FALSE] + 0.5)

set.seed(42)
genes_for_plot <- sample(top_genes, min(23040, length(top_genes)))

log2_long <- as.data.frame(log2_counts[genes_for_plot, , drop = FALSE]) %>%
  mutate(gene = rownames(.)) %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "log2_expression")

log2_long <- merge(
  log2_long,
  data.frame(sample = rownames(met), met),
  by = "sample"
)

# Get number of samples and batches for dynamic settings
n_samples <- length(unique(log2_long$sample))
n_batches <- length(unique(met$batch))

# Create a color palette for batches
batch_colors <- colorRampPalette(brewer.pal(9, "Set1"))(n_batches)
names(batch_colors) <- levels(met$batch)

# 2.1 Density by batch (aggregated) - WITHOUT outlier
p_density_batch <- ggplot(log2_long, aes(x = log2_expression, fill = batch)) +
  geom_density(alpha = 0.5, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = batch_colors) +
  labs(
    title = "Gene Expression Distribution by Batch (Timepoint 60, No Outlier)",
    subtitle = "log2(counts + 0.5), sampled genes",
    x = "log2(counts + 0.5)",
    y = "Density",
    fill = "Batch"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2))

# 2.2 Density by treatment (since all samples are same timepoint)
if ("Treatment" %in% colnames(log2_long)) {
  p_density_treatment <- ggplot(log2_long, aes(x = log2_expression, fill = Treatment)) +
    geom_density(alpha = 0.5, color = "black", linewidth = 0.3) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "Gene Expression Distribution by Treatment (Timepoint 60, No Outlier)",
      x = "log2(counts + 0.5)",
      y = "Density",
      fill = "Treatment"
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom")
} else {
  p_density_treatment <- ggplot() + theme_void() +
    labs(title = "Treatment not found in met")
}

combined_density <- (p_density_batch | p_density_treatment) +
  plot_annotation(
    title = "Expression Distribution Quality Control (Timepoint 60, No Outlier)",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

ggsave(file.path(qc_output_dir, "01b_expression_density_NO_outlier.png"),
       combined_density, width = 14, height = 7, dpi = 300)

message("Expression density plots saved (WITHOUT outlier)")

##-----------------------------------------------------------------------
# 2.5 ALL SAMPLES DENSITY PLOT WITH BATCH INFO IN LEGEND
##-----------------------------------------------------------------------

message("\n=== Generating All-Samples Density Plot (with batch in legend) ===")

# Create sample labels that include batch info: "SampleName [B#]"
sample_labels <- paste0(rownames(met), " [B", met$batch, "]")
names(sample_labels) <- rownames(met)

# Create colors for each sample based on their batch
sample_colors_labeled <- batch_colors[as.character(met$batch)]
names(sample_colors_labeled) <- sample_labels

# Add sample_label to log2_long for plotting
log2_long$sample_label <- sample_labels[log2_long$sample]

# Sort sample labels by batch so legend is organized
sample_order <- order(met$batch, rownames(met))
sample_labels_ordered <- sample_labels[sample_order]
log2_long$sample_label <- factor(log2_long$sample_label, levels = sample_labels_ordered)

# Calculate number of legend columns based on sample count
legend_ncol <- ceiling(n_samples / 25)

# Create the all-samples density plot
p_density_all_samples <- ggplot(log2_long, 
                                aes(x = log2_expression, 
                                    color = sample_label, 
                                    group = sample)) +
  geom_density(linewidth = 0.5, alpha = 0.8) +
  scale_color_manual(values = sample_colors_labeled, 
                     name = "Sample [Batch]",
                     guide = guide_legend(ncol = legend_ncol, 
                                          byrow = FALSE,
                                          keywidth = 1.2,
                                          keyheight = 0.4,
                                          override.aes = list(linewidth = 1.5))) +
  labs(
    title = "Gene Expression Distribution: All Samples Individual Curves (TP60, No Outlier)",
    subtitle = paste0("Each line = one sample (", n_samples, " total), colored by batch (", 
                      n_batches, " batches)"),
    x = "log2(counts + 0.5)",
    y = "Density"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom",
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 10, face = "bold"),
    legend.key.size = unit(0.25, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.box.spacing = unit(0.1, "cm")
  )

# Calculate dynamic figure height based on number of samples
legend_rows <- ceiling(n_samples / legend_ncol)
fig_height <- 8 + (legend_rows * 0.18)

ggsave(file.path(qc_output_dir, "01c_expression_density_all_samples.png"),
       p_density_all_samples, 
       width = 18, 
       height = min(fig_height, 30),
       dpi = 300,
       limitsize = FALSE)

message(sprintf("All-samples density plot saved (%d samples, %d batches)", n_samples, n_batches))

##-----------------------------------------------------------------------
# 2.6 BATCH COLOR KEY
##-----------------------------------------------------------------------

message("\n=== Creating Batch Color Key ===")

batch_key_df <- data.frame(
  batch = factor(levels(met$batch), levels = levels(met$batch)),
  x = 1,
  y = seq_along(levels(met$batch))
)

samples_per_batch <- table(met$batch)
batch_key_df$label <- paste0(batch_key_df$batch, " (n=", samples_per_batch[as.character(batch_key_df$batch)], ")")

p_batch_key <- ggplot(batch_key_df, aes(x = x, y = reorder(batch, -y), fill = batch)) +
  geom_tile(width = 0.6, height = 0.8, color = "black", linewidth = 0.5) +
  geom_text(aes(x = x + 0.5, label = label), hjust = 0, size = 4, fontface = "bold") +
  scale_fill_manual(values = batch_colors) +
  labs(title = "Batch Color Key (Timepoint 60, No Outlier)",
       subtitle = "Reference for density plots") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 10, hjust = 0),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  xlim(0.5, 4)

ggsave(file.path(qc_output_dir, "00_batch_color_key.png"),
       p_batch_key, width = 5, height = max(3, n_batches * 0.6), dpi = 300)

message("Batch color key saved")

##-----------------------------------------------------------------------
# 2.7 SAMPLE-BATCH MAPPING TABLE
##-----------------------------------------------------------------------

batch_sample_table <- data.frame(
  sample = rownames(met),
  RNA.Batch = if("RNA.Batch" %in% colnames(met)) met$RNA.Batch else NA,
  Library.Batch = if("Library.Batch" %in% colnames(met)) met$Library.Batch else NA,
  sample_label = sample_labels
) %>%
  arrange(RNA.Batch, sample)

write.csv(batch_sample_table, 
          file.path(qc_output_dir, "sample_batch_mapping.csv"),
          row.names = FALSE)

message("Sample-batch mapping saved to CSV")

##-----------------------------------------------------------------------
# 2.8 OUTLIER DETECTION (≥2 SD from mean)
##-----------------------------------------------------------------------

message("\n=== Detecting Expression Outliers (≥2 SD) ===")

# Calculate mean expression per sample
sample_mean_expression <- log2_long %>%
  group_by(sample) %>%
  summarise(
    mean_log2_expr = mean(log2_expression, na.rm = TRUE),
    median_log2_expr = median(log2_expression, na.rm = TRUE),
    sd_log2_expr = sd(log2_expression, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate overall mean and SD across all samples
overall_mean <- mean(sample_mean_expression$mean_log2_expr)
overall_sd <- sd(sample_mean_expression$mean_log2_expr)

# Calculate Z-scores (number of SDs from mean)
sample_mean_expression$z_score <- (sample_mean_expression$mean_log2_expr - overall_mean) / overall_sd

# Identify outliers (|z_score| >= 2)
sample_mean_expression$is_outlier <- abs(sample_mean_expression$z_score) >= 2

# Add metadata information
sample_mean_expression <- sample_mean_expression %>%
  left_join(
    data.frame(sample = rownames(met), met),
    by = "sample"
  )

# Filter to outliers only
outlier_samples <- sample_mean_expression %>%
  filter(is_outlier == TRUE) %>%
  arrange(desc(abs(z_score)))

# Save outlier samples to CSV
if (nrow(outlier_samples) > 0) {
  write.csv(outlier_samples,
            file.path(qc_output_dir, "expression_outlier_samples.csv"),
            row.names = FALSE)
  
  message(sprintf("Found %d outlier sample(s) with |Z-score| >= 2:", nrow(outlier_samples)))
  print(outlier_samples[, c("sample", "mean_log2_expr", "z_score", "batch")], 
        row.names = FALSE)
  
  # Create visualization of outliers
  p_outlier <- ggplot(sample_mean_expression, 
                      aes(x = reorder(sample, mean_log2_expr), 
                          y = mean_log2_expr,
                          fill = is_outlier)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_hline(yintercept = overall_mean, linetype = "dashed", color = "blue", linewidth = 1) +
    geom_hline(yintercept = overall_mean + 2*overall_sd, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_hline(yintercept = overall_mean - 2*overall_sd, linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "red"),
                      labels = c("Normal", "Outlier (≥2 SD)"),
                      name = "") +
    labs(
      title = "Sample Expression Outlier Detection (Timepoint 60, No Outlier)",
      subtitle = sprintf("Mean ± 2SD shown; %d outlier(s) detected", nrow(outlier_samples)),
      x = "Sample",
      y = "Mean log2(expression)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      legend.position = "top"
    )
  
  ggsave(file.path(qc_output_dir, "01d_expression_outlier_detection.png"),
         p_outlier, 
         width = max(12, n_samples * 0.2), 
         height = 7, 
         dpi = 300,
         limitsize = FALSE)
  
  message("Outlier detection plot saved: 01d_expression_outlier_detection.png")
  
} else {
  message("No outlier samples detected (all samples within 2 SD of mean)")
}

# Save all samples with Z-scores for reference
write.csv(sample_mean_expression,
          file.path(qc_output_dir, "all_samples_expression_metrics.csv"),
          row.names = FALSE)

message(sprintf("Overall mean expression: %.3f", overall_mean))
message(sprintf("Overall SD: %.3f", overall_sd))
message(sprintf("2 SD range: [%.3f, %.3f]", 
                overall_mean - 2*overall_sd, 
                overall_mean + 2*overall_sd))

##-----------------------------------------------------------------------
# 3. PCA ANALYSIS
##-----------------------------------------------------------------------

message("\n=== PCA Analysis ===")

keep <- rowSums(counts >= 10) >= 3
counts_filtered <- counts[keep, , drop = FALSE]
message(sprintf("Genes after filtering: %d", nrow(counts_filtered)))

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = met,
  design = ~ 1
)

vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)

ntop <- 23040
rv <- rowVars(vst_mat)
select_genes <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca_result <- prcomp(t(vst_mat[select_genes, , drop = FALSE]))

var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)

pca_data <- as.data.frame(pca_result$x[, 1:5, drop = FALSE])
pca_data$sample <- rownames(pca_data)
pca_data <- merge(
  pca_data,
  data.frame(sample = rownames(met), met),
  by = "sample"
)

cor_matrix <- cor(vst_mat, method = "pearson")

message(sprintf("PCA complete - PC1: %.1f%%, PC2: %.1f%%, PC3: %.1f%%",
                var_explained[1], var_explained[2], var_explained[3]))
##-----------------------------------------------------------------------
# 3.5 COMBINED PCA PLOT (Batch + Treatment + Strain)
##-----------------------------------------------------------------------

message("\n=== Generating Combined PCA Plot (Batch, Treatment, Strain) ===")

# Create separate folder for combined PCA
combined_pca_dir <- file.path(qc_output_dir, "PCA_Batch_Sex_Treatment")
dir.create(combined_pca_dir, showWarnings = FALSE, recursive = TRUE)

# Define deep color palettes (needed here before Section 4)
deep_colors_set1 <- c("#FDE725", "#D55E00", "#6DCD59", "#FF0000", 
                      "#4477AA", "#0000FF", "#5B01A5", "#66CCEE", "#EE3377")

# Extract strain and sex abbreviation from sample names (e.g., B6.F, D2.M)
extract_strain_sex_label <- function(sample_name) {
  # Determine strain
  if (grepl("B6", sample_name, ignore.case = TRUE) || grepl("BL6", sample_name, ignore.case = TRUE)) {
    strain <- "B6"
  } else if (grepl("D2", sample_name, ignore.case = TRUE) || grepl("DBA", sample_name, ignore.case = TRUE)) {
    strain <- "D2"
  } else {
    strain <- "Unknown"
  }
  
  # Determine sex (look for .F. or .M. pattern)
  if (grepl("F", sample_name, ignore.case = TRUE) || grepl("female", sample_name, ignore.case = TRUE)) {
    sex <- ".F"
  } else if (grepl("M", sample_name, ignore.case = TRUE) || grepl("male", sample_name, ignore.case = TRUE)) {
    sex <- ".M"
  } else {
    sex <- ""
  }
  
  return(paste0(strain, sex))
}

pca_data$strain_sex_label <- sapply(pca_data$sample, extract_strain_sex_label)

# Create deep batch colors
n_batches_combined <- length(unique(pca_data$batch))
if (n_batches_combined <= 9) {
  batch_colors_combined <- deep_colors_set1[1:n_batches_combined]
} else {
  batch_colors_combined <- colorRampPalette(deep_colors_set1)(n_batches_combined)
}
names(batch_colors_combined) <- levels(pca_data$batch)

# Define filled shapes for treatment
# 21 = filled circle (vehicle), 24 = filled triangle (THC), 22 = filled square (None)
treatment_shapes <- c("vehicle" = 21, "THC" = 24, "None" = 22)

# Create the combined PCA plot
p_combined_pca <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  # Add ellipses by strain (drawn first, behind points) - both dashed
  stat_ellipse(aes(group = Strain), 
               level = 0.95, 
               linewidth = 1.2, 
               color = "black",
               linetype = "dashed",
               show.legend = FALSE) +
  # Add points with fill by batch, shape by treatment (no outline)
  geom_point(aes(fill = batch, shape = Treatment), 
             size = 3.5, 
             alpha = 0.9) +
  # Add strain.sex labels with repel (plain text, no boxes)
  geom_text_repel(aes(label = strain_sex_label), 
                  size = 3,
                  max.overlaps = 100,
                  segment.size = 0.3,
                  segment.alpha = 0.6,
                  force = 2,
                  fontface = "bold") +
  # Fill scale for batch (inside color)
  scale_fill_manual(values = batch_colors_combined, 
                    name = "RNA Batch (Fill)",
                    guide = guide_legend(order = 1, 
                                         ncol = 1,
                                         override.aes = list(shape = 21, size = 4))) +
  # Shape scale for treatment
  scale_shape_manual(values = treatment_shapes, 
                     name = "Treatment (Shape)",
                     guide = guide_legend(order = 2, 
                                          ncol = 1,
                                          override.aes = list(size = 4))) +
  labs(
    title = "PCA: Batch (Fill) × Treatment (Shape) × Strain (Ellipse)",
    subtitle = sprintf("Timepoint 60, No Outlier | PC1: %.1f%%, PC2: %.1f%% | Dashed ellipses = Strain (95%% CI) | Labels = Strain.Sex", 
                       var_explained[1], var_explained[2]),
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5),
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    legend.spacing.y = unit(0.3, "cm"),
    legend.box = "vertical",
    legend.box.spacing = unit(0.5, "cm"),
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_line(color = "gray92")
  )

# Save the combined PCA plot
ggsave(file.path(combined_pca_dir, "PCA_Batch_Sex_Strain_Treatment_Combined.png"),
       p_combined_pca, width = 13, height = 9, dpi = 300)


##-----------------------------------------------------------------------
# 4. INDIVIDUAL PCA PLOTS (Saved separately in PCA_by_Variable folder)
##-----------------------------------------------------------------------

message("\n=== Generating Individual PCA Plots with Deep Colors ===")

# Define very deep, saturated color palettes
deep_colors_set1 <- c("#8B0000", "#00008B", "#006400", "#8B008B", 
                      "#FF4500", "#000080", "#8B4513", "#2F4F4F", "#8B0000")

deep_colors_set2 <- c("#4B0082", "#8B0000", "#006400", "#00008B", 
                      "#8B4500", "#4B0082", "#8B008B", "#2F4F4F", "#191970")

deep_colors_set3 <- c("#800000", "#000080", "#2F4F4F", "#8B4513", 
                      "#556B2F", "#8B0000", "#4B0082", "#00008B", "#8B008B")

# Create very deep batch colors
n_batches_pca <- length(unique(pca_data$batch))
if (n_batches_pca <= 9) {
  deep_batch_colors <- deep_colors_set1[1:n_batches_pca]
} else {
  deep_batch_colors <- colorRampPalette(deep_colors_set1)(n_batches_pca)
}
names(deep_batch_colors) <- levels(pca_data$batch)

create_deep_pca_plot <- function(data, color_var, var_exp, title, color_palette = NULL) {
  p <- ggplot(data, aes_string(x = "PC1", y = "PC2", color = color_var)) +
    geom_point(size = 3.5, alpha = 0.9) +  # Larger points, higher alpha for intensity
    geom_text_repel(aes(label = sample), 
                    size = 1.8,
                    max.overlaps = 100,
                    segment.size = 0.3,
                    segment.alpha = 0.6,
                    alpha = 0.9,
                    force = 2) +
    stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.8, alpha = 0.6,
                 show.legend = FALSE) +
    labs(
      title = title,
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      legend.text = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      panel.grid.major = element_line(color = "gray85"),
      panel.grid.minor = element_line(color = "gray92")
    )
  
  # Apply custom color palette if provided, otherwise use defaults
  if (!is.null(color_palette)) {
    p <- p + scale_color_manual(values = color_palette)
  } else {
    n_levels <- length(unique(data[[color_var]]))
    if (n_levels <= 9) {
      p <- p + scale_color_manual(values = deep_colors_set1[1:n_levels])
    } else {
      p <- p + scale_color_manual(values = colorRampPalette(deep_colors_set1)(n_levels))
    }
  }
  p
}

# 4.1 PCA by RNA.Batch - very deep colors
message("  Creating PCA plot for RNA.Batch")
p_pca_batch <- create_deep_pca_plot(pca_data, "batch", var_explained, 
                                    "PCA: Colored by RNA Batch (TP60, No Outlier)",
                                    color_palette = deep_batch_colors)
ggsave(file.path(pca_by_variable_dir, "PCA_RNA_Batch.png"),
       p_pca_batch, width = 10, height = 8, dpi = 300)

# 4.2 PCA by Treatment - very deep colors
if ("Treatment" %in% colnames(pca_data)) {
  message("  Creating PCA plot for Treatment")
  n_treatment <- length(unique(pca_data$Treatment))
  deep_treatment_colors <- c("vehicle" = "#00008B", "THC" = "#8B0000", "None" = "#006400")
  # Filter to only include levels that exist
  deep_treatment_colors <- deep_treatment_colors[names(deep_treatment_colors) %in% levels(pca_data$Treatment)]
  
  p_pca_treatment <- create_deep_pca_plot(pca_data, "Treatment", var_explained, 
                                          "PCA: Colored by Treatment (TP60, No Outlier)",
                                          color_palette = deep_treatment_colors)
  ggsave(file.path(pca_by_variable_dir, "PCA_Treatment.png"),
         p_pca_treatment, width = 10, height = 8, dpi = 300)
}

# 4.3 PCA by Sex - very deep colors
if ("Sex" %in% colnames(pca_data)) {
  message("  Creating PCA plot for Sex")
  deep_sex_colors <- c("female" = "#8B008B", "male" = "#00008B")
  
  p_pca_sex <- create_deep_pca_plot(pca_data, "Sex", var_explained, 
                                    "PCA: Colored by Sex (TP60, No Outlier)",
                                    color_palette = deep_sex_colors)
  ggsave(file.path(pca_by_variable_dir, "PCA_Sex.png"),
         p_pca_sex, width = 10, height = 8, dpi = 300)
}

# 4.4 PCA by Strain - very deep colors
if ("Strain" %in% colnames(pca_data)) {
  message("  Creating PCA plot for Strain")
  n_strain <- length(unique(pca_data$Strain))
  deep_strain_colors <- colorRampPalette(c("#8B0000", "#00008B", "#006400", 
                                           "#8B008B", "#8B4513", "#2F4F4F"))(n_strain)
  names(deep_strain_colors) <- levels(pca_data$Strain)
  
  p_pca_strain <- create_deep_pca_plot(pca_data, "Strain", var_explained, 
                                       "PCA: Colored by Strain (TP60, No Outlier)",
                                       color_palette = deep_strain_colors)
  ggsave(file.path(pca_by_variable_dir, "PCA_Strain.png"),
         p_pca_strain, width = 10, height = 8, dpi = 300)
}

# 4.5 PCA by Library.Batch - very deep colors
if ("Library.Batch" %in% colnames(pca_data)) {
  message("  Creating PCA plot for Library.Batch")
  n_lib_batches <- length(unique(pca_data$Library.Batch))
  deep_lib_colors <- colorRampPalette(c("#4B0082", "#8B0000", "#0000FF", 
                                        "#FF0000", "#8B4513", "#2F4F4F"))(n_lib_batches)
  names(deep_lib_colors) <- levels(pca_data$Library.Batch)
  
  p_pca_libbatch <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Library.Batch)) +
    geom_point(size = 3.5, alpha = 0.9) +
    geom_text_repel(aes(label = sample), size = 1.8, max.overlaps = 100,
                    segment.size = 0.3, segment.alpha = 0.6, alpha = 0.9, force = 2) +
    stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.8, alpha = 0.6, 
                 show.legend = FALSE) +
    scale_color_manual(values = deep_lib_colors) +
    labs(title = "PCA: Colored by Library Batch (TP60, No Outlier)",
         x = paste0("PC1 (", var_explained[1], "%)"),
         y = paste0("PC2 (", var_explained[2], "%)")) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      legend.text = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      panel.grid.major = element_line(color = "gray85"),
      panel.grid.minor = element_line(color = "gray92")
    )
  
  ggsave(file.path(pca_by_variable_dir, "PCA_Library_Batch.png"),
         p_pca_libbatch, width = 10, height = 8, dpi = 300)
}

message(sprintf("Individual PCA plots with deep colors saved in: %s", pca_by_variable_dir))

##-----------------------------------------------------------------------
# 4.6 INDIVIDUAL GROUP PCA PLOTS (group1-group7)
##-----------------------------------------------------------------------

message("\n=== Generating Individual Group PCA Plots ===")

for (i in 1:7) {
  group_col <- paste0("Group", i)
  
  if (group_col %in% colnames(pca_data)) {
    message(sprintf("Creating PCA plot for %s", group_col))
    
    n_group_levels <- length(unique(pca_data[[group_col]]))
    
    # Choose color palette based on number of levels
    if (n_group_levels <= 9) {
      group_colors <- brewer.pal(min(9, n_group_levels), "Set1")[1:n_group_levels]
    } else {
      group_colors <- colorRampPalette(brewer.pal(9, "Set1"))(n_group_levels)
    }
    names(group_colors) <- levels(pca_data[[group_col]])
    
    p_group <- ggplot(pca_data, aes_string(x = "PC1", y = "PC2", color = group_col)) +
      geom_point(size = 3, alpha = 0.8) +
      geom_text_repel(aes(label = sample), 
                      size = 2,
                      max.overlaps = 100,
                      segment.size = 0.2,
                      segment.alpha = 0.5,
                      alpha = 0.8,
                      force = 2) +
      stat_ellipse(level = 0.95, linetype = "dashed", alpha = 0.5,
                   show.legend = FALSE) +
      scale_color_manual(values = group_colors) +
      labs(
        title = sprintf("PCA: Colored by %s (TP60, No Outlier)", group_col),
        subtitle = sprintf("%d levels: %s", 
                           n_group_levels, 
                           paste(levels(pca_data[[group_col]]), collapse = ", ")),
        x = paste0("PC1 (", var_explained[1], "%)"),
        y = paste0("PC2 (", var_explained[2], "%)"),
        color = group_col
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10),
        legend.position = "bottom",
        legend.text = element_text(size = 9)
      ) +
      guides(color = guide_legend(nrow = ceiling(n_group_levels / 5)))
    
    # Save individual group PCA plot
    ggsave(file.path(group_pca_dir, sprintf("PCA_%s.png", group_col)),
           p_group, width = 10, height = 8, dpi = 300)
    
    message(sprintf("  %s PCA plot saved", group_col))
  } else {
    message(sprintf("  %s not found in data, skipping", group_col))
  }
}

message(sprintf("Individual group PCA plots saved in: %s", group_pca_dir))

##-----------------------------------------------------------------------
# 5. CORRELATION HEATMAP
##-----------------------------------------------------------------------

message("\n=== Generating Correlation Heatmap ===")

annotation_vars <- intersect(c("batch", "RNA.Batch", "Library.Batch", 
                               "Treatment", "Sex", "Strain",
                               "Group1", "Group2", "Group3", "Group4", 
                               "Group5", "Group6", "Group7"),
                             colnames(met))
annotation_df <- met[, annotation_vars, drop = FALSE]

ann_colors <- list()
if ("batch" %in% colnames(annotation_df)) {
  ann_colors$batch <- batch_colors
}
if ("RNA.Batch" %in% colnames(annotation_df)) {
  ann_colors$RNA.Batch <- batch_colors
}
if ("Library.Batch" %in% colnames(annotation_df)) {
  n_lib <- length(unique(annotation_df$Library.Batch))
  ann_colors$Library.Batch <- setNames(colorRampPalette(brewer.pal(9, "Set2"))(n_lib),
                                       levels(annotation_df$Library.Batch))
}
if ("Treatment" %in% colnames(annotation_df)) {
  ann_colors$Treatment <- c("vehicle" = "#3182bd", "THC" = "#e6550d", "None" = "green")
}
if ("Sex" %in% colnames(annotation_df)) {
  ann_colors$Sex <- c("female" = "#ff69b4", "male" = "#4169e1")
}
if ("Strain" %in% colnames(annotation_df)) {
  n_strain <- length(levels(annotation_df$Strain))
  ann_colors$Strain <- setNames(brewer.pal(min(8, n_strain), "Dark2")[1:n_strain],
                                levels(annotation_df$Strain))
}

# Add colors for group1-group7
color_palettes <- list("Set1", "Set2", "Set3", "Paired", "Dark2", "Accent", "Pastel1")
for (i in 1:7) {
  group_col <- paste0("Group", i)
  if (group_col %in% colnames(annotation_df)) {
    n_levels <- length(levels(annotation_df[[group_col]]))
    palette <- color_palettes[[((i-1) %% length(color_palettes)) + 1]]
    if (n_levels <= 9) {
      ann_colors[[group_col]] <- setNames(
        brewer.pal(min(n_levels, 9), palette)[1:n_levels],
        levels(annotation_df[[group_col]])
      )
    } else {
      ann_colors[[group_col]] <- setNames(
        colorRampPalette(brewer.pal(9, palette))(n_levels),
        levels(annotation_df[[group_col]])
      )
    }
  }
}

fontsize_samples <- max(3, min(8, 200 / n_samples))

png(file.path(qc_output_dir, "03_correlation_heatmap.png"), 
    width = 18, height = 16, units = "in", res = 300)
pheatmap(
  cor_matrix,
  main = "Sample Correlation Heatmap (Timepoint 60, No Outlier)",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_df,
  annotation_row = annotation_df,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 6,
  fontsize_row = fontsize_samples,
  fontsize_col = fontsize_samples,
  color = colorRampPalette(c("#0000FF", "#f7f7f7", "#FF0000"))(100),
  breaks = seq(min(cor_matrix), 1, length.out = 101),
  border_color = NA
)
dev.off()

message("Correlation heatmap saved")

##-----------------------------------------------------------------------
# 6. PC-VARIABLE ASSOCIATION (ANOVA)
##-----------------------------------------------------------------------

message("\n=== PC-Variable Association Analysis (ANOVA) ===")

calc_pc_associations <- function(pca_data) {
  test_vars <- intersect(c("RNA.Batch", "Library.Batch", 
                           "Treatment", "Sex", "Strain",
                           "Group1", "Group2", "Group3", "Group4",
                           "Group5", "Group6", "Group7"),
                         colnames(pca_data))
  
  results <- data.frame(
    Variable = character(),
    PC1_pvalue = numeric(),
    PC2_pvalue = numeric(),
    PC3_pvalue = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (var in test_vars) {
    if (length(unique(pca_data[[var]])) > 1) {
      p1 <- summary(aov(pca_data$PC1 ~ pca_data[[var]]))[[1]][["Pr(>F)"]][1]
      p2 <- summary(aov(pca_data$PC2 ~ pca_data[[var]]))[[1]][["Pr(>F)"]][1]
      p3 <- summary(aov(pca_data$PC3 ~ pca_data[[var]]))[[1]][["Pr(>F)"]][1]
      
      results <- rbind(results, data.frame(
        Variable = var,
        PC1_pvalue = p1,
        PC2_pvalue = p2,
        PC3_pvalue = p3
      ))
    }
  }
  results
}

pc_associations <- calc_pc_associations(pca_data)

assoc_long <- pc_associations %>%
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "pvalue") %>%
  mutate(
    neg_log_p = -log10(pvalue + 1e-300),
    PC = gsub("_pvalue", "", PC)
  )

p_associations <- ggplot(assoc_long, aes(x = PC, y = Variable, fill = neg_log_p)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1e", pvalue)), size = 3) +
  scale_fill_gradient2(low = "white", mid = "yellow", high = "red",
                       midpoint = 2, name = "-log10(p)") +
  labs(
    title = "PC-Variable Association (ANOVA) - Timepoint 60, No Outlier",
    subtitle = "Lower p-values indicate stronger association with PCs",
    x = "Principal Component",
    y = "Variable"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 10)
  )

ggsave(file.path(qc_output_dir, "04_PC_associations_anova.png"),
       p_associations, width = 10, height = max(6, nrow(pc_associations) * 0.4), dpi = 300)

write.csv(pc_associations,
          file.path(qc_output_dir, "PC_variable_associations_anova.csv"),
          row.names = FALSE)

message("PC-variable associations (ANOVA) saved")

##-----------------------------------------------------------------------
# 7. PC2 CORRELATION ANALYSIS (Spearman/Pearson)
##-----------------------------------------------------------------------

message("\n=== PC2 Correlation Analysis ===")

# Define variable types
categorical_vars <- c("Treatment", "Sex", "Strain", "RNA.Batch", "Library.Batch",
                      "Group1", "Group2", "Group3", "Group4", "Group5", "Group6", "Group7")
continuous_vars <- c("RIN", "Conc")

# Get variables that actually exist in metadata
categorical_vars <- intersect(categorical_vars, colnames(met))
continuous_vars <- intersect(continuous_vars, colnames(met))

message(sprintf("Categorical variables found: %s", paste(categorical_vars, collapse = ", ")))
message(sprintf("Continuous variables found: %s", paste(continuous_vars, collapse = ", ")))

# Initialize results dataframe
pc2_results <- data.frame(
  Variable = character(),
  Type = character(),
  Method = character(),
  Correlation = numeric(),
  R_squared = numeric(),
  P_value = numeric(),
  N_levels = integer(),
  Encoding = character(),
  stringsAsFactors = FALSE
)

pc2_values <- pca_data$PC2
names(pc2_values) <- pca_data$sample

# --- Categorical variables (Spearman with encoding) ---
for (var in categorical_vars) {
  if (var %in% colnames(pca_data)) {
    var_values <- pca_data[[var]]
  } else {
    next
  }
  
  # Convert to factor if not already
  if (!is.factor(var_values)) {
    var_values <- factor(var_values)
  }
  
  n_levels <- length(levels(var_values))
  
  # Encode as numeric (1, 2, 3, ...)
  var_encoded <- as.numeric(var_values)
  
  # Create encoding key for reference
  encoding_key <- paste(levels(var_values), "=", seq_along(levels(var_values)), collapse = "; ")
  
  # Spearman correlation
  cor_test <- cor.test(pc2_values, var_encoded, method = "spearman", exact = FALSE)
  
  pc2_results <- rbind(pc2_results, data.frame(
    Variable = var,
    Type = "Categorical",
    Method = "Spearman",
    Correlation = cor_test$estimate,
    R_squared = cor_test$estimate^2,
    P_value = cor_test$p.value,
    N_levels = n_levels,
    Encoding = encoding_key,
    stringsAsFactors = FALSE
  ))
}

# --- Continuous variables (Pearson) ---
for (var in continuous_vars) {
  if (var %in% colnames(pca_data)) {
    var_values <- pca_data[[var]]
  } else {
    next
  }
  
  # Check if truly numeric
  if (!is.numeric(var_values)) {
    var_values <- as.numeric(as.character(var_values))
  }
  
  # Remove NAs for correlation
  valid_idx <- !is.na(var_values) & !is.na(pc2_values)
  
  if (sum(valid_idx) < 3) {
    message(sprintf("  Warning: %s has too few valid observations, skipping", var))
    next
  }
  
  # Pearson correlation
  cor_test <- cor.test(pc2_values[valid_idx], var_values[valid_idx], method = "pearson")
  
  pc2_results <- rbind(pc2_results, data.frame(
    Variable = var,
    Type = "Continuous",
    Method = "Pearson",
    Correlation = cor_test$estimate,
    R_squared = cor_test$estimate^2,
    P_value = cor_test$p.value,
    N_levels = NA,
    Encoding = "N/A (continuous)",
    stringsAsFactors = FALSE
  ))
}

# Add significance indicators
pc2_results$Significance <- ifelse(pc2_results$P_value < 0.001, "***",
                                   ifelse(pc2_results$P_value < 0.01, "**",
                                          ifelse(pc2_results$P_value < 0.05, "*", "ns")))

# Round for readability
pc2_results$Correlation <- round(pc2_results$Correlation, 4)
pc2_results$R_squared <- round(pc2_results$R_squared, 4)

# Reorder columns
pc2_results <- pc2_results[, c("Variable", "Type", "Method", "Correlation", 
                               "R_squared", "P_value", "Significance", "N_levels", "Encoding")]

# Save to CSV
write.csv(pc2_results, 
          file.path(qc_output_dir, "PC2_correlations_with_metadata.csv"),
          row.names = FALSE)

message("\nPC2 correlation results:")
print(pc2_results[, c("Variable", "Type", "Correlation", "R_squared", "P_value", "Significance")], 
      row.names = FALSE)

##-----------------------------------------------------------------------
# 8. PC2 CORRELATION VISUALIZATION
##-----------------------------------------------------------------------

message("\n=== Generating PC2 Correlation Visualization ===")

# Bar plot of R² values
p_r2_bars <- ggplot(pc2_results, aes(x = reorder(Variable, R_squared), y = R_squared, 
                                     fill = Type)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.85) +
  geom_text(aes(label = Significance, y = R_squared + 0.01), 
            size = 4, vjust = 0) +
  geom_text(aes(label = sprintf("r=%.2f", Correlation), y = R_squared / 2),
            size = 3, color = "white", fontface = "bold") +
  scale_fill_manual(values = c("Categorical" = "#FF0000", "Continuous" = "#0000FF")) +
  labs(
    title = "PC2 Correlation with Metadata Variables (TP60, No Outlier)",
    subtitle = "R² values shown; r values inside bars; * p<0.05, ** p<0.01, *** p<0.001",
    x = "Variable",
    y = expression(R^2),
    fill = "Variable Type"
  ) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    axis.text = element_text(size = 10),
    legend.position = "top"
  ) +
  ylim(0, max(pc2_results$R_squared, na.rm = TRUE) * 1.15)

ggsave(file.path(qc_output_dir, "05a_PC2_correlations_barplot.png"),
       p_r2_bars, width = 10, height = max(7, nrow(pc2_results) * 0.35), dpi = 300)

# Heatmap version
p_heatmap <- ggplot(pc2_results, aes(x = "PC2", y = reorder(Variable, abs(Correlation)))) +
  geom_tile(aes(fill = Correlation), color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f%s", Correlation, Significance)), 
            size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#0000FF", mid = "white", high = "#FF0000",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Correlation (r)") +
  labs(
    title = "PC2 Correlation Heatmap (TP60, No Outlier)",
    subtitle = "Values show correlation coefficient; * p<0.05, ** p<0.01, *** p<0.001",
    x = "",
    y = "Metadata Variable"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(face = "bold", size = 12),
    panel.grid = element_blank()
  )

ggsave(file.path(qc_output_dir, "05b_PC2_correlations_heatmap.png"),
       p_heatmap, width = 6, height = max(7, nrow(pc2_results) * 0.35), dpi = 300)

# Combined figure
combined_pc2 <- (p_r2_bars | p_heatmap) +
  plot_annotation(
    title = "PC2 Association with Metadata Variables (TP60, No Outlier)",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

ggsave(file.path(qc_output_dir, "05_PC2_correlations_combined.png"),
       combined_pc2, width = 16, height = max(8, nrow(pc2_results) * 0.35), dpi = 300)

message("PC2 correlation visualizations saved")

# Save encoding reference
encoding_ref <- pc2_results %>%
  filter(Type == "Categorical") %>%
  dplyr::select(Variable, N_levels, Encoding)

write.csv(encoding_ref,
          file.path(qc_output_dir, "PC2_categorical_encoding_reference.csv"),
          row.names = FALSE)

##-----------------------------------------------------------------------
# 9. LIBRARY SIZE QC
##-----------------------------------------------------------------------

message("\n=== Library Size QC ===")

lib_sizes <- colSums(counts)
lib_size_df <- data.frame(
  sample = names(lib_sizes),
  total_counts = lib_sizes,
  millions = lib_sizes / 1e6
)

lib_size_df <- merge(lib_size_df,
                     data.frame(sample = rownames(met), met),
                     by = "sample")

p_libsize_batch <- ggplot(lib_size_df, aes(x = batch, y = millions, fill = batch)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = batch_colors) +
  labs(title = "Library Size by RNA Batch (TP60, No Outlier)", 
       x = "RNA Batch", 
       y = "Total Counts (Millions)") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

if ("Treatment" %in% colnames(lib_size_df)) {
  p_libsize_treatment <- ggplot(lib_size_df, aes(x = Treatment, y = millions, fill = Treatment)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "Library Size by Treatment (TP60, No Outlier)", 
         x = "Treatment", 
         y = "Total Counts (Millions)") +
    theme_bw() +
    theme(legend.position = "none")
} else {
  p_libsize_treatment <- ggplot() + theme_void() + labs(title = "Treatment not found")
}

combined_qc <- (p_libsize_batch | p_libsize_treatment) +
  plot_annotation(
    title = "Technical QC Metrics (Timepoint 60, No Outlier)",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

ggsave(file.path(qc_output_dir, "06_technical_qc_metrics.png"),
       combined_qc, width = 12, height = 6, dpi = 300)

##-----------------------------------------------------------------------
# 10. SAVE FINAL DATASETS WITH PC2
##-----------------------------------------------------------------------

message("\n=== Saving Final Datasets with PC2 ===")

# Add PC2 to metadata
met$PC2 <- pca_data$PC2[match(rownames(met), pca_data$sample)]

# Save metadata with PC2
write.csv(met, 
          file.path(qc_output_dir, "metv2_timepoint60_with_PC2_no_outlier.csv"),
          row.names = TRUE)

# Save count data
write.csv(counts,
          file.path(qc_output_dir, "counts_timepoint60_no_outlier.csv"),
          row.names = TRUE)

message("Final datasets saved:")
message("  - metv2_timepoint60_with_PC2_no_outlier.csv")
message("  - counts_timepoint60_no_outlier.csv")

##-----------------------------------------------------------------------
# 11. SUMMARY REPORT
##-----------------------------------------------------------------------

message("\n=== Generating Summary Report ===")

upper_tri <- cor_matrix[upper.tri(cor_matrix)]

summary_df <- data.frame(
  Metric = c(
    "Total Samples (Timepoint 60, No Outlier)",
    "Total Genes (raw)",
    "Genes after filtering",
    "Number of RNA Batches",
    "Number of Library Batches",
    "PC1 variance (%)",
    "PC2 variance (%)",
    "PC3 variance (%)",
    "Mean sample correlation",
    "Min sample correlation"
  ),
  Value = c(
    ncol(counts),
    nrow(counts),
    nrow(counts_filtered),
    length(unique(met$RNA.Batch)),
    ifelse("Library.Batch" %in% colnames(met), length(unique(met$Library.Batch)), "N/A"),
    var_explained[1],
    var_explained[2],
    var_explained[3],
    round(mean(upper_tri), 4),
    round(min(upper_tri), 4)
  )
)

write.csv(summary_df, file.path(qc_output_dir, "QC_summary.csv"), row.names = FALSE)
write.csv(lib_size_df, file.path(qc_output_dir, "library_sizes.csv"), row.names = FALSE)
write.csv(pca_data, file.path(qc_output_dir, "PCA_coordinates.csv"), row.names = FALSE)

