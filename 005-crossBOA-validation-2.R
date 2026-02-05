library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggh4x)  # For label alignment in circular plots
library(RColorBrewer)  # For color scheme optimization

# Create output directory
output_dir <- 'results/001-cross-sectional-top100-gene-check'
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# =============================================================================
# Step 1: Read both_significant_consistent_direction table, get top 100 seqid list
# =============================================================================
cat("Step 1: Reading both_significant_consistent_direction table...\n")
consistent_file <- 'results/1227-result/usa-iceland/both_significant_consistent_direction.csv'
consistent_df <- read.csv(consistent_file, stringsAsFactors = FALSE)

cat("Consistent direction table dimensions:", dim(consistent_df), "\n")
cat("Consistent direction table column names:", colnames(consistent_df), "\n")

# Create label_symbol to identify duplicate genes
consistent_df <- consistent_df %>%
  mutate(
    label_symbol = ifelse(!is.na(EntrezGeneSymbol) & EntrezGeneSymbol != "", 
                         EntrezGeneSymbol, 
                         ifelse(!is.na(gene_symbol) & gene_symbol != "", 
                                gene_symbol, 
                                seq_id))
  )

# For duplicate gene names, keep only the row with smallest rank_sum (highest rank)
cat("\nProcessing duplicate gene names...\n")
cat("Count before processing:", nrow(consistent_df), "\n")
consistent_df_unique <- consistent_df %>%
  group_by(label_symbol) %>%
  slice_min(rank_sum, n = 1) %>%
  ungroup()

cat("Count after processing:", nrow(consistent_df_unique), "\n")

# Get top 100 genes (based on rank_sum)
top100_genes <- consistent_df_unique %>%
  arrange(rank_sum) %>%
  head(100)

top100_seqids <- top100_genes$seq_id
top100_symbols <- top100_genes$label_symbol

cat("\nTop 100 genes rank range:", range(top100_genes$rank_sum, na.rm = TRUE), "\n")
cat("Top 100 seqid count:", length(top100_seqids), "\n")

# =============================================================================
# Step 2: Read mixed data (age, sex, individual, residual)
# =============================================================================
cat("\nStep 2: Reading mixed data...\n")
mixed_file <- 'results/1227-result/variance-decomposition-survival-analysis-merged.csv'
mixed_df <- read.csv(mixed_file, stringsAsFactors = FALSE)

cat("Mixed data table dimensions:", dim(mixed_df), "\n")
cat("Mixed data table column names:", colnames(mixed_df), "\n")

# Create gene_symbol column (for deduplication)
mixed_df <- mixed_df %>%
  mutate(seq_id = seqid)

# Check and create gene_symbol column
if ("EntrezGeneSymbol" %in% colnames(mixed_df)) {
  mixed_df <- mixed_df %>%
    mutate(
      gene_symbol = ifelse(!is.na(gene_symbol) & gene_symbol != "", 
                          gene_symbol, 
                          ifelse(!is.na(EntrezGeneSymbol) & EntrezGeneSymbol != "", 
                                 EntrezGeneSymbol, 
                                 seq_id))
    )
} else {
  mixed_df <- mixed_df %>%
    mutate(
      gene_symbol = ifelse(!is.na(gene_symbol) & gene_symbol != "", 
                          gene_symbol, 
                          seq_id)
    )
}

# Deduplicate by gene_symbol, keep row with largest age contribution
cat("\nDeduplicating by gene_symbol (keep row with largest age contribution)...\n")
cat("Count before deduplication:", nrow(mixed_df), "\n")
mixed_df <- mixed_df %>%
  filter(!is.na(prop_age)) %>%
  group_by(gene_symbol) %>%
  slice_max(prop_age, n = 1, with_ties = FALSE) %>%
  ungroup()
cat("Count after deduplication:", nrow(mixed_df), "\n")

# Use top100_symbols to filter top 100 genes
mixed_top100 <- mixed_df %>%
  filter(gene_symbol %in% top100_symbols) %>%
  mutate(
    label_symbol = gene_symbol  # Since filtered by top100_symbols, gene_symbol is label_symbol
  ) %>%
  arrange(match(gene_symbol, top100_symbols))

cat("Number of top100 matches in mixed data:", nrow(mixed_top100), "\n")

# =============================================================================
# Step 3: Read male data (age, individual, residual)
# =============================================================================
cat("\nStep 3: Reading male data...\n")
male_file <- 'results/1227-result-sex-separated/combined_protein_analysis_results_male.csv'
male_df <- read.csv(male_file, stringsAsFactors = FALSE)

cat("Male data table dimensions:", dim(male_df), "\n")
cat("Male data table column names:", colnames(male_df), "\n")

# Filter top100 seqids (need to confirm seqid column name in male data)
# Assume column name may be seqid or seq_id
seqid_col_male <- ifelse("seqid" %in% colnames(male_df), "seqid", 
                        ifelse("seq_id" %in% colnames(male_df), "seq_id", NA))

if (is.na(seqid_col_male)) {
  stop("Cannot find seqid or seq_id column in male data")
}

# Use dynamic column name
male_df$seq_id <- male_df[[seqid_col_male]]

# Create gene_symbol column (for deduplication)
if ("EntrezGeneSymbol" %in% colnames(male_df)) {
  male_df <- male_df %>%
    mutate(
      gene_symbol = ifelse(!is.na(gene_symbol) & gene_symbol != "", 
                          gene_symbol, 
                          ifelse(!is.na(EntrezGeneSymbol) & EntrezGeneSymbol != "", 
                                 EntrezGeneSymbol, 
                                 seq_id))
    )
} else {
  male_df <- male_df %>%
    mutate(
      gene_symbol = ifelse(!is.na(gene_symbol) & gene_symbol != "", 
                          gene_symbol, 
                          seq_id)
    )
}

# Deduplicate by gene_symbol, keep row with largest age contribution
cat("\nDeduplicating by gene_symbol (keep row with largest age contribution)...\n")
cat("Count before deduplication:", nrow(male_df), "\n")
male_df <- male_df %>%
  filter(!is.na(prop_age)) %>%
  group_by(gene_symbol) %>%
  slice_max(prop_age, n = 1, with_ties = FALSE) %>%
  ungroup()
cat("Count after deduplication:", nrow(male_df), "\n")

# Use top100_symbols to filter top 100 genes
male_top100 <- male_df %>%
  filter(gene_symbol %in% top100_symbols) %>%
  mutate(
    label_symbol = gene_symbol  # Since filtered by top100_symbols, gene_symbol is label_symbol
  ) %>%
  arrange(match(gene_symbol, top100_symbols))

cat("Number of top100 matches in male data:", nrow(male_top100), "\n")

# =============================================================================
# Step 4: Read female data (age, individual, residual)
# =============================================================================
cat("\nStep 4: Reading female data...\n")
female_file <- 'results/1227-result-sex-separated/combined_protein_analysis_results_female.csv'
female_df <- read.csv(female_file, stringsAsFactors = FALSE)

cat("Female data table dimensions:", dim(female_df), "\n")
cat("Female data table column names:", colnames(female_df), "\n")

# Filter top100 seqids
seqid_col_female <- ifelse("seqid" %in% colnames(female_df), "seqid", 
                          ifelse("seq_id" %in% colnames(female_df), "seq_id", NA))

if (is.na(seqid_col_female)) {
  stop("Cannot find seqid or seq_id column in female data")
}

# Use dynamic column name
female_df$seq_id <- female_df[[seqid_col_female]]

# Create gene_symbol column (for deduplication)
if ("EntrezGeneSymbol" %in% colnames(female_df)) {
  female_df <- female_df %>%
    mutate(
      gene_symbol = ifelse(!is.na(gene_symbol) & gene_symbol != "", 
                          gene_symbol, 
                          ifelse(!is.na(EntrezGeneSymbol) & EntrezGeneSymbol != "", 
                                 EntrezGeneSymbol, 
                                 seq_id))
    )
} else {
  female_df <- female_df %>%
    mutate(
      gene_symbol = ifelse(!is.na(gene_symbol) & gene_symbol != "", 
                          gene_symbol, 
                          seq_id)
    )
}

# Deduplicate by gene_symbol, keep row with largest age contribution
cat("\nDeduplicating by gene_symbol (keep row with largest age contribution)...\n")
cat("Count before deduplication:", nrow(female_df), "\n")
female_df <- female_df %>%
  filter(!is.na(prop_age)) %>%
  group_by(gene_symbol) %>%
  slice_max(prop_age, n = 1, with_ties = FALSE) %>%
  ungroup()
cat("Count after deduplication:", nrow(female_df), "\n")

# Use top100_symbols to filter top 100 genes
female_top100 <- female_df %>%
  filter(gene_symbol %in% top100_symbols) %>%
  mutate(
    label_symbol = gene_symbol  # Since filtered by top100_symbols, gene_symbol is label_symbol
  ) %>%
  arrange(match(gene_symbol, top100_symbols))

cat("Number of top100 matches in female data:", nrow(female_top100), "\n")

# =============================================================================
# Step 5: Prepare plotting data and create three stacked barplots
# =============================================================================

# Function: Prepare stacked barplot data
prepare_stacked_data <- function(df, variance_cols, col_names) {
  df %>%
    select(label_symbol, all_of(variance_cols)) %>%
    filter(!is.na(label_symbol)) %>%
    tidyr::pivot_longer(
      cols = all_of(variance_cols),
      names_to = "variance_component",
      values_to = "proportion"
    ) %>%
    mutate(
      variance_component = case_when(
        variance_component == variance_cols[1] ~ col_names[1],
        variance_component == variance_cols[2] ~ col_names[2],
        variance_component == variance_cols[3] ~ col_names[3],
        variance_component == variance_cols[4] ~ col_names[4],
        TRUE ~ variance_component
      )
    )
}

# Plot 1: Mixed data (age, sex, individual, residual)
cat("\nStep 5: Preparing mixed data for plotting...\n")
mixed_cols <- c("prop_age", "prop_sex", "prop_idno", "prop_residual")
mixed_col_names <- c("Age", "Sex", "Individual", "Residual")

# Check if columns exist
missing_mixed_cols <- setdiff(mixed_cols, colnames(mixed_top100))
if (length(missing_mixed_cols) > 0) {
  cat("Warning: Missing columns in mixed data:", paste(missing_mixed_cols, collapse = ", "), "\n")
  cat("Available columns:", paste(colnames(mixed_top100), collapse = ", "), "\n")
  # Try to find similar column names
  for (col in missing_mixed_cols) {
    similar_cols <- grep(col, colnames(mixed_top100), ignore.case = TRUE, value = TRUE)
    if (length(similar_cols) > 0) {
      cat("  Possible matches:", paste(similar_cols, collapse = ", "), "\n")
    }
  }
}

mixed_long <- mixed_top100 %>%
  select(label_symbol, all_of(intersect(mixed_cols, colnames(mixed_top100)))) %>%
  filter(!is.na(label_symbol)) %>%
  mutate(label_symbol = factor(label_symbol, levels = top100_symbols)) %>%
  tidyr::pivot_longer(
    cols = all_of(intersect(mixed_cols, colnames(mixed_top100))),
    names_to = "variance_component",
    values_to = "proportion"
  ) %>%
  mutate(
    variance_component = case_when(
      variance_component == "prop_age" ~ "Age",
      variance_component == "prop_sex" ~ "Sex",
      variance_component == "prop_idno" ~ "Individual",
      variance_component == "prop_residual" ~ "Residual",
      TRUE ~ variance_component
    ),
    variance_component = factor(
      variance_component,
      levels = c("Residual", "Individual", "Sex", "Age")
    )
  )

# Plot 2: Male data (age, individual, residual)
cat("\nPreparing male data for plotting...\n")
male_cols <- c("prop_age", "prop_idno", "prop_residual")
male_col_names <- c("Age", "Individual", "Residual")

missing_male_cols <- setdiff(male_cols, colnames(male_top100))
if (length(missing_male_cols) > 0) {
  cat("Warning: Missing columns in male data:", paste(missing_male_cols, collapse = ", "), "\n")
  cat("Available columns:", paste(colnames(male_top100), collapse = ", "), "\n")
}

male_long <- male_top100 %>%
  select(label_symbol, all_of(intersect(male_cols, colnames(male_top100)))) %>%
  filter(!is.na(label_symbol)) %>%
  mutate(label_symbol = factor(label_symbol, levels = top100_symbols)) %>%
  tidyr::pivot_longer(
    cols = all_of(intersect(male_cols, colnames(male_top100))),
    names_to = "variance_component",
    values_to = "proportion"
  ) %>%
  mutate(
    variance_component = case_when(
      variance_component == "prop_age" ~ "Age",
      variance_component == "prop_idno" ~ "Individual",
      variance_component == "prop_residual" ~ "Residual",
      TRUE ~ variance_component
    ),
    variance_component = factor(
      variance_component,
      levels = c("Residual", "Individual", "Age")
    )
  )

# Plot 3: Female data (age, individual, residual)
cat("\nPreparing female data for plotting...\n")
female_cols <- c("prop_age", "prop_idno", "prop_residual")
female_col_names <- c("Age", "Individual", "Residual")

missing_female_cols <- setdiff(female_cols, colnames(female_top100))
if (length(missing_female_cols) > 0) {
  cat("Warning: Missing columns in female data:", paste(missing_female_cols, collapse = ", "), "\n")
  cat("Available columns:", paste(colnames(female_top100), collapse = ", "), "\n")
}

female_long <- female_top100 %>%
  select(label_symbol, all_of(intersect(female_cols, colnames(female_top100)))) %>%
  filter(!is.na(label_symbol)) %>%
  mutate(label_symbol = factor(label_symbol, levels = top100_symbols)) %>%
  tidyr::pivot_longer(
    cols = all_of(intersect(female_cols, colnames(female_top100))),
    names_to = "variance_component",
    values_to = "proportion"
  ) %>%
  mutate(
    variance_component = case_when(
      variance_component == "prop_age" ~ "Age",
      variance_component == "prop_idno" ~ "Individual",
      variance_component == "prop_residual" ~ "Residual",
      TRUE ~ variance_component
    ),
    variance_component = factor(
      variance_component,
      levels = c("Residual", "Individual", "Age")
    )
  )

# Plotting function (circular stacked barplot)
draw_stacked_barplot <- function(data_long, title_text, fill_colors, output_filename, source_data = NULL) {
  n_genes <- length(unique(data_long$label_symbol))
  
  # Sort by age contribution (prop_age)
  # Extract prop_age value for each gene from data_long
  gene_prop_age <- data_long %>%
    filter(variance_component == "Age") %>%
    select(label_symbol, proportion) %>%
    rename(prop_age = proportion) %>%
    group_by(label_symbol) %>%
    summarise(prop_age = first(prop_age), .groups = "drop") %>%
    arrange(desc(prop_age))  # Sort by age contribution from large to small
  
  # Use sorted gene order
  gene_levels <- gene_prop_age$label_symbol
  
  # Create index for each gene (for polar coordinates)
  data_long <- data_long %>%
    mutate(
      label_symbol = factor(label_symbol, levels = gene_levels),
      gene_index = as.numeric(label_symbol)
    )
  
  # Calculate total proportion for each gene (to determine y-axis range)
  max_prop <- data_long %>%
    group_by(gene_index) %>%
    summarise(total_prop = sum(proportion, na.rm = TRUE)) %>%
    pull(total_prop) %>%
    max(na.rm = TRUE)
  
  # Prepare label data (use ggh4x geom_text_aimed for radial alignment)
  label_data <- data.frame(
    gene_index = 1:n_genes,
    label = gene_levels,
    y_pos = max_prop * 1.15  # Label position outside barplot
  )
  
  # Create circular stacked barplot
  p <- ggplot(data_long, aes(x = gene_index, y = proportion, fill = variance_component)) +
    geom_bar(stat = "identity", position = position_stack(reverse = FALSE), width = 0.8) +
    # Use geom_text_aimed to add radially aligned labels
    # xend set to gene index (keep x unchanged), yend=0 means label points to center (y=0), achieving radial alignment
    geom_text_aimed(
      data = label_data,
      aes(x = gene_index, y = y_pos, label = label, xend = gene_index, yend = 0),
      inherit.aes = FALSE,
      size = 7,  # Font size (unit is mm, approximately 16-18pt)
      color = "black",
      hjust = 0.5,
      vjust = 0.5
    ) +
    scale_fill_manual(
      values = fill_colors,
      name = "Variance\nComponent"
    ) +
    scale_x_continuous(
      breaks = NULL,  # Don't show default labels
      expand = expansion(add = c(0.5, 0.5))  # Add spacing at both ends to avoid connecting first and last gene bars
    ) +
    scale_y_continuous(
      limits = c(0, max_prop * 1.25),  # Increase space to accommodate labels
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = "",
      y = "",
      title = title_text
    ) +
    coord_polar(theta = "x", start = 0, direction = 1) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),  # Hide default labels
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(
        color = "black", 
        size = 20, 
        face = "bold", 
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      legend.position = "right",
      legend.title = element_text(size = 25, face = "plain", color = "black"),
      legend.text = element_text(size = 25, color = "black"),
      legend.key.size = unit(1.5, "cm"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
  
  output_file <- file.path(output_dir, output_filename)
  # Square image
  ggsave(output_file, plot = p, width = 16, height = 16, dpi = 300, bg = "white", limitsize = FALSE)
  cat(sprintf("Saved: %s (width=16, height=16, square)\n", output_file))
  
  return(p)
}

# Define color scheme (referring to image: soft, modern color palette)
# Use soft light tones, similar to light yellow, light cyan, light purple in the image
colors_mixed <- c(
  "Age" = "#FFE5B4",      # Light yellow/cream - age contribution (referring to leftmost color in image)
  "Sex" = "#B8E6D1",      # Light cyan/mint green - sex contribution (referring to middle color in image)
  "Individual" = "#D4C5E8", # Light purple/lavender - individual contribution (referring to rightmost color in image)
  "Residual" = "#E8E8E8"   # Light gray - residual (neutral color)
)

colors_sex_separated <- c(
  "Age" = "#FFE5B4",      # Light yellow/cream - age contribution
  "Individual" = "#D4C5E8", # Light purple/lavender - individual contribution
  "Residual" = "#E8E8E8"   # Light gray - residual
)

# Plot three figures
cat("\nStarting to plot stacked barplots...\n")

p1 <- draw_stacked_barplot(
  mixed_long,
  "Mixed (Both Sexes): Variance Components",
  colors_mixed,
  "top100_variance_mixed_stacked_barplot.png"
)

p2 <- draw_stacked_barplot(
  male_long,
  "Male: Variance Components",
  colors_sex_separated,
  "top100_variance_male_stacked_barplot.png"
)

p3 <- draw_stacked_barplot(
  female_long,
  "Female: Variance Components",
  colors_sex_separated,
  "top100_variance_female_stacked_barplot.png"
)

# Print plots
print(p1)
print(p2)
print(p3)

# =============================================================================
# Step 6: Statistics on gene count and proportion under age contribution thresholds
# =============================================================================
cat("\nStep 6: Statistics on gene count and proportion under age contribution thresholds...\n")

# Function: Calculate statistics
calculate_age_contribution_stats <- function(df, top100_symbols, dataset_name) {
  # All genes statistics
  all_genes <- df %>%
    filter(!is.na(prop_age))
  
  n_all <- nrow(all_genes)
  n_all_gt02 <- sum(all_genes$prop_age > 0.2, na.rm = TRUE)
  n_all_gt01 <- sum(all_genes$prop_age > 0.1, na.rm = TRUE)
  n_all_gt005 <- sum(all_genes$prop_age > 0.05, na.rm = TRUE)
  
  prop_all_gt02 <- ifelse(n_all > 0, n_all_gt02 / n_all, 0)
  prop_all_gt01 <- ifelse(n_all > 0, n_all_gt01 / n_all, 0)
  prop_all_gt005 <- ifelse(n_all > 0, n_all_gt005 / n_all, 0)
  
  # Top100 genes statistics (use gene_symbol or label_symbol to match)
  if ("gene_symbol" %in% colnames(df)) {
    top100_genes <- df %>%
      filter(!is.na(prop_age) & gene_symbol %in% top100_symbols)
  } else if ("label_symbol" %in% colnames(df)) {
    top100_genes <- df %>%
      filter(!is.na(prop_age) & label_symbol %in% top100_symbols)
  } else {
    # If neither exists, return empty data frame
    top100_genes <- df %>%
      filter(FALSE)
  }
  
  n_top100 <- nrow(top100_genes)
  n_top100_gt02 <- sum(top100_genes$prop_age > 0.2, na.rm = TRUE)
  n_top100_gt01 <- sum(top100_genes$prop_age > 0.1, na.rm = TRUE)
  n_top100_gt005 <- sum(top100_genes$prop_age > 0.05, na.rm = TRUE)
  
  prop_top100_gt02 <- ifelse(n_top100 > 0, n_top100_gt02 / n_top100, 0)
  prop_top100_gt01 <- ifelse(n_top100 > 0, n_top100_gt01 / n_top100, 0)
  prop_top100_gt005 <- ifelse(n_top100 > 0, n_top100_gt005 / n_top100, 0)
  
  # Return results (according to user-required row name format)
  # -num rows fill count, -prop and -top100 rows fill proportion
  result <- data.frame(
    row.names = c(
      paste0(dataset_name, "-all-num"),
      paste0(dataset_name, "-all-prop"),
      paste0(dataset_name, "-top100")
    ),
    "GreaterThan0.2" = c(n_all_gt02, prop_all_gt02, prop_top100_gt02),
    "GreaterThan0.1" = c(n_all_gt01, prop_all_gt01, prop_top100_gt01),
    "GreaterThan0.05" = c(n_all_gt005, prop_all_gt005, prop_top100_gt005),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("\n%s statistics:\n", dataset_name))
  cat(sprintf("  All genes: total=%d, >0.2=%d (%.4f), >0.1=%d (%.4f), >0.05=%d (%.4f)\n",
              n_all, n_all_gt02, prop_all_gt02, n_all_gt01, prop_all_gt01, 
              n_all_gt005, prop_all_gt005))
  cat(sprintf("  Top100 genes: total=%d, >0.2=%d (%.4f), >0.1=%d (%.4f), >0.05=%d (%.4f)\n",
              n_top100, n_top100_gt02, prop_top100_gt02, n_top100_gt01, prop_top100_gt01,
              n_top100_gt005, prop_top100_gt005))
  
  return(result)
}

# Calculate statistics for each dataset
cat("\nCalculating mixed data statistics...\n")
stats_mixed <- calculate_age_contribution_stats(mixed_df, top100_symbols, "mixed")

cat("\nCalculating male data statistics...\n")
stats_male <- calculate_age_contribution_stats(male_df, top100_symbols, "male")

cat("\nCalculating female data statistics...\n")
stats_female <- calculate_age_contribution_stats(female_df, top100_symbols, "female")

# Merge all statistics
stats_combined <- rbind(stats_mixed, stats_male, stats_female)

# Save table
output_stats_file <- file.path(output_dir, "age_contribution_threshold_stats.csv")
write.csv(stats_combined, output_stats_file, row.names = TRUE)
cat(sprintf("\nStatistics table saved to: %s\n", output_stats_file))

# Print table
cat("\nStatistics table:\n")
print(stats_combined)

cat("\nAnalysis completed!\n")
cat("Output directory:", output_dir, "\n")

