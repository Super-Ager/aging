############################################################
## Protein Variance Decomposition Visualization
## Based on linear mixed-effects model results
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(ggrepel)
library(ggalluvial)

# Source function for protein trajectory plotting
source("005-plot-protein-trajectory.R")

options(warn = -1)

# ----------------------------------------------------------
# Paths
# ----------------------------------------------------------
df_combined_path <- "/home/ww/Project/Protein-Composition/result/1227-result/combined_protein_analysis_results.csv"
df_male_path <- "/home/ww/Project/Protein-Composition/result/1227-result/combined_protein_analysis_results_male.csv"
df_female_path <- "/home/ww/Project/Protein-Composition/result/1227-result/combined_protein_analysis_results_female.csv"
spearman_path <- "/home/ww/Project/Protein-Composition/result/1227-result-no-sex-separation/table-spearman/spearman_results.csv"
soma_path <- "/home/ww/Project/Longi_OA/results/tony_model/soma_df_Age_Sex_F_Dataset.csv"
anno_path <- "/home/ww/Project/AnzhenData/annoinfo_df_SOMA.csv"

out_root <- "/home/ww/Project/Protein-Composition/result/00-1-蛋白水平方差分解的统计"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_root, "analysis_log.txt")
legend_file <- file.path(out_root, "figure_legends.txt")
sink(log_file, append = FALSE)
# Initialize legend file
cat("Figure Legends\n", file = legend_file, append = FALSE)
cat("===============\n\n", file = legend_file, append = TRUE)

cat("========================================\n")
cat("Protein Variance Decomposition Analysis\n")
cat("========================================\n\n")

# ----------------------------------------------------------
# Load data
# ----------------------------------------------------------
cat("Loading data...\n")
df_combined <- read.csv(df_combined_path, stringsAsFactors = FALSE)
df_male <- read.csv(df_male_path, stringsAsFactors = FALSE)
df_female <- read.csv(df_female_path, stringsAsFactors = FALSE)
df_spearman <- read.csv(spearman_path, stringsAsFactors = FALSE)

# Load data for trajectory plotting
soma_df <- read.csv(soma_path, stringsAsFactors = FALSE) %>%
  filter(Dataset %in% 1:4) %>%
  rename(visit = Dataset)

anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
anno_map_symbol  <- setNames(anno_df$EntrezGeneSymbol, anno_df$AptName)
anno_map_uniprot <- setNames(anno_df$UniProt, anno_df$AptName)

num_categories <- length(unique(anno_df$EntrezGeneSymbol))
cat("Number of EntrezGeneSymbol categories:", num_categories, "\n")
cat(sprintf("df_combined: %d proteins\n", nrow(df_combined)))
cat(sprintf("df_male: %d proteins\n", nrow(df_male)))
cat(sprintf("df_female: %d proteins\n", nrow(df_female)))
cat(sprintf("df_spearman: %d proteins\n", nrow(df_spearman)))
cat("\n")

# ----------------------------------------------------------
# Color schemes
# ----------------------------------------------------------
colors_mixed <- c(
  "Age" = "#FFE5B4",
  "Sex" = "#B8E6D1",
  "Individual" = "#D4C5E8",
  "Residual" = "#E8E8E8"
)

colors_sex_separated <- c(
  "Age" = "#FFE5B4",
  "Individual" = "#D4C5E8",
  "Residual" = "#E8E8E8"
)

# Unified figure settings (except 01)
fig_height <- 5
# Font sizes: title=12, axis labels=10, annotations=8
font_title <- 12
font_axis <- 10
font_annotation <- 8

# ----------------------------------------------------------
# Create output directories
# ----------------------------------------------------------
out_combined <- file.path(out_root, "combined_analysis")
out_sex_separated <- file.path(out_root, "sex_separated_analysis")

dir.create(out_combined, recursive = TRUE, showWarnings = FALSE)
dir.create(out_sex_separated, recursive = TRUE, showWarnings = FALSE)













# ==========================================================
# Analysis based on df_combined
# ==========================================================
cat("========================================\n")
cat("Analysis based on df_combined\n")
cat("========================================\n\n")

# Task 1: Variance proportion distribution
# cat("Task 1: Variance proportion distribution plots...\n")
# df_variance <- df_combined %>%
#   select(prop_age, prop_sex, prop_idno, prop_residual) %>%
#   pivot_longer(cols = everything(), names_to = "Component", values_to = "Proportion") %>%
#   mutate(Component = case_when(
#     Component == "prop_age" ~ "Age",
#     Component == "prop_sex" ~ "Sex",
#     Component == "prop_idno" ~ "Individual",
#     Component == "prop_residual" ~ "Residual"
#   ),
#   Component = factor(Component, levels = c("Age", "Sex", "Individual", "Residual")))

# # Create separate plots with different y-axis limits
# df_variance_row1 <- df_variance %>% filter(Component %in% c("Age", "Sex"))
# df_variance_row2 <- df_variance %>% filter(Component %in% c("Individual", "Residual"))

# p1_row1 <- ggplot(df_variance_row1, aes(x = Proportion, fill = Component)) +
#   geom_histogram(bins = 50, alpha = 1, position = "identity") +
#   facet_wrap(~ Component, ncol = 2) +
#   scale_fill_manual(values = colors_mixed) +
#   labs(x = "Variance Proportion", y = "Count") +
#   theme_bw() +
#   theme(legend.position = "none") +
#   ylim(0, 10000)

# p1_row2 <- ggplot(df_variance_row2, aes(x = Proportion, fill = Component)) +
#   geom_histogram(bins = 50, alpha = 1, position = "identity") +
#   facet_wrap(~ Component, ncol = 2) +
#   scale_fill_manual(values = colors_mixed) +
#   labs(x = "Variance Proportion", y = "Count") +
#   theme_bw() +
#   theme(legend.position = "none") +
#   ylim(0, 1250)

# p1_combined <- ggarrange(p1_row1, p1_row2, nrow = 2, heights = c(1, 1))

# ggsave(file.path(out_combined, "01_variance_proportion_distribution.png"),
#        p1_combined, width = 7, height = 5.5, dpi = 300)

# # Task 1: Variance proportion distribution
cat("Task 1: Variance proportion distribution plots...\n")
df_variance <- df_combined %>%
  select(prop_age, prop_sex, prop_idno, prop_residual) %>%
  pivot_longer(cols = everything(), names_to = "Component", values_to = "Proportion") %>%
  mutate(
    Proportion = Proportion * 100,  # Convert to percentage
    Component = case_when(
      Component == "prop_age" ~ "Age",
      Component == "prop_sex" ~ "Sex",
      Component == "prop_idno" ~ "Individual",
      Component == "prop_residual" ~ "Residual"
    ),
    Component = factor(Component, levels = c("Age", "Sex", "Individual", "Residual"))
  )

# Define custom y-axis range for each Component (may need adjustment after converting to percentage)
y_limits <- list(
  Age = c(0, 5000),       # Age y-axis range
  Sex = c(0, 10000),        # Sex y-axis range
  Individual = c(0, 800), # Individual y-axis range
  Residual = c(0, 800)    # Residual y-axis range
)

# Create plots for each Component separately
p_age <- ggplot(df_variance %>% filter(Component == "Age"), aes(x = Proportion, fill = Component)) +
  geom_histogram(bins = 50, alpha = 1, position = "identity") +
  facet_wrap(~ Component) +
  scale_fill_manual(values = colors_mixed) +
  labs(x = "Variance Explained (%) by Age", y = "Number of proteins") +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(y_limits$Age) +
  scale_x_continuous(labels = scales::percent_format(scale = 1))  # Add percentage sign

p_sex <- ggplot(df_variance %>% filter(Component == "Sex"), aes(x = Proportion, fill = Component)) +
  geom_histogram(bins = 50, alpha = 1, position = "identity") +
  facet_wrap(~ Component) +
  scale_fill_manual(values = colors_mixed) +
  labs(x = "Variance Explained (%) by Sex", y = "Number of proteins") +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(y_limits$Sex) +
  scale_x_continuous(labels = scales::percent_format(scale = 1))

p_individual <- ggplot(df_variance %>% filter(Component == "Individual"), aes(x = Proportion, fill = Component)) +
  geom_histogram(bins = 50, alpha = 1, position = "identity") +
  facet_wrap(~ Component) +
  scale_fill_manual(values = colors_mixed) +
  labs(x = "Variance Explained (%) by Individual", y = "Number of proteins") +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(y_limits$Individual) +
  scale_x_continuous(labels = scales::percent_format(scale = 1))

p_residual <- ggplot(df_variance %>% filter(Component == "Residual"), aes(x = Proportion, fill = Component)) +
  geom_histogram(bins = 50, alpha = 1, position = "identity") +
  facet_wrap(~ Component) +
  scale_fill_manual(values = colors_mixed) +
  labs(x = "Variance Explained (%) by Residual", y = "Number of proteins") +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(y_limits$Residual) +
  scale_x_continuous(labels = scales::percent_format(scale = 1))

# Combine all plots (2 rows x 2 columns layout)
p1_combined <- ggarrange(
  p_age, p_sex, 
  p_individual, p_residual, 
  nrow = 2, ncol = 2, 
  heights = c(1, 1),
  widths = c(1, 1)
)

ggsave(file.path(out_combined, "01_variance_proportion_distribution.png"),
       p1_combined, width = 7, height = 5.5, dpi = 300)




# Task 2: Count proteins with maximum variance for each component
cat("Task 2: Counting proteins with maximum variance...\n")
df_combined <- df_combined %>%
  mutate(
    max_component = case_when(
      prop_age >= prop_sex & prop_age >= prop_idno & prop_age >= prop_residual ~ "Age",
      prop_sex >= prop_idno & prop_sex >= prop_residual ~ "Sex",
      prop_idno >= prop_residual ~ "Individual",
      TRUE ~ "Residual"
    )
  )

count_summary <- df_combined %>%
  count(max_component) %>%
  mutate(
    percentage = n / nrow(df_combined) * 100,
    label = sprintf("%d\n(%.1f%%)", n, percentage)
  )

p2 <- ggplot(count_summary, aes(x = max_component, y = n, fill = max_component)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = label), vjust = -0.5, size = font_annotation-4) +
  scale_fill_manual(values = colors_mixed) +
  labs(x = "", y = "Number of Proteins") +
  theme_bw() +
  theme(axis.text = element_text(size = font_axis-1, color = "black"),
        axis.title = element_text(size = font_title-1),
        legend.position = "none") +
  ylim(0, 10000)

ggsave(file.path(out_combined, "02_max_variance_count.png"),
       p2, width = 3.5, height = fig_height/2, dpi = 300)

# Write legend to file
cat("02_max_variance_count.png: Proteins with Maximum Variance by Component\n", file = legend_file, append = TRUE)

cat("\nProteins with maximum variance:\n")
print(count_summary)
cat("\n")


# # Task 3: Stacked barplot for proteins with Age as maximum variance
# cat("Task 3: Stacked barplot for Age-dominant proteins...\n")
df_age_max <- df_combined %>%
  filter(max_component == "Age") %>%
  arrange(desc(prop_age)) %>%
  select(seqid, gene_symbol, prop_age, prop_sex, prop_idno, prop_residual)

# df_age_max_long <- df_age_max %>%
#   pivot_longer(cols = c(prop_age, prop_sex, prop_idno, prop_residual),
#                names_to = "Component", values_to = "Proportion") %>%
#   mutate(Component = case_when(
#     Component == "prop_age" ~ "Age",
#     Component == "prop_sex" ~ "Sex",
#     Component == "prop_idno" ~ "Individual",
#     Component == "prop_residual" ~ "Residual"
#   ),
#   Component = factor(Component, levels = c("Residual", "Individual", "Sex", "Age"))) %>%
#   mutate(seqid = factor(seqid, levels = unique(df_age_max$seqid)))

# # Create two-line labels for x-axis
# age_labels <- paste0(df_age_max$gene_symbol, "\n(", df_age_max$seqid, ")")

# p3 <- ggplot(df_age_max_long, aes(x = seqid, y = Proportion, fill = Component)) +
#   geom_bar(stat = "identity", position = "stack", width = 0.8) +
#   scale_fill_manual(values = colors_mixed) +
#   scale_x_discrete(labels = age_labels) +
#   labs(x = "Protein", y = "Variance Proportion", fill = "Component") +
#   theme_bw() +
#   theme(axis.text = element_text(size = font_axis, color = "black"),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.title = element_text(size = font_title),
#         legend.text = element_text(size = font_axis),
#         legend.title = element_text(size = font_title),
#         legend.position = "bottom")

# ggsave(file.path(out_combined, "03_age_dominant_stacked.png"),
#        p3, width = nrow(df_age_max) * 0.7, height = fig_height, dpi = 300, limitsize = FALSE)

# # Write legend to file
# cat("03_age_dominant_stacked.png: Variance Proportions for Age-Dominant Proteins\n", file = legend_file, append = TRUE)


# Correctly implement horizontal stacked barplot
cat("Task 3 alternative: Stacked barplot for Age-dominant proteins (one per gene_symbol, top 20)...\n")
df_age_max_by_gene <- df_combined %>%
  filter(max_component == "Age") %>%
  arrange(desc(prop_age)) %>%
  group_by(gene_symbol) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(desc(prop_age)) %>%
  head(20) %>%
  select(seqid, gene_symbol, prop_age, prop_sex, prop_idno, prop_residual)

df_age_max_by_gene_long <- df_age_max_by_gene %>%
  pivot_longer(cols = c(prop_age, prop_sex, prop_idno, prop_residual),
               names_to = "Component", values_to = "Proportion") %>%
  mutate(
    Proportion = Proportion * 100,  # Convert to percentage
    Component = case_when(
      Component == "prop_age" ~ "Age",
      Component == "prop_sex" ~ "Sex",
      Component == "prop_idno" ~ "Individual",
      Component == "prop_residual" ~ "Residual"
    ),
    Component = factor(Component, levels = c("Residual", "Individual", "Sex", "Age"))) %>%
  # Reverse factor order so top1 gene appears at the top
  mutate(gene_symbol = factor(gene_symbol, levels = rev(unique(df_age_max_by_gene$gene_symbol))))

age_labels_by_gene <- df_age_max_by_gene$gene_symbol

# Correctly implement horizontal stacked barplot
p3_by_gene <- ggplot(df_age_max_by_gene_long, aes(x = Proportion, y = gene_symbol, fill = Component)) +
  geom_col(position = "stack", width = 0.8) +  # Horizontal stacking core
  scale_fill_manual(values = colors_mixed) +  # Use original color scheme
  scale_x_continuous(
    labels = scales::percent_format(scale = 1)  # Horizontal axis shows percentage
  ) +
  labs(
    x = "Variance Explained (%)",  # Modify x-axis label
    y = NULL,  # Don't show y-axis label
    fill = "Component"  # Legend title
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = font_axis, color = "black"),
    axis.text.y = element_text(hjust = 0),
    axis.title = element_text(size = font_title),
    # Place legend on the right
    legend.position = "right"
  )

# Adjust save size to fit horizontal barplot
ggsave(file.path(out_combined, "03_age_dominant_stacked_by_gene.png"),
       p3_by_gene, 
       width = 3.5 + 1,  # Increase width to accommodate legend
       height = nrow(df_age_max_by_gene) * 0.4,
       dpi = 300, 
       limitsize = FALSE)
cat(sprintf("Total Age-dominant proteins: %d\n", nrow(df_age_max)))
cat(sprintf("Unique Age-dominant genes (one per gene): %d\n", nrow(df_age_max_by_gene)))
cat("Top 10 Age-dominant proteins (seqid):\n")
print(head(df_age_max$seqid, 10))
cat("\n")

# # Task 3 alternative: One protein per gene_symbol (top 20 unique genes)
# cat("Task 3 alternative: Stacked barplot for Age-dominant proteins (one per gene_symbol, top 20)...\n")
# df_age_max_by_gene <- df_combined %>%
#   filter(max_component == "Age") %>%
#   arrange(desc(prop_age)) %>%
#   group_by(gene_symbol) %>%
#   slice_head(n = 1) %>%
#   ungroup() %>%
#   arrange(desc(prop_age)) %>%
#   head(20) %>%
#   select(seqid, gene_symbol, prop_age, prop_sex, prop_idno, prop_residual)

# df_age_max_by_gene_long <- df_age_max_by_gene %>%
#   pivot_longer(cols = c(prop_age, prop_sex, prop_idno, prop_residual),
#                names_to = "Component", values_to = "Proportion") %>%
#   mutate(Component = case_when(
#     Component == "prop_age" ~ "Age",
#     Component == "prop_sex" ~ "Sex",
#     Component == "prop_idno" ~ "Individual",
#     Component == "prop_residual" ~ "Residual"
#   ),
#   Component = factor(Component, levels = c("Residual", "Individual", "Sex", "Age"))) %>%
#   mutate(gene_symbol = factor(gene_symbol, levels = unique(df_age_max_by_gene$gene_symbol)))

# # age_labels_by_gene <- paste0(df_age_max_by_gene$gene_symbol, "\n(", df_age_max_by_gene$seqid, ")")
# age_labels_by_gene <- df_age_max_by_gene$gene_symbol

# p3_by_gene <- ggplot(df_age_max_by_gene_long, aes(x = gene_symbol, y = Proportion, fill = Component)) +
#   geom_bar(stat = "identity", position = "stack", width = 0.8) +
#   scale_fill_manual(values = colors_mixed) +
#   scale_x_discrete(labels = age_labels_by_gene) +
#   labs(x = "Protein", y = "Variance Proportion", fill = "Component") +
#   theme_bw() +
#   theme(axis.text = element_text(size = font_axis, color = "black"),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.title = element_text(size = font_title),
#         legend.text = element_text(size = font_axis),
#         legend.title = element_text(size = font_title),
#         legend.position = "bottom")

# ggsave(file.path(out_combined, "03_age_dominant_stacked_by_gene.png"),
#        p3_by_gene, width = nrow(df_age_max_by_gene) * 0.7, height = fig_height/2, dpi = 300, limitsize = FALSE)

# # Write legend to file
# cat("03_age_dominant_stacked_by_gene.png: Variance Proportions for Age-Dominant Proteins (Top Genes, One per Gene)\n", file = legend_file, append = TRUE)



# Task 4: Stacked barplot for Sex-dominant proteins (top 10)
cat("Task 4: Stacked barplot for Sex-dominant proteins (top 10)...\n")
df_sex_max <- df_combined %>%
  filter(max_component == "Sex") %>%
  arrange(desc(prop_sex)) %>%
  head(10) %>%
  select(seqid, gene_symbol, prop_age, prop_sex, prop_idno, prop_residual)

df_sex_max_long <- df_sex_max %>%
  pivot_longer(cols = c(prop_age, prop_sex, prop_idno, prop_residual),
               names_to = "Component", values_to = "Proportion") %>%
  mutate(Component = case_when(
    Component == "prop_age" ~ "Age",
    Component == "prop_sex" ~ "Sex",
    Component == "prop_idno" ~ "Individual",
    Component == "prop_residual" ~ "Residual"
  ),
  Component = factor(Component, levels = c("Residual", "Individual", "Age", "Sex"))) %>%
  mutate(seqid = factor(seqid, levels = unique(df_sex_max$seqid)))

# Create two-line labels for x-axis
sex_labels <- paste0(df_sex_max$gene_symbol, "\n(", df_sex_max$seqid, ")")

p4 <- ggplot(df_sex_max_long, aes(x = seqid, y = Proportion, fill = Component)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = colors_mixed) +
  scale_x_discrete(labels = sex_labels) +
  labs(x = "Protein", y = "Variance Proportion", fill = "Component") +
  theme_bw() +
  theme(axis.text = element_text(size = font_axis, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = font_title),
        legend.text = element_text(size = font_axis),
        legend.title = element_text(size = font_title),
        legend.position = "bottom")

ggsave(file.path(out_combined, "04_sex_dominant_stacked_top10.png"),
       p4, width = nrow(df_sex_max) * 0.7, height = fig_height, dpi = 300)

# Write legend to file
cat("04_sex_dominant_stacked_top10.png: Variance Proportions for Sex-Dominant Proteins (Top 10)\n", file = legend_file, append = TRUE)

# Task 4 alternative: One protein per gene_symbol (top 20 unique genes)
cat("Task 4 alternative: Stacked barplot for Sex-dominant proteins (one per gene_symbol, top 20)...\n")
df_sex_max_by_gene <- df_combined %>%
  filter(max_component == "Sex") %>%
  arrange(desc(prop_sex)) %>%
  group_by(gene_symbol) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(desc(prop_sex)) %>%
  head(20) %>%
  select(seqid, gene_symbol, prop_age, prop_sex, prop_idno, prop_residual)

df_sex_max_by_gene_long <- df_sex_max_by_gene %>%
  pivot_longer(cols = c(prop_age, prop_sex, prop_idno, prop_residual),
               names_to = "Component", values_to = "Proportion") %>%
  mutate(Component = case_when(
    Component == "prop_age" ~ "Age",
    Component == "prop_sex" ~ "Sex",
    Component == "prop_idno" ~ "Individual",
    Component == "prop_residual" ~ "Residual"
  ),
  Component = factor(Component, levels = c("Residual", "Individual", "Age", "Sex"))) %>%
  mutate(gene_symbol = factor(gene_symbol, levels = unique(df_sex_max_by_gene$gene_symbol)))

sex_labels_by_gene <- paste0(df_sex_max_by_gene$gene_symbol, "\n(", df_sex_max_by_gene$seqid, ")")

p4_by_gene <- ggplot(df_sex_max_by_gene_long, aes(x = gene_symbol, y = Proportion, fill = Component)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = colors_mixed) +
  scale_x_discrete(labels = sex_labels_by_gene) +
  labs(x = "Protein", y = "Variance Proportion", fill = "Component") +
  theme_bw() +
  theme(axis.text = element_text(size = font_axis, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = font_title),
        legend.text = element_text(size = font_axis),
        legend.title = element_text(size = font_title),
        legend.position = "bottom")

ggsave(file.path(out_combined, "04_sex_dominant_stacked_by_gene.png"),
       p4_by_gene, width = nrow(df_sex_max_by_gene) * 0.7, height = fig_height, dpi = 300, limitsize = FALSE)

# Write legend to file
cat("04_sex_dominant_stacked_by_gene.png: Variance Proportions for Sex-Dominant Proteins (Top Genes, One per Gene)\n", file = legend_file, append = TRUE)

cat("Top 10 Sex-dominant proteins (seqid):\n")
print(df_sex_max$seqid)
cat(sprintf("Unique Sex-dominant genes (one per gene): %d\n", nrow(df_sex_max_by_gene)))
cat("\n")

# Task 5: Stacked barplot for Individual-dominant proteins (top 10)
cat("Task 5: Stacked barplot for Individual-dominant proteins (top 10)...\n")
df_ind_max <- df_combined %>%
  filter(max_component == "Individual") %>%
  arrange(desc(prop_idno)) %>%
  head(10) %>%
  select(seqid, gene_symbol, prop_age, prop_sex, prop_idno, prop_residual)

df_ind_max_long <- df_ind_max %>%
  pivot_longer(cols = c(prop_age, prop_sex, prop_idno, prop_residual),
               names_to = "Component", values_to = "Proportion") %>%
  mutate(Component = case_when(
    Component == "prop_age" ~ "Age",
    Component == "prop_sex" ~ "Sex",
    Component == "prop_idno" ~ "Individual",
    Component == "prop_residual" ~ "Residual"
  ),
  Component = factor(Component, levels = c("Residual", "Age", "Sex", "Individual"))) %>%
  mutate(seqid = factor(seqid, levels = unique(df_ind_max$seqid)))

# Create two-line labels for x-axis
ind_labels <- paste0(df_ind_max$gene_symbol, "\n(", df_ind_max$seqid, ")")

p5 <- ggplot(df_ind_max_long, aes(x = seqid, y = Proportion, fill = Component)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = colors_mixed) +
  scale_x_discrete(labels = ind_labels) +
  labs(x = "Protein", y = "Variance Proportion", fill = "Component") +
  theme_bw() +
  theme(axis.text = element_text(size = font_axis, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = font_title),
        legend.text = element_text(size = font_axis),
        legend.title = element_text(size = font_title),
        legend.position = "bottom")

ggsave(file.path(out_combined, "05_individual_dominant_stacked_top10.png"),
       p5, width = nrow(df_ind_max) * 0.7, height = fig_height, dpi = 300)

# Write legend to file
cat("05_individual_dominant_stacked_top10.png: Variance Proportions for Individual-Dominant Proteins (Top 10)\n", file = legend_file, append = TRUE)









# Task 5 alternative: One protein per gene_symbol (top 20 unique genes)
cat("Task 5 alternative: Stacked barplot for Individual-dominant proteins (one per gene_symbol, top 20)...\n")
df_ind_max_by_gene <- df_combined %>%
  filter(max_component == "Individual") %>%
  arrange(desc(prop_idno)) %>%
  group_by(gene_symbol) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(desc(prop_idno)) %>%
  head(20) %>%
  select(seqid, gene_symbol, prop_age, prop_sex, prop_idno, prop_residual)

df_ind_max_by_gene_long <- df_ind_max_by_gene %>%
  pivot_longer(cols = c(prop_age, prop_sex, prop_idno, prop_residual),
               names_to = "Component", values_to = "Proportion") %>%
  mutate(Component = case_when(
    Component == "prop_age" ~ "Age",
    Component == "prop_sex" ~ "Sex",
    Component == "prop_idno" ~ "Individual",
    Component == "prop_residual" ~ "Residual"
  ),
  Component = factor(Component, levels = c("Residual", "Age", "Sex", "Individual"))) %>%
  mutate(gene_symbol = factor(gene_symbol, levels = unique(df_ind_max_by_gene$gene_symbol)))

ind_labels_by_gene <- paste0(df_ind_max_by_gene$gene_symbol, "\n(", df_ind_max_by_gene$seqid, ")")

p5_by_gene <- ggplot(df_ind_max_by_gene_long, aes(x = gene_symbol, y = Proportion, fill = Component)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = colors_mixed) +
  scale_x_discrete(labels = ind_labels_by_gene) +
  labs(x = "Protein", y = "Variance Proportion", fill = "Component") +
  theme_bw() +
  theme(axis.text = element_text(size = font_axis, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = font_title),
        legend.text = element_text(size = font_axis),
        legend.title = element_text(size = font_title),
        legend.position = "bottom")

ggsave(file.path(out_combined, "05_individual_dominant_stacked_by_gene.png"),
       p5_by_gene, width = nrow(df_ind_max_by_gene) * 0.7, height = fig_height, dpi = 300, limitsize = FALSE)

# Write legend to file
cat("05_individual_dominant_stacked_by_gene.png: Variance Proportions for Individual-Dominant Proteins (Top Genes, One per Gene)\n", file = legend_file, append = TRUE)

# Task: Scatter plot for residual variance proportion vs SCC (r_23)
cat("Task: Scatter plot for residual variance proportion vs SCC...\n")
df_residual_scc <- df_combined %>%
  select(seqid, prop_residual) %>%
  inner_join(
    df_spearman %>% select(seqid, r_23),
    by = "seqid"
  ) %>%
  filter(!is.na(prop_residual) & !is.na(r_23))

if (nrow(df_residual_scc) > 0) {

  cor_residual_scc <- cor.test(
    df_residual_scc$prop_residual,
    df_residual_scc$r_23
  )

  pval_text <- ifelse(
    cor_residual_scc$p.value < 0.001,
    "p < 0.001",
    sprintf("p = %.3f", cor_residual_scc$p.value)
  )

  label_text <- sprintf(
    "r = %.2f\n%s",
    cor_residual_scc$estimate,
    pval_text
  )

  # data-driven annotation position
  x_anno <- quantile(df_residual_scc$prop_residual, 1, na.rm = TRUE)
  y_anno <- quantile(df_residual_scc$r_23, 1, na.rm = TRUE)

  p_residual_scc <- ggplot(
    df_residual_scc,
    aes(x = prop_residual, y = r_23)
  ) +
    geom_point(
      size = 1,
      alpha = 0.7,
      color = "grey30"
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = "black",
      linewidth = 0.6,
      alpha = 0.15
    ) +
    annotate(
      "text",
      x = x_anno,
      y = y_anno,
      label = label_text,
      hjust = 1,
      vjust = 1,
      size = font_annotation/2,
      color = "black"
    ) +
    labs(
      x = "Residual variance proportion",
      y = "Spearman R"
    ) +
    theme_classic() +
    theme(
      axis.text  = element_text(size = font_axis, color = "black"),
      axis.title = element_text(size = font_title, color = "black"),
      axis.line  = element_line(color = "black", linewidth = 0.6),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )

  ggsave(
    file.path(out_combined, "residual_prop_vs_scc.png"),
    p_residual_scc,
    width = 3.5,
    height = fig_height/2,
    dpi = 300
  )

  cat(
    "residual_prop_vs_scc.png: Residual variance proportion vs Spearman correlation (r_23)\n",
    file = legend_file,
    append = TRUE
  )

  cat(
    sprintf(
      "Correlation: r = %.3f, p = %.2e\n",
      cor_residual_scc$estimate,
      cor_residual_scc$p.value
    )
  )

} else {
  cat("No data for residual variance proportion vs SCC\n")
  p_residual_scc <- NULL
}















# Task: Plot protein trajectory for seq.3309.2
cat("Task: Plotting protein trajectory for seq.3309.2...\n")
target_seqid <- "seq.3309.2"
p_trajectory <- plot_protein_trajectory_sankey(
  seqid = target_seqid,
  soma_df = soma_df,
  anno_map_symbol = anno_map_symbol,
  anno_map_uniprot = anno_map_uniprot,
  font_title = font_title-2,
  font_axis = font_axis-2,
  font_annotation = font_annotation-5
)

if (is.null(p_trajectory)) {
  cat(sprintf("Failed to plot trajectory for %s\n", target_seqid))
  p_trajectory <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = "No data available", size = 10) +
    theme_void()
}

ggsave(file.path(out_combined, "example-trajectory.png"),
       p_trajectory, width = 3.5, height = fig_height/2, dpi = 300, limitsize = FALSE)

# # Create combined figure for paper
# cat("Creating combined figure for paper...\n")
# if (!is.null(p_residual_scc)) {
#   p_combined_paper <- ggarrange(
#     ggarrange(p1_combined, p2, ncol = 2, widths = c(1.2, 1), labels = c("A", "B")),
#     ggarrange(p3_by_gene, p4_by_gene, ncol = 2, widths = c(1, 1), labels = c("C", "D")),
#     ggarrange(p5_by_gene, NULL, ncol = 2, widths = c(1, 0), labels = c("E", "")),
#     ggarrange(p_residual_scc, p_trajectory, ncol = 2, widths = c(1, 1), labels = c("F", "G")),
#     nrow = 4, heights = c(1, 1, 1, 1)
#   )
# } else {
#   p_combined_paper <- ggarrange(
#     ggarrange(p1_combined, p2, ncol = 2, widths = c(1.2, 1), labels = c("A", "B")),
#     ggarrange(p3_by_gene, p4_by_gene, ncol = 2, widths = c(1, 1), labels = c("C", "D")),
#     ggarrange(p5_by_gene, NULL, ncol = 2, widths = c(1, 0), labels = c("E", "")),
#     ggarrange(NULL, p_trajectory, ncol = 2, widths = c(0, 1), labels = c("", "F")),
#     nrow = 4, heights = c(1, 1, 1, 1)
#   )
# }

# ggsave(file.path(out_combined, "combined_figure_paper.png"),
#        p_combined_paper, width = 16, height = 20, dpi = 600, limitsize = FALSE)

# cat("Top 10 Individual-dominant proteins (seqid):\n")
# print(df_ind_max$seqid)
# cat(sprintf("Unique Individual-dominant genes (one per gene): %d\n", nrow(df_ind_max_by_gene)))
# cat("\n")












# Task 6: Stacked barplot for Residual-dominant proteins (top 10)
cat("Task 6: Stacked barplot for Residual-dominant proteins (top 10)...\n")
df_res_max <- df_combined %>%
  filter(max_component == "Residual") %>%
  arrange(desc(prop_residual)) %>%
  head(10) %>%
  select(seqid, gene_symbol, prop_age, prop_sex, prop_idno, prop_residual)

df_res_max_long <- df_res_max %>%
  pivot_longer(cols = c(prop_age, prop_sex, prop_idno, prop_residual),
               names_to = "Component", values_to = "Proportion") %>%
  mutate(Component = case_when(
    Component == "prop_age" ~ "Age",
    Component == "prop_sex" ~ "Sex",
    Component == "prop_idno" ~ "Individual",
    Component == "prop_residual" ~ "Residual"
  ),
  Component = factor(Component, levels = c("Age", "Sex", "Individual", "Residual"))) %>%
  mutate(seqid = factor(seqid, levels = unique(df_res_max$seqid)))

# Create two-line labels for x-axis
res_labels <- paste0(df_res_max$gene_symbol, "\n(", df_res_max$seqid, ")")

p6 <- ggplot(df_res_max_long, aes(x = seqid, y = Proportion, fill = Component)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = colors_mixed) +
  scale_x_discrete(labels = res_labels) +
  labs(x = "Protein", y = "Variance Proportion", fill = "Component") +
  theme_bw() +
  theme(axis.text = element_text(size = font_axis, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = font_title),
        legend.text = element_text(size = font_axis),
        legend.title = element_text(size = font_title),
        legend.position = "bottom")

ggsave(file.path(out_combined, "06_residual_dominant_stacked_top10.png"),
       p6, width = nrow(df_res_max) * 0.7, height = fig_height, dpi = 300)

# Write legend to file
cat("06_residual_dominant_stacked_top10.png: Variance Proportions for Residual-Dominant Proteins (Top 10)\n", file = legend_file, append = TRUE)

# Task 6 alternative: One protein per gene_symbol (top 20 unique genes)
cat("Task 6 alternative: Stacked barplot for Residual-dominant proteins (one per gene_symbol, top 20)...\n")
df_res_max_by_gene <- df_combined %>%
  filter(max_component == "Residual") %>%
  arrange(desc(prop_residual)) %>%
  group_by(gene_symbol) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(desc(prop_residual)) %>%
  head(20) %>%
  select(seqid, gene_symbol, prop_age, prop_sex, prop_idno, prop_residual)

df_res_max_by_gene_long <- df_res_max_by_gene %>%
  pivot_longer(cols = c(prop_age, prop_sex, prop_idno, prop_residual),
               names_to = "Component", values_to = "Proportion") %>%
  mutate(Component = case_when(
    Component == "prop_age" ~ "Age",
    Component == "prop_sex" ~ "Sex",
    Component == "prop_idno" ~ "Individual",
    Component == "prop_residual" ~ "Residual"
  ),
  Component = factor(Component, levels = c("Age", "Sex", "Individual", "Residual"))) %>%
  mutate(gene_symbol = factor(gene_symbol, levels = unique(df_res_max_by_gene$gene_symbol)))

res_labels_by_gene <- paste0(df_res_max_by_gene$gene_symbol, "\n(", df_res_max_by_gene$seqid, ")")

p6_by_gene <- ggplot(df_res_max_by_gene_long, aes(x = gene_symbol, y = Proportion, fill = Component)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = colors_mixed) +
  scale_x_discrete(labels = res_labels_by_gene) +
  labs(x = "Protein", y = "Variance Proportion", fill = "Component") +
  theme_bw() +
  theme(axis.text = element_text(size = font_axis, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = font_title),
        legend.text = element_text(size = font_axis),
        legend.title = element_text(size = font_title),
        legend.position = "bottom")

ggsave(file.path(out_combined, "06_residual_dominant_stacked_by_gene.png"),
       p6_by_gene, width = nrow(df_res_max_by_gene) * 0.7, height = fig_height, dpi = 300, limitsize = FALSE)

# Write legend to file
cat("06_residual_dominant_stacked_by_gene.png: Variance Proportions for Residual-Dominant Proteins (Top Genes, One per Gene)\n", file = legend_file, append = TRUE)

cat("Top 10 Residual-dominant proteins (seqid):\n")
print(df_res_max$seqid)
cat(sprintf("Unique Residual-dominant genes (one per gene): %d\n", nrow(df_res_max_by_gene)))
cat("\n")

# Task 7: Output results to log file
cat("\n========================================\n")
cat("Summary Statistics (Tasks 2-6)\n")
cat("========================================\n")
cat("\nProteins with maximum variance by component:\n")
print(count_summary)
cat("\nTop 10 Age-dominant proteins:\n")
print(head(df_age_max %>% select(seqid, gene_symbol), 10))
cat("\nTop 10 Sex-dominant proteins:\n")
print(df_sex_max %>% select(seqid, gene_symbol))
cat("\nTop 10 Individual-dominant proteins:\n")
print(df_ind_max %>% select(seqid, gene_symbol))
cat("\nTop 10 Residual-dominant proteins:\n")
print(df_res_max %>% select(seqid, gene_symbol))
cat("\n")

# Task 8: Volcano plot for age_Estimate and age_pvalue
cat("Task 8: Volcano plot for age effects...\n")
df_combined <- df_combined %>%
  mutate(
    age_log100p = -log(age_pvalue) / log(100),
    age_sig = age_pvalue < 0.01
  )

# Get top 10 positive and negative age_Estimate
top10_age_pos <- df_combined %>%
  filter(age_Estimate > 0) %>%
  arrange(desc(abs(age_Estimate))) %>%
  head(10) %>%
  mutate(label = ifelse(duplicated(gene_symbol) | duplicated(gene_symbol, fromLast = TRUE),
                        paste0(gene_symbol, " (", seqid, ")"), gene_symbol))

top10_age_neg <- df_combined %>%
  filter(age_Estimate < 0) %>%
  arrange(desc(abs(age_Estimate))) %>%
  head(10) %>%
  mutate(label = ifelse(duplicated(gene_symbol) | duplicated(gene_symbol, fromLast = TRUE),
                        paste0(gene_symbol, " (", seqid, ")"), gene_symbol))

top10_age_all <- bind_rows(top10_age_pos, top10_age_neg)

p8 <- ggplot(df_combined, aes(x = age_Estimate, y = age_log100p, color = age_sig)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_hline(yintercept = -log(0.01) / log(100), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_text_repel(data = top10_age_all, aes(label = label),
                  size = font_annotation - 4, color = "black", max.overlaps = Inf,
                  box.padding = 1, point.padding = 0.1) +
  scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "red2"),
                     labels = c("FALSE" = "p >= 0.01", "TRUE" = "p < 0.01"),
                     name = "Significance") +
  labs(x = "Age Estimate", y = "-log100(p-value)") +
  theme_bw() +
  theme(axis.text = element_text(size = font_axis, color = "black"),
        axis.title = element_text(size = font_title),
        legend.text = element_text(size = font_axis),
        legend.title = element_text(size = font_title),
        legend.position = "bottom")

ggsave(file.path(out_combined, "08_age_volcano.png"),
       p8, width = 8, height = fig_height, dpi = 300)

# Write legend to file
cat("08_age_volcano.png: Volcano Plot: Age Effects\n", file = legend_file, append = TRUE)

# Task 9: Volcano plot for sex_Estimate and sex_pvalue
cat("Task 9: Volcano plot for sex effects...\n")
df_combined <- df_combined %>%
  mutate(
    sex_log10p = -log10(sex_pvalue),
    sex_sig = sex_pvalue < 0.01
  )

# Get top 10 positive and negative sex_Estimate
top10_sex_pos <- df_combined %>%
  filter(sex_Estimate > 0) %>%
  arrange(desc(abs(sex_Estimate))) %>%
  head(10) %>%
  mutate(label = ifelse(duplicated(gene_symbol) | duplicated(gene_symbol, fromLast = TRUE),
                        paste0(gene_symbol, " (", seqid, ")"), gene_symbol))

top10_sex_neg <- df_combined %>%
  filter(sex_Estimate < 0) %>%
  arrange(desc(abs(sex_Estimate))) %>%
  head(10) %>%
  mutate(label = ifelse(duplicated(gene_symbol) | duplicated(gene_symbol, fromLast = TRUE),
                        paste0(gene_symbol, " (", seqid, ")"), gene_symbol))

top10_sex_all <- bind_rows(top10_sex_pos, top10_sex_neg)

p9 <- ggplot(df_combined, aes(x = sex_Estimate, y = sex_log10p, color = sex_sig)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_text_repel(data = top10_sex_all, aes(label = label),
                  size = font_annotation - 4, color = "black", max.overlaps = Inf,
                  box.padding = 0.5, point.padding = 0.3) +
  scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "red2"),
                     labels = c("FALSE" = "p >= 0.01", "TRUE" = "p < 0.01"),
                     name = "Significance") +
  labs(x = "Sex Estimate", y = "-log10(p-value)") +
  theme_bw() +
  theme(axis.text = element_text(size = font_axis, color = "black"),
        axis.title = element_text(size = font_title),
        legend.text = element_text(size = font_axis),
        legend.title = element_text(size = font_title),
        legend.position = "bottom")

ggsave(file.path(out_combined, "09_sex_volcano.png"),
       p9, width = 8, height = fig_height, dpi = 300)

# Write legend to file
cat("09_sex_volcano.png: Volcano Plot: Sex Effects\n", file = legend_file, append = TRUE)

# Task 10: Scatter plot for age_Estimate vs prop_age (p < 0.01)
cat("Task 10: Scatter plot for age_Estimate vs prop_age (p < 0.01)...\n")
df_age_sig <- df_combined %>%
  filter(age_pvalue < 0.01) %>%
  mutate(age_Estimate_abs = abs(age_Estimate))

if (nrow(df_age_sig) > 0) {
  cor_age <- cor.test(df_age_sig$age_Estimate_abs, df_age_sig$prop_age)
  lm_age <- lm(prop_age ~ age_Estimate_abs, data = df_age_sig)
  r2_age <- summary(lm_age)$r.squared
  pval_age <- cor_age$p.value
  
  pval_text <- ifelse(pval_age < 0.001, "p < 0.001", sprintf("p = %.3f", pval_age))
  label_text <- sprintf("r = %.3f\nR² = %.3f\n%s",
                        cor_age$estimate, r2_age, pval_text)
  
  p10 <- ggplot(df_age_sig, aes(x = age_Estimate_abs, y = prop_age)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    annotate("text", x = Inf, y = Inf, label = label_text,
             hjust = 1.1, vjust = 1.5, size = font_annotation) +
    labs(x = "|Age Estimate|", y = "Age Variance Proportion") +
    theme_bw() +
    theme(axis.text = element_text(size = font_axis, color = "black"),
          axis.title = element_text(size = font_title))
  
  ggsave(file.path(out_combined, "10_age_estimate_vs_prop_age.png"),
         p10, width = 6, height = fig_height, dpi = 300)
  
  # Write legend to file
  cat("10_age_estimate_vs_prop_age.png: Age Estimate vs Age Variance Proportion (p < 0.01)\n", file = legend_file, append = TRUE)
  
  cat(sprintf("Proteins with age_pvalue < 0.01: %d\n", nrow(df_age_sig)))
  cat(sprintf("Minimum age variance proportion: %.6f\n", min(df_age_sig$prop_age, na.rm = TRUE)))
  cat(sprintf("Correlation: r = %.3f, p = %.2e\n", cor_age$estimate, pval_age))
} else {
  cat("No proteins with age_pvalue < 0.01\n")
}
cat("\n")

# Task 11: Scatter plot for sex_Estimate vs prop_sex (p < 0.01)
cat("Task 11: Scatter plot for sex_Estimate vs prop_sex (p < 0.01)...\n")
df_sex_sig <- df_combined %>%
  filter(sex_pvalue < 0.01)

if (nrow(df_sex_sig) > 0) {
  cor_sex <- cor.test(df_sex_sig$sex_Estimate, df_sex_sig$prop_sex)
  
  pval_text <- ifelse(cor_sex$p.value < 0.001, "p < 0.001", sprintf("p = %.3f", cor_sex$p.value))
  label_text <- sprintf("r = %.3f\n%s",
                        cor_sex$estimate, pval_text)
  
  p11 <- ggplot(df_sex_sig, aes(x = sex_Estimate, y = prop_sex)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    annotate("text", x = Inf, y = Inf, label = label_text,
             hjust = 1.1, vjust = 1.5, size = font_annotation) +
    labs(x = "Sex Estimate", y = "Sex Variance Proportion") +
    theme_bw() +
    theme(axis.text = element_text(size = font_axis, color = "black"),
          axis.title = element_text(size = font_title))
  
  ggsave(file.path(out_combined, "11_sex_estimate_vs_prop_sex.png"),
         p11, width = 6, height = fig_height, dpi = 300)
  
  # Write legend to file
  cat("11_sex_estimate_vs_prop_sex.png: Sex Estimate vs Sex Variance Proportion (p < 0.01)\n", file = legend_file, append = TRUE)
  
  cat(sprintf("Proteins with sex_pvalue < 0.01: %d\n", nrow(df_sex_sig)))
  cat(sprintf("Minimum sex variance proportion: %.6f\n", min(df_sex_sig$prop_sex, na.rm = TRUE)))
  cat(sprintf("Correlation: r = %.3f, p = %.2e\n", cor_sex$estimate, cor_sex$p.value))
} else {
  cat("No proteins with sex_pvalue < 0.01\n")
}
cat("\n")



# ==========================================================
# Analysis based on df_male and df_female
# ==========================================================
cat("\n========================================\n")
cat("Analysis based on df_male and df_female\n")
cat("========================================\n\n")

# Merge male and female data
df_sex_merged <- df_male %>%
  select(seqid, prop_age_male = prop_age, prop_idno_male = prop_idno, prop_residual_male = prop_residual) %>%
  inner_join(
    df_female %>%
      select(seqid, prop_age_female = prop_age, prop_idno_female = prop_idno, prop_residual_female = prop_residual),
    by = "seqid"
  ) %>%
  inner_join(
    df_combined %>%
      filter(age_pvalue < 0.01) %>%
      select(seqid),
    by = "seqid"
  )

cat(sprintf("Proteins with age_pvalue < 0.01 in combined data: %d\n", nrow(df_sex_merged)))
cat("\n")

# Task 1: Scatter plot for age variance proportion (male vs female)
cat("Task 1: Scatter plot for age variance proportion (male vs female)...\n")
top10_age_diff <- NULL
if (nrow(df_sex_merged) > 0) {
  cor_age_sex <- cor.test(df_sex_merged$prop_age_male, df_sex_merged$prop_age_female)
  
  # Calculate differences
  df_sex_merged <- df_sex_merged %>%
    mutate(age_diff = abs(prop_age_male - prop_age_female))
  
  top10_age_diff <- df_sex_merged %>%
    arrange(desc(age_diff)) %>%
    head(10) %>%
    left_join(df_combined %>% select(seqid, gene_symbol), by = "seqid")
  
  pval_text <- ifelse(cor_age_sex$p.value < 0.001, "p < 0.001", sprintf("p = %.3f", cor_age_sex$p.value))
  label_text <- sprintf("r = %.3f\n%s",
                        cor_age_sex$estimate, pval_text)
  
  p12 <- ggplot(df_sex_merged, aes(x = log(prop_age_male + 1), y = log(prop_age_female + 1))) +
    geom_point(alpha = 0.6, color = "steelblue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    geom_smooth(method = "lm", se = TRUE, color = "red2") +
    annotate("text", x = Inf, y = Inf, label = label_text,
             hjust = 1.1, vjust = 1.5, size = font_annotation) +
    labs(x = "log(Age Variance Proportion + 1) (Male)", y = "log(Age Variance Proportion + 1) (Female)") +
    theme_bw() +
    theme(axis.text = element_text(size = font_axis, color = "black"),
          axis.title = element_text(size = font_title))
  
  ggsave(file.path(out_sex_separated, "01_age_prop_male_vs_female.png"),
         p12, width = 6, height = fig_height, dpi = 300)
  
  # Write legend to file
  cat("01_age_prop_male_vs_female.png: Age Variance Proportion: Male vs Female\n", file = legend_file, append = TRUE)
  
  cat(sprintf("Correlation: r = %.3f, p = %.2e\n", cor_age_sex$estimate, cor_age_sex$p.value))
  cat("Top 10 proteins with largest age variance proportion difference:\n")
  print(top10_age_diff %>% select(seqid, gene_symbol, age_diff))
  cat("\n")
}

# Task 2: Scatter plot for individual variance proportion (male vs female)
cat("Task 2: Scatter plot for individual variance proportion (male vs female)...\n")
top10_idno_diff <- NULL
if (nrow(df_sex_merged) > 0) {
  cor_idno_sex <- cor.test(df_sex_merged$prop_idno_male, df_sex_merged$prop_idno_female)
  
  # Calculate differences
  df_sex_merged <- df_sex_merged %>%
    mutate(idno_diff = abs(prop_idno_male - prop_idno_female))
  
  top10_idno_diff <- df_sex_merged %>%
    arrange(desc(idno_diff)) %>%
    head(10) %>%
    left_join(df_combined %>% select(seqid, gene_symbol), by = "seqid")
  
  pval_text <- ifelse(cor_idno_sex$p.value < 0.001, "p < 0.001", sprintf("p = %.3f", cor_idno_sex$p.value))
  label_text <- sprintf("r = %.3f\n%s",
                        cor_idno_sex$estimate, pval_text)
  
  p13 <- ggplot(df_sex_merged, aes(x = prop_idno_male, y = prop_idno_female)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    geom_smooth(method = "lm", se = TRUE, color = "red2") +
    annotate("text", x = Inf, y = Inf, label = label_text,
             hjust = 1.1, vjust = 1.5, size = font_annotation) +
    labs(x = "Individual Variance Proportion (Male)", y = "Individual Variance Proportion (Female)") +
    theme_bw() +
    theme(axis.text = element_text(size = font_axis, color = "black"),
          axis.title = element_text(size = font_title))
  
  ggsave(file.path(out_sex_separated, "02_individual_prop_male_vs_female.png"),
         p13, width = 6, height = fig_height, dpi = 300)
  
  # Write legend to file
  cat("02_individual_prop_male_vs_female.png: Individual Variance Proportion: Male vs Female\n", file = legend_file, append = TRUE)
  
  cat(sprintf("Correlation: r = %.3f, p = %.2e\n", cor_idno_sex$estimate, cor_idno_sex$p.value))
  cat("Top 10 proteins with largest individual variance proportion difference:\n")
  print(top10_idno_diff %>% select(seqid, gene_symbol, idno_diff))
  cat("\n")
}

# Task 3: Scatter plot for residual variance proportion (male vs female)
cat("Task 3: Scatter plot for residual variance proportion (male vs female)...\n")
top10_residual_diff <- NULL
if (nrow(df_sex_merged) > 0) {
  cor_residual_sex <- cor.test(df_sex_merged$prop_residual_male, df_sex_merged$prop_residual_female)
  
  # Calculate differences
  df_sex_merged <- df_sex_merged %>%
    mutate(residual_diff = abs(prop_residual_male - prop_residual_female))
  
  top10_residual_diff <- df_sex_merged %>%
    arrange(desc(residual_diff)) %>%
    head(10) %>%
    left_join(df_combined %>% select(seqid, gene_symbol), by = "seqid")
  
  pval_text <- ifelse(cor_residual_sex$p.value < 0.001, "p < 0.001", sprintf("p = %.3f", cor_residual_sex$p.value))
  label_text <- sprintf("r = %.3f\n%s",
                        cor_residual_sex$estimate, pval_text)
  
  p14 <- ggplot(df_sex_merged, aes(x = prop_residual_male, y = prop_residual_female)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    geom_smooth(method = "lm", se = TRUE, color = "red2") +
    annotate("text", x = Inf, y = Inf, label = label_text,
             hjust = 1.1, vjust = 1.5, size = font_annotation) +
    labs(x = "Residual Variance Proportion (Male)", y = "Residual Variance Proportion (Female)") +
    theme_bw() +
    theme(axis.text = element_text(size = font_axis, color = "black"),
          axis.title = element_text(size = font_title))
  
  ggsave(file.path(out_sex_separated, "03_residual_prop_male_vs_female.png"),
         p14, width = 6, height = fig_height, dpi = 300)
  
  # Write legend to file
  cat("03_residual_prop_male_vs_female.png: Residual Variance Proportion: Male vs Female\n", file = legend_file, append = TRUE)
  
  cat(sprintf("Correlation: r = %.3f, p = %.2e\n", cor_residual_sex$estimate, cor_residual_sex$p.value))
  cat("Top 10 proteins with largest residual variance proportion difference:\n")
  print(top10_residual_diff %>% select(seqid, gene_symbol, residual_diff))
  cat("\n")
}

# Output sex-separated results to log
cat("\n========================================\n")
cat("Sex-Separated Analysis Summary\n")
cat("========================================\n")
if (!is.null(top10_age_diff) && nrow(top10_age_diff) > 0) {
  cat("\nTop 10 proteins with largest age variance proportion difference:\n")
  print(top10_age_diff %>% select(seqid, gene_symbol, age_diff))
}
if (!is.null(top10_idno_diff) && nrow(top10_idno_diff) > 0) {
  cat("\nTop 10 proteins with largest individual variance proportion difference:\n")
  print(top10_idno_diff %>% select(seqid, gene_symbol, idno_diff))
}
if (!is.null(top10_residual_diff) && nrow(top10_residual_diff) > 0) {
  cat("\nTop 10 proteins with largest residual variance proportion difference:\n")
  print(top10_residual_diff %>% select(seqid, gene_symbol, residual_diff))
}

cat("\n========================================\n")
cat("Analysis completed successfully!\n")
cat("========================================\n")

sink()

cat("All analyses completed. Results saved to:\n")
cat(sprintf("  Combined analysis: %s\n", out_combined))
cat(sprintf("  Sex-separated analysis: %s\n", out_sex_separated))
cat(sprintf("  Log file: %s\n", log_file))

