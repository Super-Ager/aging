library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggh4x)  # For label alignment in circular plots
library(RColorBrewer)  # For color scheme optimization

# Create output directory
output_dir <- '/home/ww/Project/Protein-Composition/result/001-横断面top100gene的检查'
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# =============================================================================
# Step 1: Read both_significant_consistent_direction table, get top 100 seqid list
# =============================================================================
cat("Step 1: Reading both_significant_consistent_direction table...\n")
consistent_file <- '/home/ww/Project/Protein-Composition/result/1227-result/usa-iceland/both_significant_consistent_direction.csv'
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
mixed_file <- '/home/ww/Project/Protein-Composition/result/1227-result/方差分解和生存分析合并.csv'
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
    label_symbol = gene_symbol  # 因为已经用top100_symbols过滤，gene_symbol就是label_symbol
  ) %>%
  arrange(match(gene_symbol, top100_symbols))

cat("混合数据中匹配到top100的数量:", nrow(mixed_top100), "\n")

# =============================================================================
# 步骤3: 读取男性数据（年龄、个体、残差）
# =============================================================================
cat("\n步骤3: 读取男性数据...\n")
male_file <- '/home/ww/Project/Protein-Composition/result/1227-result-sex-separated/combined_protein_analysis_results_male.csv'
male_df <- read.csv(male_file, stringsAsFactors = FALSE)

cat("男性数据表格维度:", dim(male_df), "\n")
cat("男性数据表格列名:", colnames(male_df), "\n")

# 筛选top100的seqid（需要确认男性数据中seqid的列名）
# 假设列名可能是seqid或seq_id
seqid_col_male <- ifelse("seqid" %in% colnames(male_df), "seqid", 
                        ifelse("seq_id" %in% colnames(male_df), "seq_id", NA))

if (is.na(seqid_col_male)) {
  stop("无法在男性数据中找到seqid或seq_id列")
}

# 使用动态列名
male_df$seq_id <- male_df[[seqid_col_male]]

# 创建gene_symbol列（用于去重）
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

# 按gene_symbol去重，保留年龄贡献最大的行
cat("\n按gene_symbol去重（保留年龄贡献最大的行）...\n")
cat("去重前数量:", nrow(male_df), "\n")
male_df <- male_df %>%
  filter(!is.na(prop_age)) %>%
  group_by(gene_symbol) %>%
  slice_max(prop_age, n = 1, with_ties = FALSE) %>%
  ungroup()
cat("去重后数量:", nrow(male_df), "\n")

# 使用top100_symbols来筛选top100的基因
male_top100 <- male_df %>%
  filter(gene_symbol %in% top100_symbols) %>%
  mutate(
    label_symbol = gene_symbol  # 因为已经用top100_symbols过滤，gene_symbol就是label_symbol
  ) %>%
  arrange(match(gene_symbol, top100_symbols))

cat("男性数据中匹配到top100的数量:", nrow(male_top100), "\n")

# =============================================================================
# 步骤4: 读取女性数据（年龄、个体、残差）
# =============================================================================
cat("\n步骤4: 读取女性数据...\n")
female_file <- '/home/ww/Project/Protein-Composition/result/1227-result-sex-separated/combined_protein_analysis_results_female.csv'
female_df <- read.csv(female_file, stringsAsFactors = FALSE)

cat("女性数据表格维度:", dim(female_df), "\n")
cat("女性数据表格列名:", colnames(female_df), "\n")

# 筛选top100的seqid
seqid_col_female <- ifelse("seqid" %in% colnames(female_df), "seqid", 
                          ifelse("seq_id" %in% colnames(female_df), "seq_id", NA))

if (is.na(seqid_col_female)) {
  stop("无法在女性数据中找到seqid或seq_id列")
}

# 使用动态列名
female_df$seq_id <- female_df[[seqid_col_female]]

# 创建gene_symbol列（用于去重）
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

# 按gene_symbol去重，保留年龄贡献最大的行
cat("\n按gene_symbol去重（保留年龄贡献最大的行）...\n")
cat("去重前数量:", nrow(female_df), "\n")
female_df <- female_df %>%
  filter(!is.na(prop_age)) %>%
  group_by(gene_symbol) %>%
  slice_max(prop_age, n = 1, with_ties = FALSE) %>%
  ungroup()
cat("去重后数量:", nrow(female_df), "\n")

# 使用top100_symbols来筛选top100的基因
female_top100 <- female_df %>%
  filter(gene_symbol %in% top100_symbols) %>%
  mutate(
    label_symbol = gene_symbol  # 因为已经用top100_symbols过滤，gene_symbol就是label_symbol
  ) %>%
  arrange(match(gene_symbol, top100_symbols))

cat("女性数据中匹配到top100的数量:", nrow(female_top100), "\n")

# =============================================================================
# 步骤5: 准备绘图数据并绘制三个堆叠柱状图
# =============================================================================

# 函数：准备堆叠柱状图数据
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

# 图1: 混合数据（年龄、性别、个体、残差）
cat("\n步骤5: 准备混合数据绘图...\n")
mixed_cols <- c("prop_age", "prop_sex", "prop_idno", "prop_residual")
mixed_col_names <- c("Age", "Sex", "Individual", "Residual")

# 检查列是否存在
missing_mixed_cols <- setdiff(mixed_cols, colnames(mixed_top100))
if (length(missing_mixed_cols) > 0) {
  cat("警告：混合数据中缺少列:", paste(missing_mixed_cols, collapse = ", "), "\n")
  cat("可用列:", paste(colnames(mixed_top100), collapse = ", "), "\n")
  # 尝试查找相似的列名
  for (col in missing_mixed_cols) {
    similar_cols <- grep(col, colnames(mixed_top100), ignore.case = TRUE, value = TRUE)
    if (length(similar_cols) > 0) {
      cat("  可能的匹配:", paste(similar_cols, collapse = ", "), "\n")
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

# 图2: 男性数据（年龄、个体、残差）
cat("\n准备男性数据绘图...\n")
male_cols <- c("prop_age", "prop_idno", "prop_residual")
male_col_names <- c("Age", "Individual", "Residual")

missing_male_cols <- setdiff(male_cols, colnames(male_top100))
if (length(missing_male_cols) > 0) {
  cat("警告：男性数据中缺少列:", paste(missing_male_cols, collapse = ", "), "\n")
  cat("可用列:", paste(colnames(male_top100), collapse = ", "), "\n")
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

# 图3: 女性数据（年龄、个体、残差）
cat("\n准备女性数据绘图...\n")
female_cols <- c("prop_age", "prop_idno", "prop_residual")
female_col_names <- c("Age", "Individual", "Residual")

missing_female_cols <- setdiff(female_cols, colnames(female_top100))
if (length(missing_female_cols) > 0) {
  cat("警告：女性数据中缺少列:", paste(missing_female_cols, collapse = ", "), "\n")
  cat("可用列:", paste(colnames(female_top100), collapse = ", "), "\n")
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

# 绘制函数（环形堆叠柱状图）
draw_stacked_barplot <- function(data_long, title_text, fill_colors, output_filename, source_data = NULL) {
  n_genes <- length(unique(data_long$label_symbol))
  
  # 按照年龄贡献（prop_age）排序
  # 从data_long中提取每个基因的prop_age值
  gene_prop_age <- data_long %>%
    filter(variance_component == "Age") %>%
    select(label_symbol, proportion) %>%
    rename(prop_age = proportion) %>%
    group_by(label_symbol) %>%
    summarise(prop_age = first(prop_age), .groups = "drop") %>%
    arrange(desc(prop_age))  # 按年龄贡献从大到小排序
  
  # 使用排序后的基因顺序
  gene_levels <- gene_prop_age$label_symbol
  
  # 为每个基因创建索引（用于极坐标）
  data_long <- data_long %>%
    mutate(
      label_symbol = factor(label_symbol, levels = gene_levels),
      gene_index = as.numeric(label_symbol)
    )
  
  # 计算每个基因的总比例（用于确定y轴范围）
  max_prop <- data_long %>%
    group_by(gene_index) %>%
    summarise(total_prop = sum(proportion, na.rm = TRUE)) %>%
    pull(total_prop) %>%
    max(na.rm = TRUE)
  
  # 准备标签数据（使用ggh4x的geom_text_aimed实现沿半径方向对齐）
  label_data <- data.frame(
    gene_index = 1:n_genes,
    label = gene_levels,
    y_pos = max_prop * 1.15  # 标签位置在柱状图外侧
  )
  
  # 创建环形堆叠柱状图
  p <- ggplot(data_long, aes(x = gene_index, y = proportion, fill = variance_component)) +
    geom_bar(stat = "identity", position = position_stack(reverse = FALSE), width = 0.8) +
    # 使用geom_text_aimed添加沿半径方向对齐的标签
    # xend设置为基因索引（保持x不变），yend=0表示标签指向圆心（y=0），实现沿半径方向对齐
    geom_text_aimed(
      data = label_data,
      aes(x = gene_index, y = y_pos, label = label, xend = gene_index, yend = 0),
      inherit.aes = FALSE,
      size = 7,  # 字体大小（单位是mm，约等于16-18pt）
      color = "black",
      hjust = 0.5,
      vjust = 0.5
    ) +
    scale_fill_manual(
      values = fill_colors,
      name = "Variance\nComponent"
    ) +
    scale_x_continuous(
      breaks = NULL,  # 不显示默认标签
      expand = expansion(add = c(0.5, 0.5))  # 在两端添加间隔，避免第一个和最后一个基因的柱子相连
    ) +
    scale_y_continuous(
      limits = c(0, max_prop * 1.25),  # 增加空间以容纳标签
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
      axis.text.x = element_blank(),  # 隐藏默认标签
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
  # 正方形图片
  ggsave(output_file, plot = p, width = 16, height = 16, dpi = 300, bg = "white", limitsize = FALSE)
  cat(sprintf("已保存: %s (宽度=16, 高度=16, 正方形)\n", output_file))
  
  return(p)
}

# 定义颜色方案（参考图片：柔和、现代的配色）
# 使用柔和的浅色调，类似图片中的浅黄色、浅青色、淡紫色
colors_mixed <- c(
  "Age" = "#FFE5B4",      # 浅黄色/奶油色 - 年龄贡献（参考图片最左侧颜色）
  "Sex" = "#B8E6D1",      # 浅青色/薄荷绿 - 性别贡献（参考图片中间颜色）
  "Individual" = "#D4C5E8", # 淡紫色/薰衣草色 - 个体贡献（参考图片最右侧颜色）
  "Residual" = "#E8E8E8"   # 浅灰色 - 残差（中性色）
)

colors_sex_separated <- c(
  "Age" = "#FFE5B4",      # 浅黄色/奶油色 - 年龄贡献
  "Individual" = "#D4C5E8", # 淡紫色/薰衣草色 - 个体贡献
  "Residual" = "#E8E8E8"   # 浅灰色 - 残差
)

# 绘制三个图
cat("\n开始绘制堆叠柱状图...\n")

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

# 打印图片
print(p1)
print(p2)
print(p3)

# =============================================================================
# 步骤6: 统计年龄贡献阈值下的基因数量和比例
# =============================================================================
cat("\n步骤6: 统计年龄贡献阈值下的基因数量和比例...\n")

# 函数：计算统计信息
calculate_age_contribution_stats <- function(df, top100_symbols, dataset_name) {
  # 全部基因统计
  all_genes <- df %>%
    filter(!is.na(prop_age))
  
  n_all <- nrow(all_genes)
  n_all_gt02 <- sum(all_genes$prop_age > 0.2, na.rm = TRUE)
  n_all_gt01 <- sum(all_genes$prop_age > 0.1, na.rm = TRUE)
  n_all_gt005 <- sum(all_genes$prop_age > 0.05, na.rm = TRUE)
  
  prop_all_gt02 <- ifelse(n_all > 0, n_all_gt02 / n_all, 0)
  prop_all_gt01 <- ifelse(n_all > 0, n_all_gt01 / n_all, 0)
  prop_all_gt005 <- ifelse(n_all > 0, n_all_gt005 / n_all, 0)
  
  # Top100基因统计（使用gene_symbol或label_symbol来匹配）
  if ("gene_symbol" %in% colnames(df)) {
    top100_genes <- df %>%
      filter(!is.na(prop_age) & gene_symbol %in% top100_symbols)
  } else if ("label_symbol" %in% colnames(df)) {
    top100_genes <- df %>%
      filter(!is.na(prop_age) & label_symbol %in% top100_symbols)
  } else {
    # 如果都没有，返回空数据框
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
  
  # 返回结果（按照用户要求的行名格式）
  # -num行填数量，-prop和-top100行填比例
  result <- data.frame(
    row.names = c(
      paste0(dataset_name, "-all-num"),
      paste0(dataset_name, "-all-prop"),
      paste0(dataset_name, "-top100")
    ),
    "大于0.2比例" = c(n_all_gt02, prop_all_gt02, prop_top100_gt02),
    "大于0.1比例" = c(n_all_gt01, prop_all_gt01, prop_top100_gt01),
    "大于0.05比例" = c(n_all_gt005, prop_all_gt005, prop_top100_gt005),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("\n%s统计:\n", dataset_name))
  cat(sprintf("  全部基因: 总数=%d, >0.2=%d (%.4f), >0.1=%d (%.4f), >0.05=%d (%.4f)\n",
              n_all, n_all_gt02, prop_all_gt02, n_all_gt01, prop_all_gt01, 
              n_all_gt005, prop_all_gt005))
  cat(sprintf("  Top100基因: 总数=%d, >0.2=%d (%.4f), >0.1=%d (%.4f), >0.05=%d (%.4f)\n",
              n_top100, n_top100_gt02, prop_top100_gt02, n_top100_gt01, prop_top100_gt01,
              n_top100_gt005, prop_top100_gt005))
  
  return(result)
}

# 计算各数据集的统计信息
cat("\n计算混合数据统计...\n")
stats_mixed <- calculate_age_contribution_stats(mixed_df, top100_symbols, "mixed")

cat("\n计算男性数据统计...\n")
stats_male <- calculate_age_contribution_stats(male_df, top100_symbols, "male")

cat("\n计算女性数据统计...\n")
stats_female <- calculate_age_contribution_stats(female_df, top100_symbols, "female")

# 合并所有统计结果
stats_combined <- rbind(stats_mixed, stats_male, stats_female)

# 保存表格
output_stats_file <- file.path(output_dir, "age_contribution_threshold_stats.csv")
write.csv(stats_combined, output_stats_file, row.names = TRUE)
cat(sprintf("\n统计表格已保存到: %s\n", output_stats_file))

# 打印表格
cat("\n统计表格:\n")
print(stats_combined)

cat("\n分析完成！\n")
cat("输出目录:", output_dir, "\n")

