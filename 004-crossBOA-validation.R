fig_width <- 4
fig_height <- 5
font_title <- 12
font_axis <- 10
font_annotation <- 8

library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)

# ===============================
# 0. 路径设置
# ===============================
output_dir <- "/home/ww/Project/Protein-Composition/result/000-横断面一致且显著的衰老"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

log_file <- file.path(output_dir, "analysis_log_p0.01.txt")
sink(log_file, split = TRUE)

cat("Analysis start - Significance threshold: 0.01\n")
cat("Output dir:", output_dir, "\n\n")

# ===============================
# 1. 读取数据
# ===============================
input_file <- "/home/ww/Project/Protein-Composition/result/1227-result/usa-iceland/merged_protein_tables.xlsx"
cat("Reading merged table:\n", input_file, "\n")

df <- read_excel(input_file)

cat("Dimension:", dim(df), "\n")
cat("Columns:\n")
print(colnames(df))
cat("\n")

# 确保 p-value 是数值
df <- df %>%
  mutate(
    iceland_pvalue_age = as.numeric(iceland_pvalue_age),
    usa_pvalue_age     = as.numeric(usa_pvalue_age)
  )

# ===============================
# 2. 读取注释文件并合并 gene symbol
# ===============================
anno_path <- "/home/ww/Project/AnzhenData/annoinfo_df_SOMA.csv"
cat("Reading annotation:\n", anno_path, "\n")

anno_df <- read.csv(anno_path, stringsAsFactors = FALSE) %>%
  select(
    seq_id = AptName,
    EntrezGeneSymbol,
    UniProt
  )

df <- df %>%
  inner_join(anno_df, by = "seq_id")

cat("After merge with annotation:", dim(df), "\n\n")

# ===============================
# 3. 显著性定义修改为 p < 0.01
# ===============================
df <- df %>%
  mutate(
    iceland_sig = !is.na(iceland_pvalue_age) & iceland_pvalue_age < 0.01,
    usa_sig     = !is.na(usa_pvalue_age) & usa_pvalue_age < 0.01,
    both_sig    = iceland_sig & usa_sig,
    direction_consistent = (iceland_beta_age * usa_beta_age) > 0
  )

cat("========== Significance Summary (p < 0.01) ==========\n")
cat("USA significant (p < 0.01):", sum(df$usa_sig, na.rm = TRUE), "\n")
cat("Iceland significant (p < 0.01):", sum(df$iceland_sig, na.rm = TRUE), "\n")
cat("Both significant (p < 0.01):", sum(df$both_sig, na.rm = TRUE), "\n")
cat("Both significant & consistent:",
    sum(df$both_sig & df$direction_consistent, na.rm = TRUE), "\n")
cat("Both significant but inconsistent:",
    sum(df$both_sig & !df$direction_consistent, na.rm = TRUE), "\n\n")

# 保存完整表格
write.csv(df, 
          file.path(output_dir, "usa_iceland_merged_full_table_p0.01.csv"), 
          row.names = FALSE)


plot_df <- df %>%
  mutate(
    sig_group = case_when(
      both_sig & direction_consistent ~ "Both sig & consistent",
      both_sig & !direction_consistent ~ "Both sig & inconsistent",
      usa_sig & !iceland_sig ~ "USA only",
      iceland_sig & !usa_sig ~ "Iceland only",
      TRUE ~ "Not significant"
    )
  ) %>%
  mutate(sig_group = factor(sig_group, 
                            levels = c("Both sig & consistent", 
                                       "Both sig & inconsistent",
                                       "USA only", 
                                       "Iceland only", 
                                       "Not significant")))

# ========== Significance Summary (p < 0.01) ==========
# USA significant (p < 0.01): 3179 
# Iceland significant (p < 0.01): 4429 
# Both significant (p < 0.01): 2877 
# Both significant & consistent: 1826 
# Both significant but inconsistent: 1051 
# ===============================
# 5. 重合显著且方向一致的结果
# ===============================
consistent_df <- df %>%
  filter(both_sig & direction_consistent) %>%
  # 保持原来的rank计算方式
  mutate(
    abs_ice = abs(iceland_beta_age),
    abs_usa = abs(usa_beta_age)
  ) %>%
  mutate(
    rank_ice = rank(-abs_ice, ties.method = "min"),
    rank_usa = rank(-abs_usa, ties.method = "min"),
    rank_sum = rank_ice + rank_usa
  ) %>%
  arrange(rank_sum) %>%
  mutate(rank_aging = row_number())

# 保存结果
write.csv(
  consistent_df,
  file.path(output_dir, "Both_significant_consistent_p0.01.csv"),
  row.names = FALSE
)

cat("Consistent significant proteins (p < 0.01):", nrow(consistent_df), "\n")
cat("Unique gene symbols:", length(unique(consistent_df$EntrezGeneSymbol)), "\n\n")

# ===============================
# 6. 绘制显著且一致的top50图
# ===============================

# 获取top50蛋白
top50 <- consistent_df %>% filter(rank_aging <= 10)
cat("[Top 1-50 Consistent Proteins]\n")
cat("Protein count:", nrow(top50), "\n")
cat("Unique gene symbols:", length(unique(top50$EntrezGeneSymbol)), "\n\n")

# 创建top50标注图
p_top50 <- ggplot(plot_df, 
                  aes(x = iceland_beta_age, y = usa_beta_age, 
                      color = sig_group)) +
  # 背景参考线
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", 
             linewidth = 0.3, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", 
             linewidth = 0.3, alpha = 0.7) +
  
  # 绘制所有蛋白点（包括不显著的）
  geom_point(alpha = 0.4, size = 1, shape = 16) +
  
  # 标注top50蛋白
  # geom_text_repel(
  #   data = top50,
  #   aes(x = iceland_beta_age, y = usa_beta_age, label = EntrezGeneSymbol),
  #   size = font_annotation / 3,
  #   color = "black",
  #   max.overlaps = Inf,
  #   min.segment.length = 0.1,
  #   segment.color = "gray50",
  #   segment.size = 0.3,
  #   box.padding = 0.3,
  #   point.padding = 0.3,
  #   force = 5,
  #   max.time = 2,
  #   max.iter = 10000
  # ) +
  
  # 高亮显示top50点
  # geom_point(data = top50,
  #            aes(x = iceland_beta_age, y = usa_beta_age),
  #            color = "black",
  #            size = 2,
  #            shape = 1,
  #            stroke = 0.8) +
  
  # 颜色方案 - 与主图保持一致
  scale_color_manual(
    values = c(
      "Both sig & consistent" = "#E41A1C",     # 红色
      "Both sig & inconsistent" = "#FF7F00",   # 橙色
      "USA only" = "#4DAF4A",                  # 绿色
      "Iceland only" = "#377EB8",              # 蓝色
      "Not significant" = "#B3B3B3"            # 灰色
    ),
    name = "Significance",
    labels = c(
      "Both sig. \n& consistent",
      "Both sig. \n& inconsistent", 
      "USA sig. only",
      "Iceland sig. only",
      "Not sig."
    )
  ) +
  
  # 坐标轴标签
  labs(
    x = expression("Iceland β"[age]),
    y = expression("USA β"[age])
  ) +
  
  # 主题设置
  theme_classic() +
  theme(
    axis.title = element_text(size = font_axis, face = "bold"),
    axis.text = element_text(size = font_axis - 1, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    legend.title = element_text(size = font_axis+1),
    legend.text = element_text(size = font_annotation+1),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.key = element_rect(fill = "white"),
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0.2, "cm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  ) +
  
  # 坐标轴比例
  coord_fixed(ratio = 3) +
  xlim(min(plot_df$iceland_beta_age, na.rm = TRUE) * 1.1,
       max(plot_df$iceland_beta_age, na.rm = TRUE) * 1.1) +
  ylim(min(plot_df$usa_beta_age, na.rm = TRUE) * 1.1,
       max(plot_df$usa_beta_age, na.rm = TRUE) * 1.1) +
  
  # 图例垂直排列，文字在点下面
  guides(color = guide_legend(
    ncol = 1,
    byrow = TRUE,
    keyheight = unit(0.9, "cm"),
    override.aes = list(size = 3, alpha = 1),
    title.position = "top",
    title.hjust = 0.5
  ))

# 保存top50图
# ggsave(
#   file.path(output_dir, "Top50_consistent_labeled_p0.01.pdf"),
#   plot = p_top50,
#   width = fig_width,
#   height = fig_height,
#   dpi = 600
# )

ggsave(
  file.path(output_dir, "Top50_consistent_labeled_p0.01.png"),
  plot = p_top50,
  width = fig_width+0.5,
  height = fig_height+0.5,
  dpi = 300
)

cat("Top50 consistent proteins labeled plot saved\n\n")

# ===============================
# 7. 绘制51-100个一致蛋白图
# ===============================
top51_100 <- consistent_df %>% 
  filter(rank_aging > 50 & rank_aging <= 100)

if (nrow(top51_100) > 0) {
  cat("[Top 51-100 Consistent Proteins]\n")
  cat("Protein count:", nrow(top51_100), "\n")
  cat("Unique gene symbols:", length(unique(top51_100$EntrezGeneSymbol)), "\n\n")
  
  # 创建top51-100标注图
  p_top51_100 <- ggplot(plot_df, 
                        aes(x = iceland_beta_age, y = usa_beta_age, 
                            color = sig_group)) +
    # 背景参考线
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", 
               linewidth = 0.3, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", 
               linewidth = 0.3, alpha = 0.7) +
    
    # 绘制所有蛋白点
    geom_point(alpha = 0.8, size = 1.5, shape = 16) +
    
    # 标注top51-100蛋白
    geom_text_repel(
      data = top51_100,
      aes(x = iceland_beta_age, y = usa_beta_age, label = EntrezGeneSymbol),
      size = font_annotation / 3,
      color = "black",
      max.overlaps = 20,
      min.segment.length = 0.1,
      segment.color = "gray50",
      segment.size = 0.3,
      box.padding = 0.3,
      point.padding = 0.3,
      force = 1
    ) +
    
    # 高亮显示top51-100点
    geom_point(data = top51_100,
               aes(x = iceland_beta_age, y = usa_beta_age),
               color = "black",
               size = 2,
               shape = 1,
               stroke = 0.8) +
    
    # 颜色方案
    scale_color_manual(
      values = c(
        "Both sig & consistent" = "#E41A1C",
        "Both sig & inconsistent" = "#FF7F00",
        "USA only" = "#4DAF4A",
        "Iceland only" = "#377EB8",
        "Not significant" = "#B3B3B3"
      ),
      name = "Significance",
      labels = c(
        "Both sig & consistent",
        "Both sig & inconsistent", 
        "USA only",
        "Iceland only",
        "Not significant"
      )
    ) +
    
    # 坐标轴标签
    labs(
      x = expression("Iceland β"[age]),
      y = expression("USA β"[age])
    ) +
    
    # 主题设置
    theme_classic() +
    theme(
      axis.title = element_text(size = font_axis, face = "bold"),
      axis.text = element_text(size = font_axis - 1, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      legend.title = element_text(size = font_axis - 1, face = "bold"),
      legend.text = element_text(size = font_annotation),
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      legend.key = element_rect(fill = "white"),
      legend.box.just = "left",
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.y = unit(0.2, "cm"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    ) +
    
    # 坐标轴比例
    coord_fixed(ratio = 1) +
    xlim(min(plot_df$iceland_beta_age, na.rm = TRUE) * 1.1,
         max(plot_df$iceland_beta_age, na.rm = TRUE) * 1.1) +
    ylim(min(plot_df$usa_beta_age, na.rm = TRUE) * 1.1,
         max(plot_df$usa_beta_age, na.rm = TRUE) * 1.1) +
    
    # 图例垂直排列
    guides(color = guide_legend(
      ncol = 1,
      byrow = TRUE,
      keyheight = unit(0.6, "cm"),
      override.aes = list(size = 3, alpha = 1),
      title.position = "top",
      title.hjust = 0.5
    ))
  
  # 保存top51-100图
  # ggsave(
  #   file.path(output_dir, "Top51_100_consistent_labeled_p0.01.pdf"),
  #   plot = p_top51_100,
  #   width = fig_width,
  #   height = fig_height,
  #   dpi = 600
  # )
  
  ggsave(
    file.path(output_dir, "Top51_100_consistent_labeled_p0.01.png"),
    plot = p_top51_100,
    width = fig_width,
    height = fig_height,
    dpi = 300
  )
  
  cat("Top51-100 consistent proteins labeled plot saved\n\n")
}

# ===============================
# 8. 绘制不一致的显著蛋白图
# ===============================
inconsistent_df <- df %>%
  filter(both_sig & !direction_consistent) %>%
  # 为不一致的蛋白添加排序
  mutate(
    discordance_score = abs(iceland_beta_age - usa_beta_age),
    rank_discordance = rank(-discordance_score, ties.method = "min")
  ) %>%
  arrange(rank_discordance) %>%
  mutate(rank_inconsistent = row_number())

if (nrow(inconsistent_df) > 0) {
  cat("[Inconsistent Significant Proteins]\n")
  cat("Protein count:", nrow(inconsistent_df), "\n")
  cat("Unique gene symbols:", length(unique(inconsistent_df$EntrezGeneSymbol)), "\n\n")
  
  # 保存不一致蛋白列表
  write.csv(inconsistent_df,
            file.path(output_dir, "Inconsistent_significant_p0.01.csv"),
            row.names = FALSE)
  
  # 选择前20个最不一致的进行标注
  top_inconsistent <- inconsistent_df %>% slice_min(rank_inconsistent, n = 20)
  
  # 创建不一致蛋白图
  p_inconsistent <- ggplot(plot_df, 
                          aes(x = iceland_beta_age, y = usa_beta_age, 
                              color = sig_group)) +
    # 背景参考线
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", 
               linewidth = 0.3, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", 
               linewidth = 0.3, alpha = 0.7) +
    
    # 绘制所有蛋白点
    geom_point(alpha = 0.8, size = 1.5, shape = 16) +
    
    # 标注不一致蛋白
    geom_text_repel(
      data = top_inconsistent,
      aes(x = iceland_beta_age, y = usa_beta_age, label = EntrezGeneSymbol),
      size = font_annotation / 3,
      color = "black",
      max.overlaps = 20,
      min.segment.length = 0.1,
      segment.color = "gray50",
      segment.size = 0.3,
      box.padding = 0.3,
      point.padding = 0.3,
      force = 1
    ) +
    
    # 高亮显示不一致蛋白点
    geom_point(data = top_inconsistent,
               aes(x = iceland_beta_age, y = usa_beta_age),
               color = "black",
               size = 2,
               shape = 1,
               stroke = 0.8) +
    
    # 颜色方案
    scale_color_manual(
      values = c(
        "Both sig & consistent" = "#E41A1C",
        "Both sig & inconsistent" = "#FF7F00",
        "USA only" = "#4DAF4A",
        "Iceland only" = "#377EB8",
        "Not significant" = "#B3B3B3"
      ),
      name = "Significance",
      labels = c(
        "Both sig & consistent",
        "Both sig & inconsistent", 
        "USA only",
        "Iceland only",
        "Not significant"
      )
    ) +
    
    # 坐标轴标签
    labs(
      x = expression("Iceland β"[age]),
      y = expression("USA β"[age])
    ) +
    
    # 主题设置
    theme_classic() +
    theme(
      axis.title = element_text(size = font_axis, face = "bold"),
      axis.text = element_text(size = font_axis - 1, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      legend.title = element_text(size = font_axis - 1, face = "bold"),
      legend.text = element_text(size = font_annotation),
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      legend.key = element_rect(fill = "white"),
      legend.box.just = "left",
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.y = unit(0.2, "cm"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    ) +
    
    # 坐标轴比例
    coord_fixed(ratio = 1) +
    xlim(min(plot_df$iceland_beta_age, na.rm = TRUE) * 1.1,
         max(plot_df$iceland_beta_age, na.rm = TRUE) * 1.1) +
    ylim(min(plot_df$usa_beta_age, na.rm = TRUE) * 1.1,
         max(plot_df$usa_beta_age, na.rm = TRUE) * 1.1) +
    
    # 图例垂直排列
    guides(color = guide_legend(
      ncol = 1,
      byrow = TRUE,
      keyheight = unit(0.6, "cm"),
      override.aes = list(size = 3, alpha = 1),
      title.position = "top",
      title.hjust = 0.5
    ))
  
  # 保存不一致蛋白图
  # ggsave(
  #   file.path(output_dir, "Inconsistent_significant_labeled_p0.01.pdf"),
  #   plot = p_inconsistent,
  #   width = fig_width,
  #   height = fig_height,
  #   dpi = 600
  # )
  
  ggsave(
    file.path(output_dir, "Inconsistent_significant_labeled_p0.01.png"),
    plot = p_inconsistent,
    width = fig_width,
    height = fig_height,
    dpi = 300
  )
  
  cat("Inconsistent significant proteins labeled plot saved\n\n")
}

# ===============================
# 9. 创建汇总表格和统计报告
# ===============================

# 创建汇总统计
summary_stats <- data.frame(
  Category = c("Total proteins",
               "USA significant (p<0.01)",
               "Iceland significant (p<0.01)",
               "Both significant (p<0.01)",
               "Both significant & consistent",
               "Both significant & inconsistent"),
  Count = c(nrow(df),
            sum(df$usa_sig, na.rm = TRUE),
            sum(df$iceland_sig, na.rm = TRUE),
            sum(df$both_sig, na.rm = TRUE),
            sum(df$both_sig & df$direction_consistent, na.rm = TRUE),
            sum(df$both_sig & !df$direction_consistent, na.rm = TRUE)),
  Percentage = c(100,
                 round(sum(df$usa_sig, na.rm = TRUE)/nrow(df)*100, 2),
                 round(sum(df$iceland_sig, na.rm = TRUE)/nrow(df)*100, 2),
                 round(sum(df$both_sig, na.rm = TRUE)/nrow(df)*100, 2),
                 round(sum(df$both_sig & df$direction_consistent, na.rm = TRUE)/nrow(df)*100, 2),
                 round(sum(df$both_sig & !df$direction_consistent, na.rm = TRUE)/nrow(df)*100, 2))
)

write.csv(summary_stats,
          file.path(output_dir, "Summary_statistics_p0.01.csv"),
          row.names = FALSE)

# 打印最终汇总
cat("========== FINAL SUMMARY (p < 0.01) ==========\n")
print(summary_stats)
cat("\n")

cat("Analysis completed successfully!\n")
cat("All results saved to:", output_dir, "\n")

sink()