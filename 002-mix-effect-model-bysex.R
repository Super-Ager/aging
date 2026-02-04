library(parallel)
library(dplyr)
library(ggplot2)
library(lme4)
library(MuMIn)
library(lmerTest)

options(warn = -1)

# -------------------------------------------------------------------------
# 路径
# -------------------------------------------------------------------------
soma_path <- '/home/ww/Project/Longi_OA/results/tony_model/soma_df_Age_Sex_F_Dataset.csv'
anno_path <- '/home/ww/Project/AnzhenData/annoinfo_df_SOMA.csv'

out_root <- '/home/ww/Project/Protein-Composition/result/1227-result-sex-separated'
plot_dir_male <- file.path(out_root, 'png-protein-visit-batch-male')
plot_dir_female <- file.path(out_root, 'png-protein-visit-batch-female')

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir_male, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir_female, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# 数据读取 + visit 重命名
# -------------------------------------------------------------------------
soma_df <- read.csv(soma_path, stringsAsFactors = FALSE) %>%
  filter(Dataset %in% c(1,2,3,4)) %>%
  rename(visit = Dataset)

anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)

anno_map_symbol  <- setNames(anno_df$EntrezGeneSymbol, anno_df$AptName)
anno_map_uniprot <- setNames(anno_df$UniProt, anno_df$AptName)

seqid_cols <- grep("^seq\\.", names(soma_df), value = TRUE)

# -------------------------------------------------------------------------
# 单蛋白分析 + 作图（按性别分组）
# -------------------------------------------------------------------------
analyze_protein_by_sex <- function(seqid, sex_value, sex_label) {
  tryCatch({
    cat(sprintf("%s - %s\n", seqid, sex_label))
    if (!seqid %in% names(soma_df)) return(NULL)

  # 按性别筛选数据
  df <- data.frame(
    idno  = soma_df$ID,
    age   = soma_df$age,
    sex   = soma_df$Sex_F,
    visit = factor(soma_df$visit, levels = c(1,2,3,4)),
    expr  = soma_df[[seqid]]
  ) %>%
    filter(sex == sex_value) %>%
    na.omit()

  if (nrow(df) < 30) return(NULL)

  # log
  df$log_expr <- log(df$expr + 1)

  # ±3 SD
  mu  <- mean(df$log_expr)
  sdv <- sd(df$log_expr)
  df  <- df[df$log_expr >= mu - 3*sdv & df$log_expr <= mu + 3*sdv, ]
  if (nrow(df) < 30) return(NULL)

  # ---------------- batch correction ----------------
  lm_batch <- lm(log_expr ~ age + visit, data = df)

  visit_mat  <- model.matrix(lm_batch)[, grep("^visit", colnames(model.matrix(lm_batch))), drop = FALSE]
  visit_beta <- coef(lm_batch)[grep("^visit", names(coef(lm_batch)))]
  batch_eff  <- as.numeric(visit_mat %*% visit_beta)

  df$log_expr_corrected <- df$log_expr - batch_eff

  # ---------------- plotting ----------------
#   plot_df <- rbind(
#     df %>% transmute(age, expr = log_expr, visit, stage = "Pre-normalization"),
#     df %>% transmute(age, expr = log_expr_corrected, visit, stage = "Post-normalization")
#   )

#   y_lim <- range(plot_df$expr)

#   p <- ggplot(plot_df, aes(age, expr, color = visit)) +
#     geom_point(alpha = 0.25, size = 1) +
#     geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, formula = y ~ x) +
#     facet_wrap(~stage, nrow = 1) +
#     scale_y_continuous(limits = y_lim) +
#     labs(
#       x = "Age",
#       y = seqid,
#       title = paste0(seqid, " (", anno_map_symbol[seqid], ", ", anno_map_uniprot[seqid], ") - ", sex_label)
#     ) +
#     theme_bw(base_size = 12) +
#     theme(
#       legend.position = "bottom",
#       panel.grid = element_blank()
#     )

#   plot_dir <- ifelse(sex_value == 0, plot_dir_male, plot_dir_female)
#   ggsave(
#     filename = file.path(plot_dir, paste0(seqid, ".png")),
#     plot = p,
#     width = 9,
#     height = 4,
#     dpi = 150
#   )

  # ---------------- downstream model (去掉sex固定效应) ----------------
  df$z_expr <- as.numeric(scale(df$log_expr_corrected))
  df$age_cs <- as.numeric(scale(df$age))

  # 模型：z_expr ~ age_cs + (1 | idno) （没有sex）
  m <- lmer(z_expr ~ age_cs + (1 | idno), data = df)
  sm <- summary(m)

  coef_age <- sm$coefficients["age_cs", ]

  # 正确的方差分解方法（参考 Nakagawa-Schielzeth 方法）
  # ============================================================
  # 1. 固定效应方差：Var(Xβ) —— 只有age的贡献
  # ============================================================
  
  # 固定效应设计矩阵
  X <- getME(m, "X")
  
  # 固定效应系数
  beta <- fixef(m)
  
  # 固定效应线性预测子
  eta_fixed <- as.vector(X %*% beta)
  
  # 固定效应方差（σ²_fixed）
  var_fixed_total <- var(eta_fixed)
  
  # ============================================================
  # 2. age的方差贡献
  # ============================================================
  
  # age的贡献：只包含intercept和age
  beta_intercept <- beta["(Intercept)"]
  beta_age <- beta["age_cs"]
  
  # age的贡献：只包含intercept和age
  X_age_only <- cbind(1, df$age_cs)  # intercept, age
  eta_age_only <- as.vector(X_age_only %*% c(beta_intercept, beta_age))
  var_fixed_age_only <- var(eta_age_only)
  
  # 只有intercept的固定效应预测值（基线）
  eta_intercept_only <- rep(beta_intercept, nrow(df))
  var_fixed_intercept_only <- var(eta_intercept_only)  # 应该为0（常数）
  
  # age的增量方差贡献
  var_age <- var_fixed_age_only - var_fixed_intercept_only
  
  # 由于只有age一个固定效应，var_age应该等于var_fixed_total
  # 但为了数值稳定性，我们使用计算出的var_age
  var_age <- max(var_age, var_fixed_total)
  
  # ============================================================
  # 3. 随机效应方差
  # ============================================================
  
  vc <- VarCorr(m)
  var_idno <- as.numeric(attr(vc$idno, "stddev"))^2
  
  # ============================================================
  # 4. 残差方差
  # ============================================================
  
  var_res <- sigma(m)^2
  
  # ============================================================
  # 5. 总方差（Nakagawa–Schielzeth 定义）
  # ============================================================
  
  var_tot <- var_fixed_total + var_idno + var_res
  
  # ============================================================
  # 6. 计算R²用于输出
  # ============================================================
  
  R2_full <- r.squaredGLMM(m)
  R2m_full <- R2_full[1, "R2m"]
  R2c_full <- R2_full[1, "R2c"]
  
  # 验证：prop_age 应该等于 var_fixed_total / var_tot = R2m
  # 验证：prop_age + prop_idno 应该等于 (var_fixed_total + var_idno) / var_tot = R2c

  data.frame(
    seqid = seqid,
    gene_symbol = anno_map_symbol[seqid],
    uniprot = anno_map_uniprot[seqid],
    sex = sex_label,

    n_samples = nrow(df),
    n_individuals = length(unique(df$idno)),

    age_Estimate = coef_age["Estimate"],
    age_StdError = coef_age["Std. Error"],
    age_tvalue   = coef_age["t value"],
    age_pvalue   = coef_age["Pr(>|t|)"],

    variance_age = var_age,
    variance_idno = var_idno,
    variance_residual = var_res,
    variance_total = var_tot,

    prop_age = var_age / var_tot,
    prop_idno = var_idno / var_tot,
    prop_residual = var_res / var_tot,

    R2m = R2m_full,
    R2c = R2c_full,

    stringsAsFactors = FALSE
  )
  }, error = function(e) {
    cat(sprintf("错误: %s - %s: %s\n", seqid, sex_label, conditionMessage(e)))
    return(NULL)
  })
}

# -------------------------------------------------------------------------
# 并行运行（分别处理男性和女性）
# -------------------------------------------------------------------------
n_cores <- max(1, min(detectCores() - 2, 8))  # 限制最大核心数，避免资源问题

# 处理男性（sex = 0）
cat("\n========== 开始处理男性数据 ==========\n")
cat(sprintf("使用 %d 个核心进行并行处理\n", n_cores))
res_male <- tryCatch({
  mclapply(seqid_cols, function(seqid) {
    analyze_protein_by_sex(seqid, sex_value = 0, sex_label = "Male")
  }, mc.cores = n_cores, mc.preschedule = FALSE)  # mc.preschedule = FALSE 可以减少SIGPIPE错误
}, error = function(e) {
  cat(sprintf("并行处理出错: %s\n", conditionMessage(e)))
  cat("尝试使用单核处理...\n")
  lapply(seqid_cols, function(seqid) {
    analyze_protein_by_sex(seqid, sex_value = 0, sex_label = "Male")
  })
})
res_male <- res_male[!sapply(res_male, is.null)]
if (length(res_male) > 0) {
  df_male <- do.call(rbind, res_male)
} else {
  df_male <- data.frame()
  cat("警告: 没有成功处理的男性数据\n")
}

# 处理女性（sex = 1）
cat("\n========== 开始处理女性数据 ==========\n")
cat(sprintf("使用 %d 个核心进行并行处理\n", n_cores))
res_female <- tryCatch({
  mclapply(seqid_cols, function(seqid) {
    analyze_protein_by_sex(seqid, sex_value = 1, sex_label = "Female")
  }, mc.cores = n_cores, mc.preschedule = FALSE)
}, error = function(e) {
  cat(sprintf("并行处理出错: %s\n", conditionMessage(e)))
  cat("尝试使用单核处理...\n")
  lapply(seqid_cols, function(seqid) {
    analyze_protein_by_sex(seqid, sex_value = 1, sex_label = "Female")
  })
})
res_female <- res_female[!sapply(res_female, is.null)]
if (length(res_female) > 0) {
  df_female <- do.call(rbind, res_female)
} else {
  df_female <- data.frame()
  cat("警告: 没有成功处理的女性数据\n")
}

# -------------------------------------------------------------------------
# 保存结果
# -------------------------------------------------------------------------
if (nrow(df_male) > 0) {
  write.csv(
    df_male,
    file.path(out_root, "combined_protein_analysis_results_male.csv"),
    row.names = FALSE
  )
} else {
  cat("警告: 男性数据为空，未保存文件\n")
}

if (nrow(df_female) > 0) {
  write.csv(
    df_female,
    file.path(out_root, "combined_protein_analysis_results_female.csv"),
    row.names = FALSE
  )
} else {
  cat("警告: 女性数据为空，未保存文件\n")
}

cat("\n========== 完成 ==========\n")
cat(sprintf("男性结果: %d 个蛋白质\n", nrow(df_male)))
cat(sprintf("女性结果: %d 个蛋白质\n", nrow(df_female)))
cat(sprintf("结果保存到: %s\n", out_root))

