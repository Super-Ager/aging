library(parallel)
library(dplyr)
library(survival)
library(lme4)  # 用于混合效应模型
library(glmnet)  # 用于Lasso回归
library(rms)  # 用于ANOVA和R2计算
library(ggplot2)
library(gridExtra)
library(scales)

options(warn = -1)

# -------------------------------------------------------------------------
# 路径
# -------------------------------------------------------------------------
soma_path <- '/home/ww/Project/Longi_OA/results/tony_model/soma_df_Age_Sex_F_Dataset.csv'
anno_path <- '/home/ww/Project/AnzhenData/annoinfo_df_SOMA.csv'
outcome_path <- '/home/ww/Project/Protein-Composition/data/CMCS_outcome_final.csv'

out_root <- '/home/ww/Project/Protein-Composition/result/004-longevity-target'
plot_dir <- file.path(out_root, 'plots')
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# 数据读取
# -------------------------------------------------------------------------
cat("正在读取数据...\n")
soma_df <- read.csv(soma_path, stringsAsFactors = FALSE) %>%
  filter(Dataset %in% c(1,2,3,4)) %>%
  rename(visit = Dataset)

anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
anno_map_symbol  <- setNames(anno_df$EntrezGeneSymbol, anno_df$AptName)
anno_map_uniprot <- setNames(anno_df$UniProt, anno_df$AptName)

# 读取死亡信息
outcome_df <- read.csv(outcome_path, stringsAsFactors = FALSE)

# 只保留出现在soma_df中的ID
outcome_df <- outcome_df %>%
  filter(ID %in% soma_df$ID)

# 转换死亡日期
death_date_original <- outcome_df$death_Date
outcome_df$death_Date <- as.Date(outcome_df$death_Date)

# 检查转换失败的数量
conversion_failed <- !is.na(death_date_original) & death_date_original != "" & is.na(outcome_df$death_Date)
if (sum(conversion_failed) > 0) {
  cat(sprintf("  警告：%d 个death_Date转换失败，尝试其他格式...\n", sum(conversion_failed)))
  failed_indices <- which(conversion_failed)
  failed_dates <- death_date_original[failed_indices]
  
  for (fmt in c("%Y-%m-%d", "%Y/%m/%d", "%Y-%m-%d %H:%M:%S", "%Y/%m/%d %H:%M:%S")) {
    if (sum(conversion_failed) == 0) break
    try_dates <- as.Date(failed_dates, format = fmt)
    success <- !is.na(try_dates)
    if (any(success)) {
      outcome_df$death_Date[failed_indices[success]] <- try_dates[success]
      conversion_failed[failed_indices[success]] <- FALSE
    }
  }
}

# 最后更新日期
last_update_date <- as.Date("2024/8/21", format = "%Y/%m/%d")

seqid_cols <- grep("^seq\\.", names(soma_df), value = TRUE)
cat(sprintf("找到 %d 个seqid\n", length(seqid_cols)))

# 获取所有ID（用于最终表格）
all_ids <- unique(soma_df$ID)
all_ids <- sort(all_ids)

# -------------------------------------------------------------------------
# 辅助函数：计算Cox模型的R2（基于部分似然比统计量）
# -------------------------------------------------------------------------
calculate_cox_r2 <- function(cox_model) {
  if (is.null(cox_model)) {
    return(0)
  }
  
  if (is.null(cox_model$loglik) || length(cox_model$loglik) < 2) {
    return(0)
  }
  
  if (is.null(cox_model$nevent) || is.na(cox_model$nevent) || cox_model$nevent <= 0) {
    return(0)
  }
  
  loglik_initial <- cox_model$loglik[1]
  loglik_final <- cox_model$loglik[2]
  
  if (is.na(loglik_initial) || is.na(loglik_final)) {
    return(0)
  }
  
  logtest <- -2 * (loglik_initial - loglik_final)
  nevent <- cox_model$nevent
  
  if (!is.na(logtest) && logtest > 0 && nevent > 0) {
    r2 <- 1 - exp(-logtest / nevent)
    r2 <- max(0, min(1, r2))
    if (is.na(r2)) {
      r2 <- 0
    }
  } else {
    r2 <- 0
  }
  
  return(r2)
}

# -------------------------------------------------------------------------
# 辅助函数：计算C-index及其置信区间
# -------------------------------------------------------------------------
calculate_cindex_ci <- function(cox_model) {
  tryCatch({
    cindex_obj <- concordance(cox_model)
    cindex <- as.numeric(cindex_obj$concordance)
    
    # 计算置信区间（使用标准误）
    if (!is.null(cindex_obj$var)) {
      se <- sqrt(cindex_obj$var)
      ci_lower <- cindex - 1.96 * se
      ci_upper <- cindex + 1.96 * se
      ci_lower <- max(0, min(1, ci_lower))
      ci_upper <- max(0, min(1, ci_upper))
    } else {
      ci_lower <- NA
      ci_upper <- NA
    }
    
    return(list(cindex = cindex, ci_lower = ci_lower, ci_upper = ci_upper))
  }, error = function(e) {
    return(list(cindex = NA, ci_lower = NA, ci_upper = NA))
  })
}

# -------------------------------------------------------------------------
# 单蛋白轨迹建模和死亡分析
# -------------------------------------------------------------------------
analyze_protein_gap_survival <- function(seqid) {
  cat(sprintf("分析 %s...\n", seqid))
  
  # 准备数据
  df <- data.frame(
    idno  = soma_df$ID,
    age   = soma_df$age,
    sex   = soma_df$Sex_F,
    visit = factor(soma_df$visit, levels = c(1,2,3,4)),
    expr  = soma_df[[seqid]]
  ) %>% na.omit()
  
  if (nrow(df) < 30) return(NULL)
  
  # log转换
  df$log_expr <- log(df$expr + 1)
  
  # ±3 SD过滤
  mu  <- mean(df$log_expr)
  sdv <- sd(df$log_expr)
  df  <- df[df$log_expr >= mu - 3*sdv & df$log_expr <= mu + 3*sdv, ]
  if (nrow(df) < 30) return(NULL)
  
  # 批次效应矫正
  lm_batch <- lm(log_expr ~ age + visit, data = df)
  visit_mat  <- model.matrix(lm_batch)[, grep("^visit", colnames(model.matrix(lm_batch))), drop = FALSE]
  visit_beta <- coef(lm_batch)[grep("^visit", names(coef(lm_batch)))]
  batch_eff  <- as.numeric(visit_mat %*% visit_beta)
  df$log_expr_corrected <- df$log_expr - batch_eff
  
  # =======================================================================
  # 步骤1：筛选有visit 1、2、3的个体，得到df_V123
  # =======================================================================
  idno_V123 <- df %>%
    filter(visit %in% c(1, 2, 3)) %>%
    group_by(idno) %>%
    summarise(
      has_visit1 = any(visit == 1),
      has_visit2 = any(visit == 2),
      has_visit3 = any(visit == 3),
      .groups = 'drop'
    ) %>%
    filter(has_visit1 & has_visit2 & has_visit3) %>%
    pull(idno)
  
  if (length(idno_V123) < 20) return(NULL)
  
  # 提取这些个体的visit 1、2、3数据
  df_V123 <- df %>%
    filter(idno %in% idno_V123 & visit %in% c(1, 2, 3)) %>%
    select(idno, visit, age, sex, log_expr_corrected) %>%
    arrange(idno, visit)
  
  if (nrow(df_V123) < 20) return(NULL)
  
  # =======================================================================
  # 步骤2：横断面模型（使用visit 2的数据）
  # =======================================================================
  # 提取visit 2的数据用于拟合横断面模型
  df_V2 <- df_V123 %>%
    filter(visit == 2) %>%
    select(idno, age, sex, log_expr_corrected)
  
  if (nrow(df_V2) < 10) return(NULL)
  
  # 步骤2.1：使用visit 2计算年龄的均值和标准差（用于标准化visit 2和visit 3）
  age_mean <- mean(df_V2$age, na.rm = TRUE)
  age_sd <- sd(df_V2$age, na.rm = TRUE)
  
  # 步骤2.2：使用上述均值和标准差对visit 2的age进行标准化
  df_V2$age_cs <- (df_V2$age - age_mean) / age_sd
  
  # 步骤2.3：使用visit 2的数据（标准化后的age和sex）计算回归线
  crosssectional_model <- tryCatch({
    lm(log_expr_corrected ~ age_cs + sex, data = df_V2)
  }, error = function(e) {
    cat(sprintf("  横断面模型拟合失败: %s\n", e$message))
    return(NULL)
  })
  
  if (is.null(crosssectional_model)) return(NULL)
  
  # 步骤2.4：计算每个人在visit 2时候距离回归线的gap（残差）
  df_V2$predicted_V2 <- predict(crosssectional_model, newdata = df_V2)
  df_V2$level_gap <- df_V2$log_expr_corrected - df_V2$predicted_V2
  level_gap_map <- setNames(df_V2$level_gap, df_V2$idno)
  
  # =======================================================================
  # 步骤3：准备visit 3的数据和特征
  # =======================================================================
  # 提取visit 1、2、3的数据用于建模和预测
  df_V1 <- df_V123 %>%
    filter(visit == 1) %>%
    select(idno, visit1_level = log_expr_corrected)
  
  df_V2_data <- df_V123 %>%
    filter(visit == 2) %>%
    select(idno, visit2_level = log_expr_corrected)
  
  # 提取visit 3的数据，并使用visit 2的age均值和标准差进行标准化
  df_V3 <- df_V123 %>%
    filter(visit == 3) %>%
    select(idno, age, sex, actual_visit3 = log_expr_corrected) %>%
    mutate(age_cs = (age - age_mean) / age_sd)  # 使用visit 2的参数标准化visit 3的age
  
  if (nrow(df_V3) < 20) return(NULL)
  
  # 合并所有特征数据（确保只使用有visit 1、2、3的个体）
  df_V3 <- df_V3 %>%
    left_join(df_V1, by = "idno") %>%
    left_join(df_V2_data, by = "idno") %>%
    filter(!is.na(visit1_level) & !is.na(visit2_level))
  
  if (nrow(df_V3) < 20) return(NULL)
  
  # 步骤3.1：计算每个人在visit 3时候，在回归线上的数值（使用标准化后的age和sex）
  df_V3$predicted_V3_cross <- predict(crosssectional_model, newdata = df_V3)
  
  # 步骤3.2：加上visit 2的level gap作为visit 3的预测值
  df_V3$level_gap <- level_gap_map[df_V3$idno]
  df_V3$predicted_visit3_cross <- df_V3$predicted_V3_cross + df_V3$level_gap
  
  # 计算横断面模型的预测准确性
  cor_cross <- tryCatch({
    cor.test(df_V3$actual_visit3, df_V3$predicted_visit3_cross, use = "complete.obs")
  }, error = function(e) {
    return(list(estimate = NA, p.value = NA))
  })
  cor_cross_r <- ifelse(is.null(cor_cross$estimate), NA, cor_cross$estimate)
  cor_cross_p <- ifelse(is.null(cor_cross$p.value), NA, cor_cross$p.value)
  
  ss_res_cross <- sum((df_V3$actual_visit3 - df_V3$predicted_visit3_cross)^2, na.rm = TRUE)
  ss_tot_cross <- sum((df_V3$actual_visit3 - mean(df_V3$actual_visit3, na.rm = TRUE))^2, na.rm = TRUE)
  r2_cross <- ifelse(ss_tot_cross > 0, 1 - (ss_res_cross / ss_tot_cross), NA)
  
  # 计算横断面模型的gap
  df_V3$gap_cross <- df_V3$actual_visit3 - df_V3$predicted_visit3_cross
  
  # =======================================================================
  # 步骤4：Lasso预测模型（使用visit 1的水平、visit 2的水平和性别）
  # =======================================================================
  # 准备特征矩阵和响应变量（只使用df_V3中的个体，确保有visit 1、2、3）
  # 特征：visit1_level, visit2_level, sex
  # 输出：actual_visit3
  X <- as.matrix(df_V3[, c("visit1_level", "visit2_level", "sex")])
  y <- df_V3$actual_visit3
  n_individuals <- nrow(df_V3)
  
  # 5-fold交叉验证，按个体划分
  n_folds <- 5
  fold_size <- ceiling(n_individuals / n_folds)
  
  # 随机打乱个体顺序（但保持可重复性）
  set.seed(123)
  shuffled_indices <- sample(1:n_individuals)
  
  # 初始化预测值
  predicted_visit3_lasso <- rep(NA, n_individuals)
  
  # 对每个fold进行交叉验证
  for (fold in 1:n_folds) {
    # 确定测试集索引
    start_idx <- (fold - 1) * fold_size + 1
    end_idx <- min(fold * fold_size, n_individuals)
    test_indices <- shuffled_indices[start_idx:end_idx]
    train_indices <- shuffled_indices[-test_indices]
    
    if (length(train_indices) < 10 || length(test_indices) < 1) next
    
    # 准备训练集和测试集
    X_train <- X[train_indices, , drop = FALSE]
    y_train <- y[train_indices]
    X_test <- X[test_indices, , drop = FALSE]
    
    # 使用Lasso进行回归（alpha=1，glmnet会自动标准化）
    tryCatch({
      # 使用Lasso（alpha=1）进行交叉验证
      cv_fit <- cv.glmnet(X_train, y_train, alpha = 1, nfolds = 5)
      
      # 预测（使用lambda.1se以提高泛化性）
      predictions <- predict(cv_fit, newx = X_test, s = "lambda.1se")
      
      # 保存预测值（test_indices已经是原始数据的索引）
      if (length(test_indices) > 0 && length(predictions) > 0) {
        predicted_visit3_lasso[test_indices] <- as.numeric(predictions)
      }
      
    }, error = function(e) {
      cat(sprintf("  Fold %d 预测失败: %s\n", fold, e$message))
    })
  }
  
  # 将预测值添加到df_V3（df_V3已经确保了有visit 1、2、3的个体）
  df_V3$predicted_visit3_linear <- predicted_visit3_lasso
  
  # 只保留有完整预测值的个体（横断面和Lasso模型都有预测值）
  df_V3 <- df_V3 %>%
    filter(!is.na(predicted_visit3_linear) & !is.na(predicted_visit3_cross) & 
           !is.na(actual_visit3))
  
  if (nrow(df_V3) < 10) return(NULL)
  
  # 计算Lasso模型的预测准确性
  cor_linear <- tryCatch({
    cor.test(df_V3$actual_visit3, df_V3$predicted_visit3_linear, use = "complete.obs")
  }, error = function(e) {
    return(list(estimate = NA, p.value = NA))
  })
  cor_linear_r <- ifelse(is.null(cor_linear$estimate), NA, cor_linear$estimate)
  cor_linear_p <- ifelse(is.null(cor_linear$p.value), NA, cor_linear$p.value)
  
  ss_res_linear <- sum((df_V3$actual_visit3 - df_V3$predicted_visit3_linear)^2, na.rm = TRUE)
  ss_tot_linear <- sum((df_V3$actual_visit3 - mean(df_V3$actual_visit3, na.rm = TRUE))^2, na.rm = TRUE)
  r2_linear <- ifelse(ss_tot_linear > 0, 1 - (ss_res_linear / ss_tot_linear), NA)
  
  # 计算Lasso模型的gap
  df_V3$gap_linear <- df_V3$actual_visit3 - df_V3$predicted_visit3_linear
  
  # 使用Lasso模型的预测值作为predicted_visit3_mixed（用于后续死亡分析）
  df_V3$predicted_visit3_mixed <- df_V3$predicted_visit3_linear
  df_V3$gap <- df_V3$gap_linear
  
  # 标准化predicted_visit3和gap（用于死亡分析）
  df_V3$predicted_visit3_sd <- as.numeric(scale(df_V3$predicted_visit3_mixed))
  df_V3$gap_sd <- as.numeric(scale(df_V3$gap))
  
  # =======================================================================
  # 合并死亡信息（df_V3已经确保了有visit 1、2、3的个体）
  # =======================================================================
  df_survival <- df_V3 %>%
    left_join(outcome_df %>% select(ID, death_Date, death), by = c("idno" = "ID"))
  
  if (nrow(df_survival) < 20 || sum(df_survival$death == 1, na.rm = TRUE) < 5) return(NULL)
  
  # 准备生存数据
  visit3_baseline_date <- as.Date("2013-01-01")  # 根据实际情况调整
  
  df_survival <- df_survival %>%
    mutate(
      time = ifelse(
        death == 1 & !is.na(death_Date),
        as.numeric(death_Date - visit3_baseline_date) / 365.25,
        as.numeric(last_update_date - visit3_baseline_date) / 365.25
      ),
      event = ifelse(death == 1 & !is.na(death_Date), 1, 0)
    ) %>%
    filter(time > 0)
  
  if (nrow(df_survival) < 20) return(NULL)
  
  # 准备协变量
  df_survival$age_visit3_sd <- as.numeric(scale(df_survival$age))
  df_survival$sex_factor <- factor(df_survival$sex, levels = c(0,1), labels = c("Male","Female"))
  
  # 初始化结果
  result <- list(
    seqid = seqid,
    n_samples = nrow(df_survival),
    n_death = sum(df_survival$death == 1, na.rm = TRUE),
    # 预测准确性（横断面模型）
    cross_cor_r = cor_cross_r,
    cross_cor_p = cor_cross_p,
    cross_r2 = r2_cross,
    # 预测准确性（线性模型）
    linear_cor_r = cor_linear_r,
    linear_cor_p = cor_linear_p,
    linear_r2 = r2_linear,
    # 模型结果
    base_model_r2 = NA,
    base_model_cindex = NA,
    base_model_cindex_ci_lower = NA,
    base_model_cindex_ci_upper = NA,
    full_model_r2 = NA,
    full_model_cindex = NA,
    full_model_cindex_ci_lower = NA,
    full_model_cindex_ci_upper = NA,
    lrt_chisq = NA,
    lrt_p = NA,
    gap_hr = NA,
    gap_hr_ci_lower = NA,
    gap_hr_ci_upper = NA,
    gap_p = NA
  )
  
  # 拟合模型并计算指标
  tryCatch({
    # 确保数据没有NA值
    df_model <- df_survival %>%
      select(time, event, age_visit3_sd, sex_factor, predicted_visit3_sd, gap_sd) %>%
      filter(!is.na(time) & !is.na(event) & !is.na(age_visit3_sd) & 
             !is.na(sex_factor) & !is.na(predicted_visit3_sd))
    
    if (nrow(df_model) < 20 || sum(df_model$event) < 5) {
      stop("数据不足或事件数太少")
    }
    
    # 基准模型：age + sex + predicted_visit3
    base_model <- coxph(
      Surv(time, event) ~ age_visit3_sd + sex_factor + predicted_visit3_sd,
      data = df_model
    )
    
    # 完整模型：age + sex + predicted_visit3 + gap
    full_model <- coxph(
      Surv(time, event) ~ age_visit3_sd + sex_factor + predicted_visit3_sd + gap_sd,
      data = df_model
    )
    
    # 计算R²
    result$base_model_r2 <- calculate_cox_r2(base_model)
    result$full_model_r2 <- calculate_cox_r2(full_model)
    
    # 计算C-index及其置信区间
    base_cindex <- calculate_cindex_ci(base_model)
    result$base_model_cindex <- base_cindex$cindex
    result$base_model_cindex_ci_lower <- base_cindex$ci_lower
    result$base_model_cindex_ci_upper <- base_cindex$ci_upper
    
    full_cindex <- calculate_cindex_ci(full_model)
    result$full_model_cindex <- full_cindex$cindex
    result$full_model_cindex_ci_lower <- full_cindex$ci_lower
    result$full_model_cindex_ci_upper <- full_cindex$ci_upper
    
    # Likelihood Ratio Test
    lrt <- anova(base_model, full_model)
    result$lrt_chisq <- lrt$Chisq[2]
    result$lrt_p <- lrt$`Pr(>|Chi|)`[2]
    
    # 提取gap的HR及其置信区间
    gap_coef <- summary(full_model)$coefficients["gap_sd", ]
    result$gap_hr <- exp(gap_coef["coef"])
    result$gap_hr_ci_lower <- exp(gap_coef["coef"] - 1.96 * gap_coef["se(coef)"])
    result$gap_hr_ci_upper <- exp(gap_coef["coef"] + 1.96 * gap_coef["se(coef)"])
    result$gap_p <- gap_coef["Pr(>|z|)"]
    
  }, error = function(e) {
    cat(sprintf("  模型拟合错误: %s\n", e$message))
  })
  
  # 准备预测值表格数据
  prediction_df <- data.frame(
    idno = all_ids,
    stringsAsFactors = FALSE
  )
  
  # 为当前蛋白添加三列
  col_actual <- paste0(seqid, "_actual")
  col_mixed <- paste0(seqid, "_mixed")
  col_cross <- paste0(seqid, "_cross")
  
  prediction_df[[col_actual]] <- NA
  prediction_df[[col_mixed]] <- NA
  prediction_df[[col_cross]] <- NA
  
  # 填充有数据的个体（df_V3已经确保了有visit 1、2、3和完整预测值）
  for (i in seq_len(nrow(df_V3))) {
    idx <- which(prediction_df$idno == df_V3$idno[i])
    if (length(idx) == 1) {
      prediction_df[[col_actual]][idx] <- df_V3$actual_visit3[i]
      prediction_df[[col_mixed]][idx] <- df_V3$predicted_visit3_mixed[i]
      prediction_df[[col_cross]][idx] <- df_V3$predicted_visit3_cross[i]
    }
  }
  
  # 返回结果和预测值表格
  return(list(
    result = data.frame(result, stringsAsFactors = FALSE),
    prediction = prediction_df
  ))
}

# -------------------------------------------------------------------------
# 主程序：并行分析
# -------------------------------------------------------------------------
cat("开始轨迹建模和死亡分析...\n")
cat(sprintf("总共需要分析 %d 个seqid\n", length(seqid_cols)))

# 测试代码：只分析指定的seqid
# seqid_cols <- c("seq.4374.45")
# seqid_cols <- seqid_cols[1:10]

# 设置并行计算
n_cores <- 14
cat(sprintf("使用 %d 个核心进行并行计算...\n", n_cores))

# 使用mclapply进行并行计算
cat("开始并行分析...\n")
result_lists <- mclapply(seq_along(seqid_cols), function(i) {
  seqid <- seqid_cols[i]
  
  result_list <- tryCatch({
    analyze_protein_gap_survival(seqid)
  }, error = function(e) {
    return(NULL)
  })
  
  # 每100个输出一次进度
  if (i %% 100 == 0) {
    cat(sprintf("进度: %d/%d (%.1f%%)\n", i, length(seqid_cols), 100 * i / length(seqid_cols)))
  }
  
  return(result_list)
}, mc.cores = n_cores)

# 处理结果
res <- list()
prediction_tables <- list()

for (i in seq_along(result_lists)) {
  result_list <- result_lists[[i]]
  
  if (is.null(result_list)) {
    next
  }
  
  # 检查结果是否有效
  if (is.null(result_list$result) || !is.data.frame(result_list$result)) {
    next
  }
  
  # 保存结果
  res[[i]] <- result_list$result
  
  # 保存预测值表格（第一次需要初始化，后续需要合并）
  if (!is.null(result_list$prediction) && is.data.frame(result_list$prediction)) {
    if (length(prediction_tables) == 0 || is.null(prediction_tables[[1]])) {
      prediction_tables[[1]] <- result_list$prediction
    } else {
      # 合并新列
      new_cols <- result_list$prediction[, -1, drop = FALSE]  # 排除idno列
      if (ncol(new_cols) > 0) {
        prediction_tables[[1]] <- cbind(prediction_tables[[1]], new_cols)
      }
    }
  }
}

if (length(res) == 0) {
  cat("没有有效结果\n")
  stop("没有有效结果")
}

# 移除NULL结果
res <- res[!sapply(res, is.null)]

survival_df <- do.call(rbind, res)
cat(sprintf("\n最终保存的结果: %d 行 x %d 列\n", nrow(survival_df), ncol(survival_df)))

# 保存结果
output_path <- file.path(out_root, "protein_gap_survival_results.csv")
write.csv(survival_df, output_path, row.names = FALSE)
cat(sprintf("结果已保存到: %s\n", output_path))

# 保存预测值表格
prediction_path <- file.path(out_root, "protein_visit3_predictions.csv")
if (length(prediction_tables) > 0 && !is.null(prediction_tables[[1]])) {
  prediction_final <- prediction_tables[[1]]
  
  # 将NA替换为NaN（用于Python兼容性）
  for (col in names(prediction_final)) {
    if (col != "idno") {
      prediction_final[[col]][is.na(prediction_final[[col]])] <- NaN
    }
  }
  
  write.csv(prediction_final, prediction_path, row.names = FALSE)
  cat(sprintf("预测值表格已保存到: %s\n", prediction_path))
}

# -------------------------------------------------------------------------
# 绘制linear模型预测值和真实值的散点图（针对Filtered_Proteins）
# -------------------------------------------------------------------------
cat("\n========================================\n")
cat("绘制linear模型预测准确性散点图\n")
cat("========================================\n\n")

# 读取Filtered_Proteins_VarianceDecomposition
variance_decomp_path <- '/home/ww/Project/Protein-Composition/result/001-横断面衰老marker的检测/Filtered_Proteins_VarianceDecomposition.csv'
if (file.exists(variance_decomp_path)) {
  cat("正在读取Filtered_Proteins_VarianceDecomposition...\n")
  variance_df <- read.csv(variance_decomp_path, stringsAsFactors = FALSE)
  
  # 检查seq_id列名（可能是seq_id或seqid）
  seq_id_col <- if ("seq_id" %in% names(variance_df)) "seq_id" else "seqid"
  filtered_seqids <- unique(variance_df[[seq_id_col]])
  cat(sprintf("  找到 %d 个筛选的蛋白\n", length(filtered_seqids)))
  
  # 创建散点图输出目录
  scatter_plot_dir <- file.path(out_root, "linear_prediction_scatter_plots")
  dir.create(scatter_plot_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("  散点图将保存到: %s\n", scatter_plot_dir))
  
  # 读取预测值表格
  if (file.exists(prediction_path)) {
    cat("正在读取预测值表格...\n")
    prediction_df <- read.csv(prediction_path, stringsAsFactors = FALSE)
    cat(sprintf("  预测值表格: %d 行 x %d 列\n", nrow(prediction_df), ncol(prediction_df)))
    
    # 读取生存分析结果获取pearson_r
    if (file.exists(output_path)) {
      cat("正在读取生存分析结果...\n")
      survival_results <- read.csv(output_path, stringsAsFactors = FALSE)
      
      # 为每个筛选的蛋白绘制散点图
      plot_count <- 0
      for (seqid in filtered_seqids) {
        # 检查是否有对应的列
        col_actual <- paste0(seqid, "_actual")
        col_mixed <- paste0(seqid, "_mixed")
        
        if (col_actual %in% names(prediction_df) && col_mixed %in% names(prediction_df)) {
          # 提取实际值和预测值
          plot_data <- prediction_df %>%
            select(idno, actual = all_of(col_actual), predicted = all_of(col_mixed)) %>%
            filter(!is.na(actual) & !is.na(predicted) & 
                   !is.nan(actual) & !is.nan(predicted) &
                   is.finite(actual) & is.finite(predicted))
          
          if (nrow(plot_data) >= 10) {
            # 计算pearson相关系数
            pearson_r <- tryCatch({
              cor_result <- cor.test(plot_data$actual, plot_data$predicted, use = "complete.obs")
              cor_result$estimate
            }, error = function(e) {
              NA
            })
            
            # 如果survival_results中有，使用那里的值
            if (seqid %in% survival_results$seqid) {
              pearson_r <- survival_results$linear_cor_r[survival_results$seqid == seqid][1]
            }
            
            # 获取gene_symbol
            gene_symbol <- if (seqid %in% names(anno_map_symbol)) {
              anno_map_symbol[seqid]
            } else {
              seqid
            }
            
            # 如果gene_symbol为空或NA，使用seqid
            if (is.na(gene_symbol) || gene_symbol == "" || is.null(gene_symbol)) {
              gene_symbol <- seqid
            }
            
            # 绘制散点图
            scatter_plot <- ggplot(plot_data, aes(x = actual, y = predicted)) +
              geom_point(color = "#E74C3C", alpha = 0.6, size = 1.5) +
              geom_abline(intercept = 0, slope = 1, color = "#2C3E50", linetype = "dashed", linewidth = 0.8, alpha = 0.8) +
              # 添加pearson_r标注
              annotate("text", 
                       x = Inf, y = -Inf, 
                       label = sprintf("Pearson r = %.3f", ifelse(is.na(pearson_r), 0, pearson_r)),
                       hjust = 1.1, vjust = -0.5,
                       size = 4, color = "#2C3E50", fontface = "bold") +
              labs(
                x = "Actual Visit 3 Level",
                y = "Predicted Visit 3 Level (Personalized Longitudinal Model)",
                title = sprintf("%s (%s)", gene_symbol, seqid)
              ) +
              theme_minimal(base_size = 12) +
              theme(
                axis.line = element_line(color = "black", linewidth = 0.6),
                panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
                panel.grid.minor = element_blank(),
                axis.title = element_text(size = 11, color = "#2C3E50"),
                axis.text = element_text(size = 10, color = "#2C3E50"),
                plot.title = element_text(size = 13, color = "#2C3E50", hjust = 0.5, margin = margin(b = 10)),
                panel.border = element_blank(),
                plot.margin = margin(15, 15, 15, 15),
                panel.background = element_rect(fill = "white", color = NA)
              )
            
            # 保存图片
            plot_filename <- sprintf("%s_%s_scatter.png", 
                                    gsub("[^A-Za-z0-9]", "_", gene_symbol),
                                    gsub("\\.", "_", seqid))
            plot_filepath <- file.path(scatter_plot_dir, plot_filename)
            
            ggsave(
              plot_filepath,
              scatter_plot,
              width = 4,
              height = 3,
              dpi = 300,
              bg = "white"
            )
            
            plot_count <- plot_count + 1
            if (plot_count %% 50 == 0) {
              cat(sprintf("  已绘制 %d 个散点图...\n", plot_count))
            }
          }
        }
      }
      
      cat(sprintf("\n✅ 共绘制 %d 个散点图，已保存到: %s\n", plot_count, scatter_plot_dir))
    } else {
      cat("  警告：未找到生存分析结果文件，无法获取pearson_r\n")
    }
  } else {
    cat("  警告：未找到预测值表格文件\n")
  }
} else {
  cat(sprintf("  警告：未找到Filtered_Proteins_VarianceDecomposition文件: %s\n", variance_decomp_path))
}

cat("完成！\n")
