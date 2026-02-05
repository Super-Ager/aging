library(parallel)
library(dplyr)
library(survival)
library(lme4)  # For mixed-effects models
library(glmnet)  # For Lasso regression
library(rms)  # For ANOVA and R2 calculation
library(ggplot2)
library(gridExtra)
library(scales)

options(warn = -1)

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------
soma_path <- 'data/soma_df_Age_Sex_F_Dataset.csv'
anno_path <- 'data/annoinfo_df_SOMA.csv'
outcome_path <- 'data/CMCS_outcome_final.csv'

out_root <- 'results/004-longevity-target'
plot_dir <- file.path(out_root, 'plots')
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Data loading
# -------------------------------------------------------------------------
cat("Loading data...\n")
soma_df <- read.csv(soma_path, stringsAsFactors = FALSE) %>%
  filter(Visit %in% c(1,2,3,4)) %>%
  rename(visit = Visit)

anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
anno_map_symbol  <- setNames(anno_df$EntrezGeneSymbol, anno_df$AptName)
anno_map_uniprot <- setNames(anno_df$UniProt, anno_df$AptName)

# Load death information
outcome_df <- read.csv(outcome_path, stringsAsFactors = FALSE)

# Keep only IDs present in soma_df
outcome_df <- outcome_df %>%
  filter(ID %in% soma_df$ID)

# Convert death date
death_date_original <- outcome_df$death_Date
outcome_df$death_Date <- as.Date(outcome_df$death_Date)

# Check number of conversion failures
conversion_failed <- !is.na(death_date_original) & death_date_original != "" & is.na(outcome_df$death_Date)
if (sum(conversion_failed) > 0) {
  cat(sprintf("  Warning: %d death_Date conversions failed, trying other formats...\n", sum(conversion_failed)))
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

# Last update date
last_update_date <- as.Date("2024/8/21", format = "%Y/%m/%d")

seqid_cols <- grep("^seq\\.", names(soma_df), value = TRUE)
cat(sprintf("Found %d seqids\n", length(seqid_cols)))

# Get all IDs (for final table)
all_ids <- unique(soma_df$ID)
all_ids <- sort(all_ids)

# -------------------------------------------------------------------------
# Helper function: Calculate R2 for Cox model (based on partial likelihood ratio statistic)
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
# Helper function: Calculate C-index and its confidence interval
# -------------------------------------------------------------------------
calculate_cindex_ci <- function(cox_model) {
  tryCatch({
    cindex_obj <- concordance(cox_model)
    cindex <- as.numeric(cindex_obj$concordance)
    
    # Calculate confidence interval (using standard error)
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
# Single protein trajectory modeling and survival analysis
# -------------------------------------------------------------------------
analyze_protein_gap_survival <- function(seqid) {
  cat(sprintf("Analyzing %s...\n", seqid))
  
  # Prepare data
  df <- data.frame(
    idno  = soma_df$ID,
    age   = soma_df$age,
    sex   = soma_df$Sex_F,
    visit = factor(soma_df$visit, levels = c(1,2,3,4)),
    expr  = soma_df[[seqid]]
  ) %>% na.omit()
  
  if (nrow(df) < 30) return(NULL)
  
  # Log transformation
  df$log_expr <- log(df$expr + 1)
  
  # ±3 SD filtering
  mu  <- mean(df$log_expr)
  sdv <- sd(df$log_expr)
  df  <- df[df$log_expr >= mu - 3*sdv & df$log_expr <= mu + 3*sdv, ]
  if (nrow(df) < 30) return(NULL)
  
  # Batch effect correction
  lm_batch <- lm(log_expr ~ age + visit, data = df)
  visit_mat  <- model.matrix(lm_batch)[, grep("^visit", colnames(model.matrix(lm_batch))), drop = FALSE]
  visit_beta <- coef(lm_batch)[grep("^visit", names(coef(lm_batch)))]
  batch_eff  <- as.numeric(visit_mat %*% visit_beta)
  df$log_expr_corrected <- df$log_expr - batch_eff
  
  # =======================================================================
  # Step 1: Filter individuals with visit 1, 2, 3, get df_V123
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
  
  # Extract visit 1, 2, 3 data for these individuals
  df_V123 <- df %>%
    filter(idno %in% idno_V123 & visit %in% c(1, 2, 3)) %>%
    select(idno, visit, age, sex, log_expr_corrected) %>%
    arrange(idno, visit)
  
  if (nrow(df_V123) < 20) return(NULL)
  
  # =======================================================================
  # Step 2: Cross-sectional model (using visit 2 data)
  # =======================================================================
  # Extract visit 2 data for fitting cross-sectional model
  df_V2 <- df_V123 %>%
    filter(visit == 2) %>%
    select(idno, age, sex, log_expr_corrected)
  
  if (nrow(df_V2) < 10) return(NULL)
  
  # Step 2.1: Calculate mean and standard deviation of age using visit 2 (for standardizing visit 2 and visit 3)
  age_mean <- mean(df_V2$age, na.rm = TRUE)
  age_sd <- sd(df_V2$age, na.rm = TRUE)
  
  # Step 2.2: Standardize visit 2 age using the above mean and standard deviation
  df_V2$age_cs <- (df_V2$age - age_mean) / age_sd
  
  # Step 2.3: Calculate regression line using visit 2 data (standardized age and sex)
  crosssectional_model <- tryCatch({
    lm(log_expr_corrected ~ age_cs + sex, data = df_V2)
  }, error = function(e) {
    cat(sprintf("  Cross-sectional model fitting failed: %s\n", e$message))
    return(NULL)
  })
  
  if (is.null(crosssectional_model)) return(NULL)
  
  # Step 2.4: Calculate gap (residual) from regression line for each individual at visit 2
  df_V2$predicted_V2 <- predict(crosssectional_model, newdata = df_V2)
  df_V2$level_gap <- df_V2$log_expr_corrected - df_V2$predicted_V2
  level_gap_map <- setNames(df_V2$level_gap, df_V2$idno)
  
  # =======================================================================
  # Step 3: Prepare visit 3 data and features
  # =======================================================================
  # Extract visit 1, 2, 3 data for modeling and prediction
  df_V1 <- df_V123 %>%
    filter(visit == 1) %>%
    select(idno, visit1_level = log_expr_corrected)
  
  df_V2_data <- df_V123 %>%
    filter(visit == 2) %>%
    select(idno, visit2_level = log_expr_corrected)
  
  # Extract visit 3 data and standardize using visit 2 age mean and standard deviation
  df_V3 <- df_V123 %>%
    filter(visit == 3) %>%
    select(idno, age, sex, actual_visit3 = log_expr_corrected) %>%
    mutate(age_cs = (age - age_mean) / age_sd)  # Standardize visit 3 age using visit 2 parameters
  
  if (nrow(df_V3) < 20) return(NULL)
  
  # Merge all feature data (ensure only individuals with visit 1, 2, 3 are used)
  df_V3 <- df_V3 %>%
    left_join(df_V1, by = "idno") %>%
    left_join(df_V2_data, by = "idno") %>%
    filter(!is.na(visit1_level) & !is.na(visit2_level))
  
  if (nrow(df_V3) < 20) return(NULL)
  
  # Step 3.1: Calculate value on regression line for each individual at visit 3 (using standardized age and sex)
  df_V3$predicted_V3_cross <- predict(crosssectional_model, newdata = df_V3)
  
  # Step 3.2: Add visit 2 level gap as visit 3 prediction
  df_V3$level_gap <- level_gap_map[df_V3$idno]
  df_V3$predicted_visit3_cross <- df_V3$predicted_V3_cross + df_V3$level_gap
  
  # Calculate cross-sectional model prediction accuracy
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
  
  # Calculate cross-sectional model gap
  df_V3$gap_cross <- df_V3$actual_visit3 - df_V3$predicted_visit3_cross
  
  # =======================================================================
  # Step 4: Lasso prediction model (using visit 1 level, visit 2 level, and sex)
  # =======================================================================
  # Prepare feature matrix and response variable (only use individuals in df_V3, ensure visit 1, 2, 3)
  # Features: visit1_level, visit2_level, sex
  # Output: actual_visit3
  X <- as.matrix(df_V3[, c("visit1_level", "visit2_level", "sex")])
  y <- df_V3$actual_visit3
  n_individuals <- nrow(df_V3)
  
  # 5-fold cross-validation, split by individual
  n_folds <- 5
  fold_size <- ceiling(n_individuals / n_folds)
  
  # Randomly shuffle individual order (but maintain reproducibility)
  set.seed(123)
  shuffled_indices <- sample(1:n_individuals)
  
  # Initialize predictions
  predicted_visit3_lasso <- rep(NA, n_individuals)
  
  # Cross-validation for each fold
  for (fold in 1:n_folds) {
    # Determine test set indices
    start_idx <- (fold - 1) * fold_size + 1
    end_idx <- min(fold * fold_size, n_individuals)
    test_indices <- shuffled_indices[start_idx:end_idx]
    train_indices <- shuffled_indices[-test_indices]
    
    if (length(train_indices) < 10 || length(test_indices) < 1) next
    
    # Prepare training and test sets
    X_train <- X[train_indices, , drop = FALSE]
    y_train <- y[train_indices]
    X_test <- X[test_indices, , drop = FALSE]
    
    # Use Lasso for regression (alpha=1, glmnet will automatically standardize)
    tryCatch({
      # Use Lasso (alpha=1) for cross-validation
      cv_fit <- cv.glmnet(X_train, y_train, alpha = 1, nfolds = 5)
      
      # Predict (using lambda.1se for better generalization)
      predictions <- predict(cv_fit, newx = X_test, s = "lambda.1se")
      
      # Save predictions (test_indices are already original data indices)
      if (length(test_indices) > 0 && length(predictions) > 0) {
        predicted_visit3_lasso[test_indices] <- as.numeric(predictions)
      }
      
    }, error = function(e) {
      cat(sprintf("  Fold %d prediction failed: %s\n", fold, e$message))
    })
  }
  
  # Add predictions to df_V3 (df_V3 already ensures individuals with visit 1, 2, 3)
  df_V3$predicted_visit3_linear <- predicted_visit3_lasso
  
  # Keep only individuals with complete predictions (both cross-sectional and Lasso models have predictions)
  df_V3 <- df_V3 %>%
    filter(!is.na(predicted_visit3_linear) & !is.na(predicted_visit3_cross) & 
           !is.na(actual_visit3))
  
  if (nrow(df_V3) < 10) return(NULL)
  
  # Calculate Lasso model prediction accuracy
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
  
  # Calculate Lasso model gap
  df_V3$gap_linear <- df_V3$actual_visit3 - df_V3$predicted_visit3_linear
  
  # Use Lasso model predictions as predicted_visit3_mixed (for subsequent survival analysis)
  df_V3$predicted_visit3_mixed <- df_V3$predicted_visit3_linear
  df_V3$gap <- df_V3$gap_linear
  
  # Standardize predicted_visit3 and gap (for survival analysis)
  df_V3$predicted_visit3_sd <- as.numeric(scale(df_V3$predicted_visit3_mixed))
  df_V3$gap_sd <- as.numeric(scale(df_V3$gap))
  
  # =======================================================================
  # Merge death information (df_V3 already ensures individuals with visit 1, 2, 3)
  # =======================================================================
  df_survival <- df_V3 %>%
    left_join(outcome_df %>% select(ID, death_Date, death), by = c("idno" = "ID"))
  
  if (nrow(df_survival) < 20 || sum(df_survival$death == 1, na.rm = TRUE) < 5) return(NULL)
  
  # Prepare survival data
  visit3_baseline_date <- as.Date("2013-01-01")  # Adjust according to actual situation
  
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
  
  # Prepare covariates
  df_survival$age_visit3_sd <- as.numeric(scale(df_survival$age))
  df_survival$sex_factor <- factor(df_survival$sex, levels = c(0,1), labels = c("Male","Female"))
  
  # Initialize results
  result <- list(
    seqid = seqid,
    n_samples = nrow(df_survival),
    n_death = sum(df_survival$death == 1, na.rm = TRUE),
    # Prediction accuracy (cross-sectional model)
    cross_cor_r = cor_cross_r,
    cross_cor_p = cor_cross_p,
    cross_r2 = r2_cross,
    # Prediction accuracy (linear model)
    linear_cor_r = cor_linear_r,
    linear_cor_p = cor_linear_p,
    linear_r2 = r2_linear,
    # Model results
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
  
  # Fit models and calculate metrics
  tryCatch({
    # Ensure data has no NA values
    df_model <- df_survival %>%
      select(time, event, age_visit3_sd, sex_factor, predicted_visit3_sd, gap_sd) %>%
      filter(!is.na(time) & !is.na(event) & !is.na(age_visit3_sd) & 
             !is.na(sex_factor) & !is.na(predicted_visit3_sd))
    
    if (nrow(df_model) < 20 || sum(df_model$event) < 5) {
      stop("Insufficient data or too few events")
    }
    
    # Base model: age + sex + predicted_visit3
    base_model <- coxph(
      Surv(time, event) ~ age_visit3_sd + sex_factor + predicted_visit3_sd,
      data = df_model
    )
    
    # Full model: age + sex + predicted_visit3 + gap
    full_model <- coxph(
      Surv(time, event) ~ age_visit3_sd + sex_factor + predicted_visit3_sd + gap_sd,
      data = df_model
    )
    
    # Calculate R²
    result$base_model_r2 <- calculate_cox_r2(base_model)
    result$full_model_r2 <- calculate_cox_r2(full_model)
    
    # Calculate C-index and its confidence interval
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
    
    # Extract HR of gap and its confidence interval
    gap_coef <- summary(full_model)$coefficients["gap_sd", ]
    result$gap_hr <- exp(gap_coef["coef"])
    result$gap_hr_ci_lower <- exp(gap_coef["coef"] - 1.96 * gap_coef["se(coef)"])
    result$gap_hr_ci_upper <- exp(gap_coef["coef"] + 1.96 * gap_coef["se(coef)"])
    result$gap_p <- gap_coef["Pr(>|z|)"]
    
  }, error = function(e) {
    cat(sprintf("  Model fitting error: %s\n", e$message))
  })
  
  # Prepare prediction table data
  prediction_df <- data.frame(
    idno = all_ids,
    stringsAsFactors = FALSE
  )
  
  # Add three columns for current protein
  col_actual <- paste0(seqid, "_actual")
  col_mixed <- paste0(seqid, "_mixed")
  col_cross <- paste0(seqid, "_cross")
  
  prediction_df[[col_actual]] <- NA
  prediction_df[[col_mixed]] <- NA
  prediction_df[[col_cross]] <- NA
  
  # Fill individuals with data (df_V3 already ensures visit 1, 2, 3 and complete predictions)
  for (i in seq_len(nrow(df_V3))) {
    idx <- which(prediction_df$idno == df_V3$idno[i])
    if (length(idx) == 1) {
      prediction_df[[col_actual]][idx] <- df_V3$actual_visit3[i]
      prediction_df[[col_mixed]][idx] <- df_V3$predicted_visit3_mixed[i]
      prediction_df[[col_cross]][idx] <- df_V3$predicted_visit3_cross[i]
    }
  }
  
  # Return results and prediction table
  return(list(
    result = data.frame(result, stringsAsFactors = FALSE),
    prediction = prediction_df
  ))
}

# -------------------------------------------------------------------------
# Main program: Parallel analysis
# -------------------------------------------------------------------------
cat("Starting trajectory modeling and survival analysis...\n")
cat(sprintf("Total of %d seqids to analyze\n", length(seqid_cols)))

# Test code: analyze only specified seqids
# seqid_cols <- c("seq.4374.45")
# seqid_cols <- seqid_cols[1:10]

# Set parallel computation
n_cores <- 14
cat(sprintf("Using %d cores for parallel computation...\n", n_cores))

# Use mclapply for parallel computation
cat("Starting parallel analysis...\n")
result_lists <- mclapply(seq_along(seqid_cols), function(i) {
  seqid <- seqid_cols[i]
  
  result_list <- tryCatch({
    analyze_protein_gap_survival(seqid)
  }, error = function(e) {
    return(NULL)
  })
  
  # Output progress every 100
  if (i %% 100 == 0) {
    cat(sprintf("Progress: %d/%d (%.1f%%)\n", i, length(seqid_cols), 100 * i / length(seqid_cols)))
  }
  
  return(result_list)
}, mc.cores = n_cores)

# Process results
res <- list()
prediction_tables <- list()

for (i in seq_along(result_lists)) {
  result_list <- result_lists[[i]]
  
  if (is.null(result_list)) {
    next
  }
  
  # Check if results are valid
  if (is.null(result_list$result) || !is.data.frame(result_list$result)) {
    next
  }
  
  # Save results
  res[[i]] <- result_list$result
  
  # Save prediction table (initialize first time, merge subsequent times)
  if (!is.null(result_list$prediction) && is.data.frame(result_list$prediction)) {
    if (length(prediction_tables) == 0 || is.null(prediction_tables[[1]])) {
      prediction_tables[[1]] <- result_list$prediction
    } else {
      # Merge new columns
      new_cols <- result_list$prediction[, -1, drop = FALSE]  # Exclude idno column
      if (ncol(new_cols) > 0) {
        prediction_tables[[1]] <- cbind(prediction_tables[[1]], new_cols)
      }
    }
  }
}

if (length(res) == 0) {
  cat("No valid results\n")
  stop("No valid results")
}

# Remove NULL results
res <- res[!sapply(res, is.null)]

survival_df <- do.call(rbind, res)
cat(sprintf("\nFinal saved results: %d rows x %d columns\n", nrow(survival_df), ncol(survival_df)))

# Save results
output_path <- file.path(out_root, "protein_gap_survival_results.csv")
write.csv(survival_df, output_path, row.names = FALSE)
cat(sprintf("Results saved to: %s\n", output_path))

# Save prediction table
prediction_path <- file.path(out_root, "protein_visit3_predictions.csv")
if (length(prediction_tables) > 0 && !is.null(prediction_tables[[1]])) {
  prediction_final <- prediction_tables[[1]]
  
  # Replace NA with NaN (for Python compatibility)
  for (col in names(prediction_final)) {
    if (col != "idno") {
      prediction_final[[col]][is.na(prediction_final[[col]])] <- NaN
    }
  }
  
  write.csv(prediction_final, prediction_path, row.names = FALSE)
  cat(sprintf("Prediction table saved to: %s\n", prediction_path))
}

# -------------------------------------------------------------------------
# Plot scatter plots of linear model predictions vs actual values (for Filtered_Proteins)
# -------------------------------------------------------------------------
cat("\n========================================\n")
cat("Plotting linear model prediction accuracy scatter plots\n")
cat("========================================\n\n")

# Read Filtered_Proteins_VarianceDecomposition
variance_decomp_path <- 'results/001-cross-sectional-aging-marker-detection/Filtered_Proteins_VarianceDecomposition.csv'
if (file.exists(variance_decomp_path)) {
  cat("Reading Filtered_Proteins_VarianceDecomposition...\n")
  variance_df <- read.csv(variance_decomp_path, stringsAsFactors = FALSE)
  
  # Check seq_id column name (may be seq_id or seqid)
  seq_id_col <- if ("seq_id" %in% names(variance_df)) "seq_id" else "seqid"
  filtered_seqids <- unique(variance_df[[seq_id_col]])
  cat(sprintf("  Found %d filtered proteins\n", length(filtered_seqids)))
  
  # Create scatter plot output directory
  scatter_plot_dir <- file.path(out_root, "linear_prediction_scatter_plots")
  dir.create(scatter_plot_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("  Scatter plots will be saved to: %s\n", scatter_plot_dir))
  
  # Read prediction table
  if (file.exists(prediction_path)) {
    cat("Reading prediction table...\n")
    prediction_df <- read.csv(prediction_path, stringsAsFactors = FALSE)
    cat(sprintf("  Prediction table: %d rows x %d columns\n", nrow(prediction_df), ncol(prediction_df)))
    
    # Read survival analysis results to get pearson_r
    if (file.exists(output_path)) {
      cat("Reading survival analysis results...\n")
      survival_results <- read.csv(output_path, stringsAsFactors = FALSE)
      
      # Plot scatter plots for each filtered protein
      plot_count <- 0
      for (seqid in filtered_seqids) {
        # Check if corresponding columns exist
        col_actual <- paste0(seqid, "_actual")
        col_mixed <- paste0(seqid, "_mixed")
        
        if (col_actual %in% names(prediction_df) && col_mixed %in% names(prediction_df)) {
          # Extract actual and predicted values
          plot_data <- prediction_df %>%
            select(idno, actual = all_of(col_actual), predicted = all_of(col_mixed)) %>%
            filter(!is.na(actual) & !is.na(predicted) & 
                   !is.nan(actual) & !is.nan(predicted) &
                   is.finite(actual) & is.finite(predicted))
          
          if (nrow(plot_data) >= 10) {
            # Calculate Pearson correlation coefficient
            pearson_r <- tryCatch({
              cor_result <- cor.test(plot_data$actual, plot_data$predicted, use = "complete.obs")
              cor_result$estimate
            }, error = function(e) {
              NA
            })
            
            # If available in survival_results, use that value
            if (seqid %in% survival_results$seqid) {
              pearson_r <- survival_results$linear_cor_r[survival_results$seqid == seqid][1]
            }
            
            # Get gene_symbol
            gene_symbol <- if (seqid %in% names(anno_map_symbol)) {
              anno_map_symbol[seqid]
            } else {
              seqid
            }
            
            # If gene_symbol is empty or NA, use seqid
            if (is.na(gene_symbol) || gene_symbol == "" || is.null(gene_symbol)) {
              gene_symbol <- seqid
            }
            
            # Plot scatter plot
            scatter_plot <- ggplot(plot_data, aes(x = actual, y = predicted)) +
              geom_point(color = "#E74C3C", alpha = 0.6, size = 1.5) +
              geom_abline(intercept = 0, slope = 1, color = "#2C3E50", linetype = "dashed", linewidth = 0.8, alpha = 0.8) +
              # Add pearson_r annotation
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
            
            # Save plot
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
              cat(sprintf("  Plotted %d scatter plots...\n", plot_count))
            }
          }
        }
      }
      
      cat(sprintf("\n✅ Total of %d scatter plots plotted, saved to: %s\n", plot_count, scatter_plot_dir))
    } else {
      cat("  Warning: Survival analysis results file not found, cannot get pearson_r\n")
    }
  } else {
    cat("  Warning: Prediction table file not found\n")
  }
} else {
  cat(sprintf("  Warning: Filtered_Proteins_VarianceDecomposition file not found: %s\n", variance_decomp_path))
}

cat("Completed!\n")
