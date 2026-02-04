library(parallel)
library(dplyr)
library(ggplot2)
library(lme4)
library(MuMIn)
library(lmerTest)

options(warn = -1)

# -------------------------------------------------------------------------
# Paths
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
# Data loading + visit renaming
# -------------------------------------------------------------------------
soma_df <- read.csv(soma_path, stringsAsFactors = FALSE) %>%
  filter(Dataset %in% c(1,2,3,4)) %>%
  rename(visit = Dataset)

anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)

anno_map_symbol  <- setNames(anno_df$EntrezGeneSymbol, anno_df$AptName)
anno_map_uniprot <- setNames(anno_df$UniProt, anno_df$AptName)

seqid_cols <- grep("^seq\\.", names(soma_df), value = TRUE)

# -------------------------------------------------------------------------
# Single protein analysis + plotting (grouped by sex)
# -------------------------------------------------------------------------
analyze_protein_by_sex <- function(seqid, sex_value, sex_label) {
  tryCatch({
    cat(sprintf("%s - %s\n", seqid, sex_label))
    if (!seqid %in% names(soma_df)) return(NULL)

  # Filter data by sex
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

  # ---------------- downstream model (remove sex fixed effect) ----------------
  df$z_expr <- as.numeric(scale(df$log_expr_corrected))
  df$age_cs <- as.numeric(scale(df$age))

  # Model: z_expr ~ age_cs + (1 | idno) (no sex)
  m <- lmer(z_expr ~ age_cs + (1 | idno), data = df)
  sm <- summary(m)

  coef_age <- sm$coefficients["age_cs", ]

  # Correct variance decomposition method (referring to Nakagawa-Schielzeth method)
  # ============================================================
  # 1. Fixed effects variance: Var(Xβ) —— only age contribution
  # ============================================================
  
  # Fixed effects design matrix
  X <- getME(m, "X")
  
  # Fixed effects coefficients
  beta <- fixef(m)
  
  # Fixed effects linear predictor
  eta_fixed <- as.vector(X %*% beta)
  
  # Fixed effects variance (σ²_fixed)
  var_fixed_total <- var(eta_fixed)
  
  # ============================================================
  # 2. Age variance contribution
  # ============================================================
  
  # Age contribution: only intercept and age
  beta_intercept <- beta["(Intercept)"]
  beta_age <- beta["age_cs"]
  
  # Age contribution: only intercept and age
  X_age_only <- cbind(1, df$age_cs)  # intercept, age
  eta_age_only <- as.vector(X_age_only %*% c(beta_intercept, beta_age))
  var_fixed_age_only <- var(eta_age_only)
  
  # Fixed effects predictions with only intercept (baseline)
  eta_intercept_only <- rep(beta_intercept, nrow(df))
  var_fixed_intercept_only <- var(eta_intercept_only)  # Should be 0 (constant)
  
  # Incremental variance contribution of age
  var_age <- var_fixed_age_only - var_fixed_intercept_only
  
  # Since there is only one fixed effect (age), var_age should equal var_fixed_total
  # But for numerical stability, we use the calculated var_age
  var_age <- max(var_age, var_fixed_total)
  
  # ============================================================
  # 3. Random effects variance
  # ============================================================
  
  vc <- VarCorr(m)
  var_idno <- as.numeric(attr(vc$idno, "stddev"))^2
  
  # ============================================================
  # 4. Residual variance
  # ============================================================
  
  var_res <- sigma(m)^2
  
  # ============================================================
  # 5. Total variance (Nakagawa–Schielzeth definition)
  # ============================================================
  
  var_tot <- var_fixed_total + var_idno + var_res
  
  # ============================================================
  # 6. Calculate R² for output
  # ============================================================
  
  R2_full <- r.squaredGLMM(m)
  R2m_full <- R2_full[1, "R2m"]
  R2c_full <- R2_full[1, "R2c"]
  
  # Verification: prop_age should equal var_fixed_total / var_tot = R2m
  # Verification: prop_age + prop_idno should equal (var_fixed_total + var_idno) / var_tot = R2c

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
    cat(sprintf("Error: %s - %s: %s\n", seqid, sex_label, conditionMessage(e)))
    return(NULL)
  })
}

# -------------------------------------------------------------------------
# Parallel execution (process male and female separately)
# -------------------------------------------------------------------------
n_cores <- max(1, min(detectCores() - 2, 8))  # Limit max cores to avoid resource issues

# Process male (sex = 0)
cat("\n========== Starting to process male data ==========\n")
cat(sprintf("Using %d cores for parallel processing\n", n_cores))
res_male <- tryCatch({
  mclapply(seqid_cols, function(seqid) {
    analyze_protein_by_sex(seqid, sex_value = 0, sex_label = "Male")
  }, mc.cores = n_cores, mc.preschedule = FALSE)  # mc.preschedule = FALSE can reduce SIGPIPE errors
}, error = function(e) {
  cat(sprintf("Parallel processing error: %s\n", conditionMessage(e)))
  cat("Trying single-core processing...\n")
  lapply(seqid_cols, function(seqid) {
    analyze_protein_by_sex(seqid, sex_value = 0, sex_label = "Male")
  })
})
res_male <- res_male[!sapply(res_male, is.null)]
if (length(res_male) > 0) {
  df_male <- do.call(rbind, res_male)
} else {
  df_male <- data.frame()
  cat("Warning: No successfully processed male data\n")
}

# Process female (sex = 1)
cat("\n========== Starting to process female data ==========\n")
cat(sprintf("Using %d cores for parallel processing\n", n_cores))
res_female <- tryCatch({
  mclapply(seqid_cols, function(seqid) {
    analyze_protein_by_sex(seqid, sex_value = 1, sex_label = "Female")
  }, mc.cores = n_cores, mc.preschedule = FALSE)
}, error = function(e) {
  cat(sprintf("Parallel processing error: %s\n", conditionMessage(e)))
  cat("Trying single-core processing...\n")
  lapply(seqid_cols, function(seqid) {
    analyze_protein_by_sex(seqid, sex_value = 1, sex_label = "Female")
  })
})
res_female <- res_female[!sapply(res_female, is.null)]
if (length(res_female) > 0) {
  df_female <- do.call(rbind, res_female)
} else {
  df_female <- data.frame()
  cat("Warning: No successfully processed female data\n")
}

# -------------------------------------------------------------------------
# Save results
# -------------------------------------------------------------------------
if (nrow(df_male) > 0) {
  write.csv(
    df_male,
    file.path(out_root, "combined_protein_analysis_results_male.csv"),
    row.names = FALSE
  )
} else {
  cat("Warning: Male data is empty, file not saved\n")
}

if (nrow(df_female) > 0) {
  write.csv(
    df_female,
    file.path(out_root, "combined_protein_analysis_results_female.csv"),
    row.names = FALSE
  )
} else {
  cat("Warning: Female data is empty, file not saved\n")
}

cat("\n========== Completed ==========\n")
cat(sprintf("Male results: %d proteins\n", nrow(df_male)))
cat(sprintf("Female results: %d proteins\n", nrow(df_female)))
cat(sprintf("Results saved to: %s\n", out_root))

