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

out_root <- '/home/ww/Project/Protein-Composition/result/1227-result'
plot_dir <- file.path(out_root, 'png-protein-visit-batch')

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

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
# Single protein analysis + plotting
# -------------------------------------------------------------------------
analyze_protein <- function(seqid) {
  cat(seqid, "\n")
  if (!seqid %in% names(soma_df)) return(NULL)

  df <- data.frame(
    idno  = soma_df$ID,
    age   = soma_df$age,
    sex   = soma_df$Sex_F,
    visit = factor(soma_df$visit, levels = c(1,2,3,4)),
    expr  = soma_df[[seqid]]
  ) %>% na.omit()

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
  plot_df <- rbind(
    df %>% transmute(age, expr = log_expr, visit, stage = "Pre-normalization"),
    df %>% transmute(age, expr = log_expr_corrected, visit, stage = "Post-normalization")
  )

  y_lim <- range(plot_df$expr)

  p <- ggplot(plot_df, aes(age, expr, color = visit)) +
    geom_point(alpha = 0.25, size = 1) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, formula = y ~ x) +
    facet_wrap(~stage, nrow = 1) +
    scale_y_continuous(limits = y_lim) +
    labs(
      x = "Age",
      y = seqid,
      title = paste0(seqid, " (", anno_map_symbol[seqid], ", ", anno_map_uniprot[seqid], ")")
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank()
    )

  ggsave(
    filename = file.path(plot_dir, paste0(seqid, ".png")),
    plot = p,
    width = 9,
    height = 4,
    dpi = 150
  )

  # ---------------- downstream model ----------------
  df$z_expr <- as.numeric(scale(df$log_expr_corrected))
  df$age_cs <- as.numeric(scale(df$age))
  df$sex_factor <- factor(df$sex, levels = c(0,1), labels = c("Male","Female"))

  m <- lmer(z_expr ~ age_cs + sex_factor + (1 | idno), data = df)
  sm <- summary(m)

  coef_age <- sm$coefficients["age_cs", ]
  coef_sex <- sm$coefficients["sex_factorFemale", ]

  # Correct variance decomposition method (referring to Nakagawa-Schielzeth method)
  # ============================================================
  # 1. Fixed effects variance: Var(Xβ) —— consistent with r.squaredGLMM
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
  # 2. Decompose variance contributions of age and sex separately
  # ============================================================
  
  # Calculate fixed effects predictions with only age (sex kept at baseline: Male=0)
  beta_intercept <- beta["(Intercept)"]
  beta_age <- beta["age_cs"]
  beta_sex <- beta["sex_factorFemale"]
  
  # Age contribution: only intercept and age, sex=0
  X_age_only <- cbind(1, df$age_cs, 0)  # intercept, age, sex=0 (Male)
  eta_age_only <- as.vector(X_age_only %*% c(beta_intercept, beta_age, 0))
  var_fixed_age_only <- var(eta_age_only)
  
  # Sex contribution: only intercept and sex, age kept at mean (age_cs=0, since centered)
  sex_indicator <- as.numeric(df$sex_factor == "Female")
  X_sex_only <- cbind(1, 0, sex_indicator)  # intercept, age=0, sex
  eta_sex_only <- as.vector(X_sex_only %*% c(beta_intercept, 0, beta_sex))
  var_fixed_sex_only <- var(eta_sex_only)
  
  # Fixed effects predictions with only intercept (baseline)
  eta_intercept_only <- rep(beta_intercept, nrow(df))
  var_fixed_intercept_only <- var(eta_intercept_only)  # Should be 0 (constant)
  
  # Calculate incremental variance contributions of age and sex
  var_age_increment <- var_fixed_age_only - var_fixed_intercept_only
  var_sex_increment <- var_fixed_sex_only - var_fixed_intercept_only
  
  # Since age and sex may interact, their incremental sum may not equal total fixed effects variance
  # Scale proportionally so that var_age + var_sex = var_fixed_total
  if (var_age_increment + var_sex_increment > 0) {
    scale_factor <- var_fixed_total / (var_age_increment + var_sex_increment)
    var_age <- var_age_increment * scale_factor
    var_sex <- var_sex_increment * scale_factor
  } else {
    # If both are 0, distribute equally
    var_age <- var_fixed_total * 0.5
    var_sex <- var_fixed_total * 0.5
  }
  
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
  
  # Verification: prop_age + prop_sex should equal var_fixed_total / var_tot = R2m
  # Verification: prop_age + prop_sex + prop_idno should equal (var_fixed_total + var_idno) / var_tot = R2c

  data.frame(
    seqid = seqid,
    gene_symbol = anno_map_symbol[seqid],
    uniprot = anno_map_uniprot[seqid],

    n_samples = nrow(df),
    n_individuals = length(unique(df$idno)),

    age_Estimate = coef_age["Estimate"],
    age_StdError = coef_age["Std. Error"],
    age_tvalue   = coef_age["t value"],
    age_pvalue   = coef_age["Pr(>|t|)"],

    sex_Estimate = coef_sex["Estimate"],
    sex_StdError = coef_sex["Std. Error"],
    sex_tvalue   = coef_sex["t value"],
    sex_pvalue   = coef_sex["Pr(>|t|)"],

    variance_age = var_age,
    variance_sex = var_sex,
    variance_idno = var_idno,
    variance_residual = var_res,
    variance_total = var_tot,

    prop_age = var_age / var_tot,
    prop_sex = var_sex / var_tot,
    prop_idno = var_idno / var_tot,
    prop_residual = var_res / var_tot,

    R2m = R2m_full,
    R2c = R2_full[1, "R2c"],

    stringsAsFactors = FALSE
  )
}

# -------------------------------------------------------------------------
# Parallel execution (can limit to first few for testing)
# -------------------------------------------------------------------------
n_cores <- max(1, detectCores() - 2)
# seqid_cols <- seqid_cols[1:10]  # For testing
res <- mclapply(seqid_cols, analyze_protein, mc.cores = n_cores)
res <- res[!sapply(res, is.null)]

final_df <- do.call(rbind, res)

write.csv(
  final_df,
  file.path(out_root, "combined_protein_analysis_results.csv"),
  row.names = FALSE
)

cat("Finished. Results saved to:\n", out_root)

