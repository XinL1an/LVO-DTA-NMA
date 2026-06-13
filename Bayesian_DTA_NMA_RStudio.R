# Bayesian diagnostic test accuracy network meta-analysis
# RStudio-ready standalone script

# How to use:
# 1. Keep this script in the same folder as the data/ directory.
# 2. Open this script in RStudio.
# 3. Install JAGS and the required R package R2jags.
# 4. Adjust the analysis switches below if needed.
# 5. Click Source.

# -----------------------------
# User settings
# -----------------------------

RUN_PRIMARY_ANALYSIS <- TRUE
RUN_META_REGRESSION <- TRUE
RUN_SENSITIVITY_ANALYSES <- TRUE
RUN_PRIOR_SENSITIVITY <- TRUE

N_CHAINS <- 4
N_ITER <- 210000
N_BURNIN <- 10000
N_THIN <- 10
META_N_ITER <- 50000
META_N_BURNIN <- 10000
META_N_THIN <- 40
RUN_LABEL <- format(Sys.time(), "run_%Y%m%d_%H%M%S")

SEED <- 42

NODE_NAMES <- c(
  "Expert Readers",
  "Non-Expert Readers",
  "Unimodal Imaging AI",
  "Clinically-Informed Multimodal AI"
)

# -----------------------------
# Setup
# -----------------------------

required_packages <- c("R2jags")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(R2jags)
})

get_project_root <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    if (!is.null(active_path) && nzchar(active_path)) {
      return(dirname(normalizePath(active_path, winslash = "/")))
    }
  }

  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- args_all[grepl("^--file=", args_all)]
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

PROJECT_ROOT <- get_project_root()
DATA_DIR <- file.path(PROJECT_ROOT, "data")
OUTPUT_ROOT <- file.path(PROJECT_ROOT, "outputs", RUN_LABEL)
dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

cat("Project root:", PROJECT_ROOT, "\n")
cat("Data directory:", DATA_DIR, "\n")
cat("Output directory:", OUTPUT_ROOT, "\n")

# -----------------------------
# Core functions
# -----------------------------

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

read_dta_input <- function(input_csv) {
  if (!file.exists(input_csv)) {
    stop("Missing input CSV: ", input_csv)
  }

  dta_data <- read.csv(input_csv, stringsAsFactors = FALSE, na.strings = c("", "NA"))

  required_cols <- c("s", "algo", "tp", "fp", "fn", "tn")
  missing_cols <- setdiff(required_cols, names(dta_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", input_csv, ": ", paste(missing_cols, collapse = ", "))
  }

  dta_data <- dta_data[!is.na(dta_data$s), , drop = FALSE]

  numeric_cols <- intersect(
    c("tp", "fp", "fn", "tn", "Cov_Thickness", "Cov_M2", "Cov_Validation", "Cov_Commercial"),
    names(dta_data)
  )
  for (col in numeric_cols) {
    dta_data[[col]] <- as.numeric(dta_data[[col]])
  }

  algo_values <- unique(dta_data$algo)
  algo_numeric <- suppressWarnings(as.numeric(algo_values))
  if (all(!is.na(algo_numeric))) {
    algo_levels <- as.character(sort(algo_numeric))
  } else {
    algo_levels <- sort(as.character(algo_values))
  }

  dta_data$study_id <- dta_data$s
  dta_data$algo_id <- dta_data$algo
  dta_data$algo <- as.numeric(factor(as.character(dta_data$algo), levels = algo_levels))
  dta_data$s <- as.numeric(as.factor(dta_data$s))
  dta_data$pos <- dta_data$tp + dta_data$fn
  dta_data$neg <- dta_data$fp + dta_data$tn

  dta_data <- dta_data[order(dta_data$s, dta_data$algo), , drop = FALSE]
  dta_data$is_baseline_arm <- ave(dta_data$algo, dta_data$s, FUN = seq_along) == 1
  dta_data$is_baseline_arm <- as.integer(dta_data$is_baseline_arm)
  rownames(dta_data) <- NULL

  dta_data
}

build_jags_data <- function(dta_data) {
  list(
    nalgo = length(unique(dta_data$algo)),
    ns = max(dta_data$s),
    nObs = nrow(dta_data),
    s = dta_data$s,
    algo = dta_data$algo,
    tp = dta_data$tp,
    tn = dta_data$tn,
    pos = dta_data$pos,
    neg = dta_data$neg,
    is_baseline_arm = dta_data$is_baseline_arm
  )
}

heterogeneity_precision_prior_block <- function(prior) {
  if (prior == "uniform") {
    return(paste(
      "  tau.study.sens <- pow(sd.study.sens, -2)",
      "  sd.study.sens ~ dunif(0, 2)",
      "  tau.study.spec <- pow(sd.study.spec, -2)",
      "  sd.study.spec ~ dunif(0, 2)",
      "  tau.rel.sens <- pow(sd.rel.sens, -2)",
      "  sd.rel.sens ~ dunif(0, 2)",
      "  tau.rel.spec <- pow(sd.rel.spec, -2)",
      "  sd.rel.spec ~ dunif(0, 2)",
      "  sd.se ~ dunif(0, 2)",
      "  sd.sp ~ dunif(0, 2)",
      sep = "\n"
    ))
  }

  if (prior == "halfnormal") {
    return(paste(
      "  tau.study.sens <- pow(sd.study.sens, -2)",
      "  tau.study.spec <- pow(sd.study.spec, -2)",
      "  tau.rel.sens <- pow(sd.rel.sens, -2)",
      "  tau.rel.spec <- pow(sd.rel.spec, -2)",
      "  sd.study.sens ~ dnorm(0, 1) T(0,)",
      "  sd.study.spec ~ dnorm(0, 1) T(0,)",
      "  sd.rel.sens ~ dnorm(0, 1) T(0,)",
      "  sd.rel.spec ~ dnorm(0, 1) T(0,)",
      "  sd.se ~ dnorm(0, 1) T(0,)",
      "  sd.sp ~ dnorm(0, 1) T(0,)",
      sep = "\n"
    ))
  }

  stop("Unknown heterogeneity prior: ", prior)
}

make_dta_nma_model <- function(heterogeneity_prior = "uniform") {
  prior_block <- heterogeneity_precision_prior_block(heterogeneity_prior)

  paste0(
    "model {\n",
    "  for(i in 1:nObs){\n",
    "    tp[i] ~ dbin(pi[i,1], pos[i])\n",
    "    tn[i] ~ dbin(pi[i,2], neg[i])\n",
    "    logit(pi[i,1]) <- mu[i,1]\n",
    "    logit(pi[i,2]) <- mu[i,2]\n",
    "    MU[i,1] <- algo.sens[algo[i]] + study.re.sens[s[i]] + (1 - is_baseline_arm[i]) * relative.re.sens[s[i], algo[i]]\n",
    "    MU[i,2] <- algo.spec[algo[i]] + study.re.spec[s[i]] + (1 - is_baseline_arm[i]) * relative.re.spec[s[i], algo[i]]\n",
    "    mu[i,1:2] ~ dmnorm(MU[i,], prec[,])\n",
    "  }\n",
    "  for(j in 1:nalgo){\n",
    "    algo.sens[j] ~ dnorm(0, 0.01)\n",
    "    algo.spec[j] ~ dnorm(0, 0.01)\n",
    "    sens[j] <- exp(algo.sens[j])/(1+exp(algo.sens[j]))\n",
    "    spec[j] <- exp(algo.spec[j])/(1+exp(algo.spec[j]))\n",
    "    DOR[j] <- exp(algo.sens[j] + algo.spec[j])\n",
    "  }\n",
    "  for(k in 1:ns){\n",
    "    study.re.sens[k] ~ dnorm(0, tau.study.sens)\n",
    "    study.re.spec[k] ~ dnorm(0, tau.study.spec)\n",
    "    for(l in 1:nalgo){\n",
    "      relative.re.sens[k,l] ~ dnorm(0, tau.rel.sens)\n",
    "      relative.re.spec[k,l] ~ dnorm(0, tau.rel.spec)\n",
    "    }\n",
    "  }\n",
    prior_block, "\n",
    "  rho ~ dunif(-0.99, 0.99)\n",
    "  var.se <- sd.se * sd.se\n",
    "  var.sp <- sd.sp * sd.sp\n",
    "  covar <- rho * sd.se * sd.sp\n",
    "  det <- var.se * var.sp - covar * covar\n",
    "  prec[1,1] <- var.sp / det\n",
    "  prec[2,2] <- var.se / det\n",
    "  prec[1,2] <- -covar / det\n",
    "  prec[2,1] <- prec[1,2]\n",
    "  tau.sq.sens <- var.se\n",
    "  tau.sens <- sd.se\n",
    "  tau.sq.spec <- var.sp\n",
    "  tau.spec <- sd.sp\n",
    "  tau.sq.logDOR <- var.se + var.sp + 2 * covar\n",
    "  tau.logDOR <- sqrt(max(0, tau.sq.logDOR))\n",
    "  ranks_sens <- rank(algo.sens[])\n",
    "  ranks_spec <- rank(algo.spec[])\n",
    "  ranks_DOR <- rank(DOR[])\n",
    "  for(l in 1:nalgo){\n",
    "    rksens[l] <- nalgo + 1 - ranks_sens[l]\n",
    "    rkspec[l] <- nalgo + 1 - ranks_spec[l]\n",
    "    rkDOR[l] <- nalgo + 1 - ranks_DOR[l]\n",
    "  }\n",
    "}\n"
  )
}

default_inits <- function(jags_data) {
  function() {
    list(
      algo.sens = rep(0, jags_data$nalgo),
      algo.spec = rep(0, jags_data$nalgo),
      sd.study.sens = runif(1, 0.1, 1),
      sd.study.spec = runif(1, 0.1, 1),
      sd.rel.sens = runif(1, 0.1, 1),
      sd.rel.spec = runif(1, 0.1, 1),
      sd.se = runif(1, 0.1, 1),
      sd.sp = runif(1, 0.1, 1),
      rho = runif(1, -0.5, 0.5)
    )
  }
}

run_dta_nma <- function(input_csv, heterogeneity_prior = "uniform") {
  dta_data <- read_dta_input(input_csv)
  jags_data <- build_jags_data(dta_data)

  params_to_save <- c(
    "sens", "spec", "DOR", "rksens", "rkspec", "rkDOR",
    "tau.sens", "tau.sq.sens", "tau.spec", "tau.sq.spec",
    "tau.logDOR", "tau.sq.logDOR", "rho",
    "sd.study.sens", "sd.study.spec", "sd.rel.sens", "sd.rel.spec",
    "sd.se", "sd.sp"
  )

  model_connection <- textConnection(make_dta_nma_model(heterogeneity_prior))
  on.exit(close(model_connection), add = TRUE)

  set.seed(SEED)
  fit <- jags(
    data = jags_data,
    inits = default_inits(jags_data),
    parameters.to.save = params_to_save,
    model.file = model_connection,
    n.chains = N_CHAINS,
    n.iter = N_ITER,
    n.burnin = N_BURNIN,
    n.thin = N_THIN,
    progress.bar = "none"
  )

  list(fit = fit, data = dta_data, jags_data = jags_data, input_csv = input_csv)
}

summary_to_data_frame <- function(fit) {
  summary_matrix <- fit$BUGSoutput$summary
  data.frame(Parameter = rownames(summary_matrix), summary_matrix, row.names = NULL, check.names = FALSE)
}

interval_text <- function(vals, digits = 2, center = "mean") {
  center_value <- if (center == "median") median(vals) else mean(vals)
  fmt <- paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)")
  sprintf(fmt, center_value, quantile(vals, 0.025), quantile(vals, 0.975))
}

calc_sucra <- function(rank_matrix) {
  nalgo <- ncol(rank_matrix)
  ((nalgo - apply(rank_matrix, 2, mean)) / (nalgo - 1)) * 100
}

format_rank <- function(sucra_vec, idx) {
  ranks <- rank(-sucra_vec, ties.method = "first")
  sprintf("#%d (%.2f%%)", ranks[idx], sucra_vec[idx])
}

make_node_metrics <- function(fit, node_names = NODE_NAMES) {
  sims <- fit$BUGSoutput$sims.list
  nalgo <- length(sims$sens[1, ])
  node_names <- node_names[seq_len(nalgo)]

  sucra_sens <- calc_sucra(sims$rksens)
  sucra_spec <- calc_sucra(sims$rkspec)
  sucra_dor <- calc_sucra(sims$rkDOR)

  rows <- vector("list", nalgo)
  for (i in seq_len(nalgo)) {
    rows[[i]] <- data.frame(
      DiagnosticNode = node_names[i],
      Sensitivity = interval_text(sims$sens[, i]),
      Specificity = interval_text(sims$spec[, i]),
      DOR = interval_text(sims$DOR[, i]),
      SensitivitySUCRA = format_rank(sucra_sens, i),
      SpecificitySUCRA = format_rank(sucra_spec, i),
      DORSUCRA = format_rank(sucra_dor, i),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

make_sensitivity_rows <- function(scenario, metrics) {
  rows <- list()
  for (i in seq_len(nrow(metrics))) {
    rows[[length(rows) + 1]] <- data.frame(
      Scenario = scenario,
      DiagnosticNode = metrics$DiagnosticNode[i],
      Metric = "Sensitivity",
      OverallNetwork = metrics$Sensitivity[i],
      SUCRARank = metrics$SensitivitySUCRA[i],
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1]] <- data.frame(
      Scenario = scenario,
      DiagnosticNode = metrics$DiagnosticNode[i],
      Metric = "Specificity",
      OverallNetwork = metrics$Specificity[i],
      SUCRARank = metrics$SpecificitySUCRA[i],
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1]] <- data.frame(
      Scenario = scenario,
      DiagnosticNode = metrics$DiagnosticNode[i],
      Metric = "DOR",
      OverallNetwork = metrics$DOR[i],
      SUCRARank = metrics$DORSUCRA[i],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

make_league_tables <- function(fit, node_names = NODE_NAMES) {
  sims <- fit$BUGSoutput$sims.list
  nalgo <- length(sims$sens[1, ])
  node_names <- node_names[seq_len(nalgo)]

  league_ss <- matrix("", nrow = nalgo, ncol = nalgo, dimnames = list(node_names, node_names))
  league_dor <- matrix("", nrow = nalgo, ncol = nalgo, dimnames = list(node_names, node_names))

  for (i in seq_len(nalgo)) {
    for (j in seq_len(nalgo)) {
      if (i == j) {
        league_ss[i, j] <- "-"
        league_dor[i, j] <- "-"
      } else if (i < j) {
        league_ss[i, j] <- interval_text(sims$sens[, j] - sims$sens[, i])
        league_dor[i, j] <- interval_text(sims$DOR[, j] / sims$DOR[, i], center = "median")
      } else {
        league_ss[i, j] <- interval_text(sims$spec[, j] - sims$spec[, i])
        league_dor[i, j] <- interval_text(sims$DOR[, j] / sims$DOR[, i], center = "median")
      }
    }
  }

  list(sens_spec = league_ss, rdor = league_dor)
}

make_heterogeneity_summary <- function(fit) {
  summary_matrix <- fit$BUGSoutput$summary
  keep <- intersect(
    c("tau.sens", "tau.sq.sens", "tau.spec", "tau.sq.spec", "tau.logDOR", "tau.sq.logDOR", "rho"),
    rownames(summary_matrix)
  )
  data.frame(Parameter = keep, summary_matrix[keep, , drop = FALSE], row.names = NULL, check.names = FALSE)
}

write_matrix_csv <- function(x, output_path) {
  out <- cbind(Comparison = rownames(x), as.data.frame(x, check.names = FALSE))
  write.csv(out, output_path, row.names = FALSE)
}

write_primary_outputs <- function(result, output_dir) {
  output_dir <- ensure_dir(output_dir)
  fit <- result$fit

  write.csv(summary_to_data_frame(fit), file.path(output_dir, "jags_summary.csv"), row.names = FALSE)
  write.csv(make_heterogeneity_summary(fit), file.path(output_dir, "heterogeneity_threshold_summary.csv"), row.names = FALSE)
  write.csv(make_node_metrics(fit), file.path(output_dir, "node_metrics.csv"), row.names = FALSE)

  league_tables <- make_league_tables(fit)
  write_matrix_csv(league_tables$sens_spec, file.path(output_dir, "league_table_sensitivity_specificity.csv"))
  write_matrix_csv(league_tables$rdor, file.path(output_dir, "league_table_rdor.csv"))

  saveRDS(fit, file.path(output_dir, "jags_fit.rds"))
  invisible(output_dir)
}

make_meta_regression_model <- function(scope) {
  cov_term_sens <- if (scope == "ai_only") "beta.sens * cov[i] * is_ai[i]" else "beta.sens * cov[i]"
  cov_term_spec <- if (scope == "ai_only") "beta.spec * cov[i] * is_ai[i]" else "beta.spec * cov[i]"

  paste0(
    "model {\n",
    "  for(i in 1:nObs){\n",
    "    tp[i] ~ dbin(pi[i,1], pos[i])\n",
    "    tn[i] ~ dbin(pi[i,2], neg[i])\n",
    "    logit(pi[i,1]) <- mu[i,1]\n",
    "    logit(pi[i,2]) <- mu[i,2]\n",
    "    MU[i,1] <- algo.sens[algo[i]] + study.re.sens[s[i]] + (1 - is_baseline_arm[i]) * relative.re.sens[s[i], algo[i]] + ", cov_term_sens, "\n",
    "    MU[i,2] <- algo.spec[algo[i]] + study.re.spec[s[i]] + (1 - is_baseline_arm[i]) * relative.re.spec[s[i], algo[i]] + ", cov_term_spec, "\n",
    "    mu[i,1:2] ~ dmnorm(MU[i,], prec[,])\n",
    "  }\n",
    "  for(j in 1:nalgo){\n",
    "    algo.sens[j] ~ dnorm(0, 0.01)\n",
    "    algo.spec[j] ~ dnorm(0, 0.01)\n",
    "    sens[j] <- exp(algo.sens[j])/(1+exp(algo.sens[j]))\n",
    "    spec[j] <- exp(algo.spec[j])/(1+exp(algo.spec[j]))\n",
    "    DOR[j] <- exp(algo.sens[j] + algo.spec[j])\n",
    "  }\n",
    "  beta.sens ~ dnorm(0, 0.01)\n",
    "  beta.spec ~ dnorm(0, 0.01)\n",
    "  for(k in 1:ns){\n",
    "    study.re.sens[k] ~ dnorm(0, tau.study.sens)\n",
    "    study.re.spec[k] ~ dnorm(0, tau.study.spec)\n",
    "    for(l in 1:nalgo){\n",
    "      relative.re.sens[k,l] ~ dnorm(0, tau.rel.sens)\n",
    "      relative.re.spec[k,l] ~ dnorm(0, tau.rel.spec)\n",
    "    }\n",
    "  }\n",
    "  tau.study.sens <- pow(sd.study.sens, -2)\n",
    "  sd.study.sens ~ dunif(0, 2)\n",
    "  tau.study.spec <- pow(sd.study.spec, -2)\n",
    "  sd.study.spec ~ dunif(0, 2)\n",
    "  tau.rel.sens <- pow(sd.rel.sens, -2)\n",
    "  sd.rel.sens ~ dunif(0, 2)\n",
    "  tau.rel.spec <- pow(sd.rel.spec, -2)\n",
    "  sd.rel.spec ~ dunif(0, 2)\n",
    "  sd.se ~ dunif(0, 2)\n",
    "  sd.sp ~ dunif(0, 2)\n",
    "  rho ~ dunif(-0.99, 0.99)\n",
    "  var.se <- sd.se * sd.se\n",
    "  var.sp <- sd.sp * sd.sp\n",
    "  covar <- rho * sd.se * sd.sp\n",
    "  det <- var.se * var.sp - covar * covar\n",
    "  prec[1,1] <- var.sp / det\n",
    "  prec[2,2] <- var.se / det\n",
    "  prec[1,2] <- -covar / det\n",
    "  prec[2,1] <- prec[1,2]\n",
    "  tau.sens <- sd.se\n",
    "  tau.spec <- sd.sp\n",
    "}\n"
  )
}

run_single_meta_regression <- function(dta_data, covariate, scope) {
  if (!covariate %in% names(dta_data)) {
    return(NULL)
  }

  dta_data$is_ai <- ifelse(dta_data$algo %in% c(3, 4), 1, 0)
  raw_cov <- dta_data[[covariate]]

  if (scope == "ai_only") {
    working_cov <- ifelse(is.na(raw_cov) & dta_data$is_ai == 0, 0, raw_cov)
  } else {
    working_cov <- raw_cov
  }

  valid_idx <- which(!is.na(working_cov))
  if (length(valid_idx) < 10) {
    return(NULL)
  }

  centered_cov <- as.numeric(scale(working_cov[valid_idx], center = TRUE, scale = FALSE))
  if (length(centered_cov) < 2 || is.na(sd(centered_cov)) || sd(centered_cov) < 1e-6) {
    return(NULL)
  }

  jd <- list(
    nalgo = length(unique(dta_data$algo)),
    ns = length(unique(dta_data$s[valid_idx])),
    nObs = length(valid_idx),
    s = as.numeric(as.factor(dta_data$s[valid_idx])),
    algo = dta_data$algo[valid_idx],
    tp = dta_data$tp[valid_idx],
    tn = dta_data$tn[valid_idx],
    pos = dta_data$pos[valid_idx],
    neg = dta_data$neg[valid_idx],
    is_baseline_arm = dta_data$is_baseline_arm[valid_idx],
    cov = centered_cov
  )
  if (scope == "ai_only") {
    jd$is_ai <- dta_data$is_ai[valid_idx]
  }

  model_connection <- textConnection(make_meta_regression_model(scope))
  on.exit(close(model_connection), add = TRUE)

  set.seed(SEED)
  fit <- jags(
    data = jd,
    inits = function() {
      list(
        algo.sens = rep(0, jd$nalgo),
        algo.spec = rep(0, jd$nalgo),
        beta.sens = 0,
        beta.spec = 0,
        sd.study.sens = runif(1, 0.1, 1),
        sd.study.spec = runif(1, 0.1, 1),
        sd.rel.sens = runif(1, 0.1, 1),
        sd.rel.spec = runif(1, 0.1, 1),
        sd.se = runif(1, 0.1, 1),
        sd.sp = runif(1, 0.1, 1),
        rho = runif(1, -0.5, 0.5)
      )
    },
    parameters.to.save = c("beta.sens", "beta.spec", "tau.sens", "tau.spec"),
    model.file = model_connection,
    n.chains = N_CHAINS,
    n.iter = META_N_ITER,
    n.burnin = META_N_BURNIN,
    n.thin = META_N_THIN,
    progress.bar = "none"
  )

  s <- fit$BUGSoutput$summary
  data.frame(
    Covariate = covariate,
    Scope = scope,
    BetaSensMean = round(s["beta.sens", "mean"], 3),
    BetaSensCrI = sprintf("%.3f to %.3f", s["beta.sens", "2.5%"], s["beta.sens", "97.5%"]),
    BetaSensRhat = round(s["beta.sens", "Rhat"], 3),
    BetaSpecMean = round(s["beta.spec", "mean"], 3),
    BetaSpecCrI = sprintf("%.3f to %.3f", s["beta.spec", "2.5%"], s["beta.spec", "97.5%"]),
    BetaSpecRhat = round(s["beta.spec", "Rhat"], 3),
    TauSens = sprintf("%.3f (%.3f to %.3f)", s["tau.sens", "mean"], s["tau.sens", "2.5%"], s["tau.sens", "97.5%"]),
    TauSpec = sprintf("%.3f (%.3f to %.3f)", s["tau.spec", "mean"], s["tau.spec", "2.5%"], s["tau.spec", "97.5%"]),
    NObs = jd$nObs,
    NStudies = jd$ns,
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# Analysis workflow
# -----------------------------

if (RUN_PRIMARY_ANALYSIS) {
  cat("\nRunning primary Bayesian DTA-NMA...\n")
  primary_result <- run_dta_nma(
    input_csv = file.path(DATA_DIR, "input_main_corrected.csv"),
    heterogeneity_prior = "uniform"
  )
  write_primary_outputs(primary_result, file.path(OUTPUT_ROOT, "primary_analysis"))
  cat("Primary analysis completed.\n")
}

if (RUN_META_REGRESSION) {
  cat("\nRunning meta-regression analyses...\n")

  meta_input <- file.path(DATA_DIR, "input_main_corrected.csv")
  meta_data <- read_dta_input(meta_input)

  covariate_registry <- data.frame(
    Covariate = c("Cov_Thickness", "Cov_M2", "Cov_Validation", "Cov_Commercial"),
    Scope = c("global", "global", "ai_only", "ai_only"),
    Label = c(
      "CT Slice Thickness (>=3mm vs <3mm)",
      "M2 Inclusion (Yes vs No)",
      "Validation Type (External vs Internal)",
      "Algorithm Origin (Commercial vs In-house)"
    ),
    stringsAsFactors = FALSE
  )

  meta_results <- list()
  for (i in seq_len(nrow(covariate_registry))) {
    res <- run_single_meta_regression(
      dta_data = meta_data,
      covariate = covariate_registry$Covariate[i],
      scope = covariate_registry$Scope[i]
    )
    if (!is.null(res)) {
      res$Label <- covariate_registry$Label[i]
      meta_results[[length(meta_results) + 1]] <- res
      cat("  Completed:", covariate_registry$Covariate[i], "\n")
    }
  }

  if (length(meta_results) == 0) {
    stop("No meta-regression model was fitted.")
  }

  meta_output_dir <- ensure_dir(file.path(OUTPUT_ROOT, "meta_regression"))
  write.csv(do.call(rbind, meta_results), file.path(meta_output_dir, "meta_regression_summary.csv"), row.names = FALSE)
  cat("Meta-regression analyses completed.\n")
}

if (RUN_SENSITIVITY_ANALYSES) {
  cat("\nRunning scenario-based sensitivity analyses...\n")

  scenario_registry <- data.frame(
    Scenario = c(
      "SA-1 External validation / anterior circulation",
      "SA-2 Commercial AI only",
      "SA-3 M2-including studies only",
      "SA-4 Removing Rai 2025 Subset 1",
      "SA-5 Excluding Sunwoo et al. 2026 US cohort"
    ),
    InputFile = c(
      "NMA_meta_regression_covariates_FINAL-external_corrected.csv",
      "NMA_meta_regression_covariates_FINAL-commercial_corrected.csv",
      "NMA_meta_regression_covariates_FINAL-m2_corrected.csv",
      "NMA_meta_regression_covariates_FINAL-no_subset1_corrected.csv",
      "NMA_meta_regression_covariates_FINAL-no_sunwoo_us_corrected.csv"
    ),
    stringsAsFactors = FALSE
  )

  sensitivity_output_dir <- ensure_dir(file.path(OUTPUT_ROOT, "sensitivity_analyses"))
  scenario_rows <- list()
  sensitivity_table_rows <- list()

  for (i in seq_len(nrow(scenario_registry))) {
    scenario_name <- scenario_registry$Scenario[i]
    input_csv <- file.path(DATA_DIR, scenario_registry$InputFile[i])
    scenario_folder <- sprintf("scenario_%02d", i)

    cat("  Running:", scenario_name, "\n")
    result <- run_dta_nma(input_csv = input_csv, heterogeneity_prior = "uniform")
    write_primary_outputs(result, file.path(sensitivity_output_dir, scenario_folder))

    metrics <- make_node_metrics(result$fit)
    metrics$Scenario <- scenario_name
    metrics$InputFile <- basename(input_csv)
    scenario_rows[[length(scenario_rows) + 1]] <- metrics[
      ,
      c("Scenario", "InputFile", setdiff(names(metrics), c("Scenario", "InputFile")))
    ]
    sensitivity_table_rows[[length(sensitivity_table_rows) + 1]] <- make_sensitivity_rows(scenario_name, metrics)
  }

  write.csv(
    do.call(rbind, scenario_rows),
    file.path(sensitivity_output_dir, "sensitivity_node_metrics.csv"),
    row.names = FALSE
  )
  cat("Sensitivity analyses completed.\n")
}

if (RUN_PRIOR_SENSITIVITY) {
  cat("\nRunning half-normal prior sensitivity analysis...\n")

  prior_result <- run_dta_nma(
    input_csv = file.path(DATA_DIR, "input_main_corrected.csv"),
    heterogeneity_prior = "halfnormal"
  )
  write_primary_outputs(prior_result, file.path(OUTPUT_ROOT, "prior_sensitivity_halfnormal"))
  prior_metrics <- make_node_metrics(prior_result$fit)
  if (!exists("sensitivity_table_rows")) {
    sensitivity_table_rows <- list()
  }
  sensitivity_table_rows[[length(sensitivity_table_rows) + 1]] <- make_sensitivity_rows(
    "SA-6 Half-normal(0,1) heterogeneity prior",
    prior_metrics
  )
  cat("Prior sensitivity analysis completed.\n")
}

if (exists("sensitivity_table_rows") && length(sensitivity_table_rows) > 0) {
  sensitivity_output_dir <- ensure_dir(file.path(OUTPUT_ROOT, "sensitivity_analyses"))
  table14 <- do.call(rbind, sensitivity_table_rows)
  write.csv(
    table14,
    file.path(sensitivity_output_dir, "Supplementary_Table_14_sensitivity_analysis_10study_200k.csv"),
    row.names = FALSE
  )
}

cat("\nAll selected analyses completed.\n")
cat("Results were written to:", OUTPUT_ROOT, "\n")
