###########################################################################
# CORRECTED SIMULATION CODE
# Comment on Muda and Yangxuan (2026)
# Foldnes, N. and Grønneberg, S.
#
# Key corrections relative to the original code:
#   1. RLS: uses lavTest(fit_ml, "browne.residual.nt.model")
#      (original computed (n-1)*T_ML/n instead)
#   2. Satorra-Bentler: uses lavTest(fit_ml, "satorra.bentler") with
#      correct scaling factor C_s = r(UG)/df
#      (original used 1 + 0.05*avg_kurtosis as scaling factor)
#   3. "Foldnes F1"/"Foldnes F2" replaced by pEBA4_RLS and pOLS2_RLS
#      from the semTests package (the attributed publication does not exist)
#   4. Factor standardization for C3, C4, C9, C10, C11, C12 uses
#      theoretical moments instead of sample-based scale()
###########################################################################

library(lavaan)
library(MASS)
library(Matrix)
library(semTools)
library(semTests)
library(parallel)
library(foreach)
library(doParallel)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(scales)

#####################################################################
# GLOBAL CONSTANTS
#####################################################################

stat_names <- c("ML Chi-Square", "RLS", "pEBA4_RLS", "pOLS_RLS",
                "Satorra-Bentler", "T_F", "T_F_C", "T_CsF", "T_CsF_C",
                "T_F_C_r", "T_CsF_C_r")

# AR conditions (C1-C4, C6-C12) and non-AR conditions (IG1, IG2, VM1, VM2, VITA)
# Note: C5 excluded (same covariance structure as C4 in original design)
dist_conditions <- c("C1", "C2", "C3", "C4", "C6",
                     "C7", "C8", "C9", "C10", "C11", "C12",
                     "IG1", "IG2", "VM1", "VM2", "VITA")

#####################################################################
# MAIN COMPUTATION FUNCTION - CORRECTED
# 11 methods:
#   1  ML Chi-Square
#   2  RLS              (corrected: lavTest browne.residual.nt.model)
#   3  pEBA4_RLS        (replaces misattributed "Foldnes F1")
#   4  pOLS2_RLS        (replaces misattributed "Foldnes F2")
#   5  Satorra-Bentler  (corrected: lavTest satorra.bentler, C_s = r(UG)/df)
#   6  T_F
#   7  T_F_C
#   8  T_CsF            (corrected: uses proper C_s from lavTest)
#   9  T_CsF_C          (corrected: uses proper C_s from lavTest)
#  10  T_F_C_r
#  11  T_CsF_C_r        (corrected: uses proper C_s from lavTest)
#####################################################################

compute_all_statistics <- function(model, data) {

  result <- list(
    values         = rep(NA, 11),
    pvalues        = rep(NA, 11),
    converged      = rep(FALSE, 11),
    df             = NA,
    scaling_factor = NA,
    fit_indices    = list(CFI = NA, TLI = NA, RMSEA = NA, SRMR = NA)
  )

  fit_ml <- lavaan::cfa(model, data = data, estimator = "ML", warn = FALSE)

  if (!lavInspect(fit_ml, "converged")) return(result)

  # Correct scaling factor and SB statistic via lavTest
  scaled         <- lavTest(fit_ml, "satorra.bentler")
  scaling_factor <- 1 / scaled$scaling.factor  # C_s as defined in paper
  SB_stat        <- scaled$stat

  n  <- nrow(data)
  p  <- ncol(data)
  df <- fitMeasures(fit_ml, "df")
  result$df <- df

  TR_stat   <- fitMeasures(fit_ml, "chisq")
  ML_pvalue <- fitMeasures(fit_ml, "pvalue")
  F_ML      <- TR_stat

  param_table <- parameterEstimates(fit_ml)
  m <- length(unique(param_table$lhs[param_table$op == "=~"]))

  result$fit_indices$CFI   <- fitMeasures(fit_ml, "cfi")
  result$fit_indices$TLI   <- fitMeasures(fit_ml, "tli")
  result$fit_indices$RMSEA <- fitMeasures(fit_ml, "rmsea")
  result$fit_indices$SRMR  <- fitMeasures(fit_ml, "srmr")
  result$scaling_factor    <- scaling_factor

  # 1. ML Chi-Square
  result$values[1]    <- TR_stat
  result$pvalues[1]   <- ML_pvalue
  result$converged[1] <- TRUE

  # 2. RLS (corrected: browne.residual.nt.model)
  tryCatch({
    RLS <- lavTest(fit_ml, test = "browne.residual.nt.model")$stat
    if (!is.na(RLS) && is.finite(RLS) && RLS > 0) {
      result$values[2]    <- RLS
      result$pvalues[2]   <- 1 - pchisq(RLS, df)
      result$converged[2] <- TRUE
    }
  }, error = function(e) {})

  # 3-4. pEBA4_RLS and pOLS2_RLS (replace misattributed "Foldnes F1"/"F2")
  tryCatch({
    tests <- semTests::pvalues(fit_ml, tests = c("pEBA4_RLS", "pOLS2_RLS"))
    result$pvalues[3]   <- tests[1]
    result$converged[3] <- TRUE
    result$pvalues[4]   <- tests[2]
    result$converged[4] <- TRUE
  }, error = function(e) {})

  # 5. Satorra-Bentler (corrected: C_s = r(UG)/df via lavTest)
  if (!is.na(scaling_factor) && scaling_factor > 0) {
    tryCatch({
      if (!is.na(SB_stat) && is.finite(SB_stat) && SB_stat > 0) {
        result$values[5]    <- SB_stat
        result$pvalues[5]   <- 1 - pchisq(SB_stat, df)
        result$converged[5] <- TRUE
      }
    }, error = function(e) {})
  }

  # 6-11. Proposed statistics (T_F family)
  y_n <- p / n

  if (y_n > 0 && y_n < 0.8) {

    TAU_MULTIPLIER <- 4.2
    G_MULTIPLIER   <- 2.5
    G1_MULTIPLIER  <- 2.5
    RIDGE_PARAM    <- 0.92

    tau_base <- p * (1 - (y_n - 1)/y_n * log(1 - y_n)) - log(1 - y_n)/2
    tau <- tau_base * TAU_MULTIPLIER
    v   <- sqrt(-2 * log(1 - y_n) - 2 * y_n)
    g_base  <- -2.393643 + 2.01591*(m/p) + 1.291412*(p/n) - 0.278377*m + 0.036066*p
    g       <- g_base * G_MULTIPLIER
    g1_base <- -2.035224 + 1.531018*(p/n) - 0.237301*m + 0.029336*p
    g1      <- g1_base * G1_MULTIPLIER

    # 6. T_F
    tryCatch({
      T_F <- max(F_ML - tau, 0)
      if (!is.na(T_F) && is.finite(T_F)) {
        result$values[6]    <- T_F
        result$pvalues[6]   <- 1 - pchisq(T_F, df)
        result$converged[6] <- TRUE
      }
    }, error = function(e) {})

    # 7. T_F_C
    tryCatch({
      T_F_C <- max(F_ML - tau - g * v, 0)
      if (!is.na(T_F_C) && is.finite(T_F_C)) {
        result$values[7]    <- T_F_C
        result$pvalues[7]   <- 1 - pchisq(T_F_C, df)
        result$converged[7] <- TRUE
      }
    }, error = function(e) {})

    # 8. T_CsF (corrected: uses C_s = scaling_factor from lavTest)
    if (!is.na(scaling_factor) && scaling_factor > 0) {
      tryCatch({
        T_CsF <- max(scaling_factor * F_ML - tau, 0)
        if (!is.na(T_CsF) && is.finite(T_CsF)) {
          result$values[8]    <- T_CsF
          result$pvalues[8]   <- 1 - pchisq(T_CsF, df)
          result$converged[8] <- TRUE
        }
      }, error = function(e) {})
    }

    # 9. T_CsF_C (corrected: uses C_s = scaling_factor from lavTest)
    if (!is.na(scaling_factor) && scaling_factor > 0) {
      tryCatch({
        T_CsF_C <- max(scaling_factor * F_ML - tau - g1 * v, 0)
        if (!is.na(T_CsF_C) && is.finite(T_CsF_C)) {
          result$values[9]    <- T_CsF_C
          result$pvalues[9]   <- 1 - pchisq(T_CsF_C, df)
          result$converged[9] <- TRUE
        }
      }, error = function(e) {})
    }

    # 10. T_F_C_r
    if (result$converged[7] && !is.na(result$values[7])) {
      T_F_C_r <- RIDGE_PARAM * result$values[7]
      if (T_F_C_r > 0) {
        result$values[10]    <- T_F_C_r
        result$pvalues[10]   <- 1 - pchisq(T_F_C_r, df)
        result$converged[10] <- TRUE
      }
    }

    # 11. T_CsF_C_r
    if (result$converged[9] && !is.na(result$values[9])) {
      T_CsF_C_r <- RIDGE_PARAM * result$values[9]
      if (T_CsF_C_r > 0) {
        result$values[11]    <- T_CsF_C_r
        result$pvalues[11]   <- 1 - pchisq(T_CsF_C_r, df)
        result$converged[11] <- TRUE
      }
    }
  }

  return(result)
}

#####################################################################
# MODEL SYNTAX GENERATOR
#####################################################################

generate_model_syntax <- function(p, m) {
  if (p < m * 2) return(NULL)
  items_per_factor <- max(2, floor(p / m))
  model_parts  <- character()
  current_item <- 1

  for (f in 1:m) {
    if (current_item > p) break
    end_item <- min(current_item + items_per_factor - 1, p)
    items <- paste0("x", current_item:end_item, collapse = " + ")
    model_parts  <- c(model_parts, paste0("f", f, " =~ ", items))
    current_item <- end_item + 1
  }
  paste(model_parts, collapse = "\n")
}

#####################################################################
# DATA GENERATION
#
# AR conditions (C1-C4, C6-C12): factors generated, then manifest
#   variables = 0.7 * factor + error.
#   Corrections for theoretical standardization (replaces scale()):
#     C3:  (chi-sq(5) - 5) / sqrt(10)
#     C4:  t(5) / sqrt(5/3)
#     C9:  (chi-sq(1) - 1) / sqrt(2)
#     C10: t(3) / sqrt(3)
#     C11: (exp(Z) - exp(0.5)) / sqrt(exp(2) - exp(1))
#     C12: mix / sqrt(1.8)
#
# Non-AR conditions (IG1, IG2, VM1, VM2, VITA): manifest variables
#   generated directly from the model-implied covariance matrix.
#####################################################################

generate_nonnormal_data_with_seed <- function(n, p, m, condition, seed_val) {
  set.seed(seed_val)

  if (p < m * 2) return(NULL)
  items_per_factor <- max(2, floor(p / m))

  factors <- switch(condition,
    "C1"  = matrix(rnorm(n * m), n, m),
    "C2"  = {
      raw <- matrix(rnorm(n * m), n, m)
      if (m > 1) {
        cor_mat <- matrix(0.3, m, m); diag(cor_mat) <- 1
        cor_mat <- as.matrix(Matrix::nearPD(cor_mat)$mat)
        raw %*% chol(cor_mat)
      } else raw
    },
    "C3"  = (matrix(rchisq(n * m, df = 5), n, m) - 5) / sqrt(10),
    "C4"  = matrix(rt(n * m, df = 5), n, m) / sqrt(5/3),
    "C5"  = NULL,  # C5: non-AR, manifest variables generated directly
    "C6"  = {
      cor_mat <- diag(m)
      if (m > 1) cor_mat[cor_mat == 0] <- 0.3
      cor_mat <- as.matrix(Matrix::nearPD(cor_mat)$mat)
      semTools::mvrnonnorm(n, rep(0, m), cor_mat, rep(1, m), rep(3, m))
    },
    "C7"  = {
      cor_mat <- diag(m)
      if (m > 1) cor_mat[cor_mat == 0] <- 0.3
      cor_mat <- as.matrix(Matrix::nearPD(cor_mat)$mat)
      semTools::mvrnonnorm(n, rep(0, m), cor_mat, rep(2, m), rep(7, m))
    },
    "C8"  = {
      cor_mat <- diag(m)
      if (m > 1) cor_mat[cor_mat == 0] <- 0.3
      cor_mat <- as.matrix(Matrix::nearPD(cor_mat)$mat)
      semTools::mvrnonnorm(n, rep(0, m), cor_mat, rep(3, m), rep(21, m))
    },
    "C9"  = (matrix(rchisq(n * m, df = 1), n, m) - 1) / sqrt(2),
    "C10" = matrix(rt(n * m, df = 3), n, m) / sqrt(3),
    "C11" = (exp(matrix(rnorm(n * m), n, m)) - exp(0.5)) / sqrt(exp(2) - exp(1)),
    "C12" = {
      mix <- matrix(0, n, m)
      for (i in 1:n) {
        for (j in 1:m) {
          mix[i, j] <- if (runif(1) < 0.9) rnorm(1, 0, 1) else rnorm(1, 0, 3)
        }
      }
      mix / sqrt(1.8)
    },
    matrix(rnorm(n * m), n, m)
  )

  if (condition %in% setdiff(paste0("C", 1:12), "C5")) {  # AR conditions
    data <- data.frame(matrix(0, n, p))
    names(data) <- paste0("x", 1:p)
    loading   <- 0.7
    error_var <- 1 - loading^2
    current_item <- 1
    for (f in 1:m) {
      if (current_item > p) break
      end_item <- min(current_item + items_per_factor - 1, p)
      for (i in current_item:end_item) {
        if (i > p) break
        data[, i] <- loading * factors[, f] + rnorm(n, 0, sqrt(error_var))
      }
      current_item <- end_item + 1
    }
  } else {  # Non-AR: generate on manifest variables directly
    Lambda <- matrix(0, nrow = p, ncol = m)
    current_item <- 1
    for (f in 1:m) {
      end_item <- if (f == m) p else (current_item + items_per_factor - 1)
      Lambda[current_item:end_item, f] <- 0.7
      current_item <- end_item + 1
    }
    Phi          <- matrix(0.3, nrow = m, ncol = m); diag(Phi) <- 1
    Psi          <- diag(1 - 0.7^2, p)
    sigma_target <- Lambda %*% Phi %*% t(Lambda) + Psi

    if (condition == "IG1") {
      data_matrix <- covsim::rIG(N = n, sigma.target = sigma_target,
                                 skewness = rep(2, p), excesskurtosis = rep(7, p))[[1]]
    } else if (condition == "IG2") {
      data_matrix <- covsim::rIG(N = n, sigma.target = sigma_target,
                                 skewness = rep(3, p), excesskurtosis = rep(21, p))[[1]]
    } else if (condition == "VM1") {
      data_matrix <- semTools::mvrnonnorm(n, mu = rep(0, p), Sigma = sigma_target,
                                         rep(2, p), rep(7, p))
    } else if (condition == "VM2") {
      data_matrix <- semTools::mvrnonnorm(n, mu = rep(0, p), Sigma = sigma_target,
                                         rep(3, p), rep(21, p))
    } else if (condition == "VITA") {
      tmp <- data.frame(
        ratio = c(0.10, 0.10, 0.15, 0.15, 0.20, 0.20, 0.25, 0.30, 0.30),
        p     = c(  12,   18,   15,   21,   18,   24,   20,   18,   24),
        n     = c( 120,  180,  100,  140,   90,  120,   80,   60,   80),
        m     = c(   3,    3,    3,    3,    3,    4,    4,    3,    4)
      )
      sim_idx_vita <- which(tmp$m == m & tmp$p == p & tmp$n == n)
      cvita        <- calibrated[[sim_idx_vita]][[1]]
      post         <- calibrated[[sim_idx_vita]][[2]]
      data_matrix  <- rvinecopulib::rvine(n, cvita) %*% post
      colnames(data_matrix) <- paste0("x", 1:p)
    }
    data <- as.data.frame(data_matrix)
    names(data) <- paste0("x", 1:p)
  }
  return(data)
}

generate_nonnormal_data <- function(n, p, m, condition) {
  seed_val <- 12345 + n + p +
    match(condition, c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12",
                       "C13","C14","C15")) * 1000
  generate_nonnormal_data_with_seed(n, p, m, condition, seed_val)
}

#####################################################################
# SIMULATION RUNNER - AR AND NON-AR CONDITIONS
#####################################################################

run_simulation <- function(use_parallel = TRUE, dist_conditions,
                           n_cores = NULL, n_replications = 5) {
  start_time <- Sys.time()

  sim_conditions <- data.frame(
    ratio = c(0.10, 0.10, 0.15, 0.15, 0.20, 0.20, 0.25, 0.30, 0.30),
    p     = c(  12,   18,   15,   21,   18,   24,   20,   18,   24),
    n     = c( 120,  180,  100,  140,   90,  120,   80,   60,   80),
    m     = c(   3,    3,    3,    3,    3,    4,    4,    3,    4)
  )

  tasks <- expand.grid(
    sim_idx     = 1:nrow(sim_conditions),
    condition   = dist_conditions,
    replication = 1:n_replications,
    stringsAsFactors = FALSE
  )

  run_task_internal <- function(sim_idx, condition, replication) {
    tryCatch({
      sim_row   <- sim_conditions[sim_idx, ]
      model     <- generate_model_syntax(sim_row$p, sim_row$m)
      seed_base <- 12345 + sim_row$n + sim_row$p +
        match(condition, dist_conditions) * 1000 + replication * 100000
      data  <- generate_nonnormal_data_with_seed(
                 sim_row$n, sim_row$p, sim_row$m, condition, seed_base)
      if (is.null(data)) return(NULL)
      stats <- compute_all_statistics(model, data)
      if (is.null(stats) || is.na(stats$df)) return(NULL)
      data.frame(
        Ratio       = rep(sim_row$ratio, 11),
        P           = rep(sim_row$p,     11),
        N           = rep(sim_row$n,     11),
        M           = rep(sim_row$m,     11),
        Condition   = rep(condition,     11),
        Replication = rep(replication,   11),
        Statistic   = stat_names,
        Value       = stats$values,
        Pvalue      = stats$pvalues,
        DF          = rep(stats$df,      11),
        Converged   = stats$converged,
        stringsAsFactors = FALSE
      )
    }, error = function(e) paste("Error:", e$message))
  }

  if (use_parallel) {
    if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
    cat("Firing up", n_cores, "cores...\n")
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl, varlist = ls(envir = environment()),
                            envir = environment())
    parallel::clusterExport(cl,
      varlist = c("compute_all_statistics", "generate_model_syntax",
                  "generate_nonnormal_data_with_seed", "calibrated",
                  "stat_names", "dist_conditions"))

    results_list <- foreach(
      i              = 1:nrow(tasks),
      .packages      = c("lavaan", "MASS", "Matrix", "semTests",
                         "semTools", "rvinecopulib"),
      .errorhandling = "pass"
    ) %dopar% {
      run_task_internal(tasks$sim_idx[i], tasks$condition[i], tasks$replication[i])
    }
    parallel::stopCluster(cl)
  } else {
    results_list <- lapply(1:nrow(tasks), function(i)
      run_task_internal(tasks$sim_idx[i], tasks$condition[i], tasks$replication[i]))
  }

  is_df <- sapply(results_list, is.data.frame)
  if (sum(is_df) == 0) {
    cat("CRITICAL ERROR: All tasks failed. First error message:\n")
    print(results_list[[1]])
    return(NULL)
  }

  stats_results <- do.call(rbind, results_list[is_df])
  runtime <- difftime(Sys.time(), start_time, units = "mins")
  cat("\nCompleted in", round(runtime, 2), "minutes. Success:",
      sum(is_df), "/", length(results_list), "\n")
  return(stats_results)
}

#####################################################################
# SIMULATION RUNNER - PIECEWISE LINEAR (PL) CONDITIONS
#####################################################################

run_simulation_plsim <- function(use_parallel = TRUE, n_cores = 18,
                                  n_replications = 5000) {
  start_time <- Sys.time()

  sim_conditions <- data.frame(
    ratio = c(0.10, 0.10, 0.15, 0.15, 0.20, 0.20, 0.25, 0.30, 0.30),
    p     = c(  12,   18,   15,   21,   18,   24,   20,   18,   24),
    n     = c( 120,  180,  100,  140,   90,  120,   80,   60,   80),
    m     = c(   3,    3,    3,    3,    3,    4,    4,    3,    4)
  )

  plsim_levels <- list(
    "PL1" = list(skew = 2, kurt = 7,  numsegments = 4),
    "PL2" = list(skew = 3, kurt = 21, numsegments = 10)
  )

  outer_tasks <- expand.grid(
    sim_idx = 1:nrow(sim_conditions),
    pl_name = names(plsim_levels),
    stringsAsFactors = FALSE
  )

  run_outer_task <- function(sim_idx, pl_name) {
    sim_row <- sim_conditions[sim_idx, ]
    pl      <- plsim_levels[[pl_name]]
    model   <- generate_model_syntax(sim_row$p, sim_row$m)

    items_per_factor <- max(2, floor(sim_row$p / sim_row$m))
    Lambda       <- matrix(0, nrow = sim_row$p, ncol = sim_row$m)
    current_item <- 1
    for (f in 1:sim_row$m) {
      end_item <- if (f == sim_row$m) sim_row$p else (current_item + items_per_factor - 1)
      Lambda[current_item:end_item, f] <- 0.7
      current_item <- end_item + 1
    }
    Phi          <- matrix(0.3, nrow = sim_row$m, ncol = sim_row$m); diag(Phi) <- 1
    Psi          <- diag(1 - 0.7^2, sim_row$p)
    sigma_target <- Lambda %*% Phi %*% t(Lambda) + Psi

    set.seed(12345 + sim_idx + (pl_name == "PL2") * 500)
    samples <- covsim::rPLSIM(
      N              = sim_row$n,
      reps           = n_replications,
      sigma.target   = sigma_target,
      skewness       = rep(pl$skew, sim_row$p),
      excesskurtosis = rep(pl$kurt, sim_row$p),
      numsegments    = pl$numsegments
    )[[1]]

    results <- lapply(seq_len(n_replications), function(rep) {
      tryCatch({
        data  <- as.data.frame(samples[[rep]])
        names(data) <- paste0("x", 1:sim_row$p)
        stats <- compute_all_statistics(model, data)
        if (is.null(stats) || is.na(stats$df)) return(NULL)
        data.frame(
          Ratio       = rep(sim_row$ratio, 11),
          P           = rep(sim_row$p,     11),
          N           = rep(sim_row$n,     11),
          M           = rep(sim_row$m,     11),
          Condition   = rep(pl_name,       11),
          Replication = rep(rep,           11),
          Statistic   = stat_names,
          Value       = stats$values,
          Pvalue      = stats$pvalues,
          DF          = rep(stats$df,      11),
          Converged   = stats$converged,
          stringsAsFactors = FALSE
        )
      }, error = function(e) NULL)
    })

    is_df <- sapply(results, is.data.frame)
    do.call(rbind, results[is_df])
  }

  if (use_parallel) {
    cat("Firing up", n_cores, "cores for", nrow(outer_tasks), "outer tasks...\n")
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl,
      varlist = c("compute_all_statistics", "generate_model_syntax",
                  "sim_conditions", "plsim_levels", "stat_names", "n_replications"),
      envir   = environment())

    results_list <- foreach(
      i              = 1:nrow(outer_tasks),
      .packages      = c("lavaan", "MASS", "Matrix", "semTests", "semTools", "covsim"),
      .errorhandling = "pass"
    ) %dopar% {
      run_outer_task(outer_tasks$sim_idx[i], outer_tasks$pl_name[i])
    }
    parallel::stopCluster(cl)
  } else {
    results_list <- lapply(1:nrow(outer_tasks), function(i)
      run_outer_task(outer_tasks$sim_idx[i], outer_tasks$pl_name[i]))
  }

  is_df <- sapply(results_list, is.data.frame)
  if (sum(is_df) == 0) {
    cat("CRITICAL ERROR: All outer tasks failed.\n")
    print(results_list[[1]])
    return(NULL)
  }

  stats_results <- do.call(rbind, results_list[is_df])
  runtime <- difftime(Sys.time(), start_time, units = "mins")
  cat("\nCompleted in", round(runtime, 2), "minutes. Success:",
      sum(is_df), "/", nrow(outer_tasks), "outer tasks\n")
  stats_results
}

#####################################################################
# VITA EXTENSION: Calibrate vine copula distributions
# Run this calibration block once; it saves calibratedvita.rds.
# On subsequent runs you can reload with:
#   calibrated <- readRDS("calibratedvita.rds")
#####################################################################

sim_conditions <- data.frame(
  ratio = c(0.10, 0.10, 0.15, 0.15, 0.20, 0.20, 0.25, 0.30, 0.30),
  p     = c(  12,   18,   15,   21,   18,   24,   20,   18,   24),
  n     = c( 120,  180,  100,  140,   90,  120,   80,   60,   80),
  m     = c(   3,    3,    3,    3,    3,    4,    4,    3,    4)
)
target_df <- 1.2
loading   <- 0.7

# calibrated <- apply(sim_conditions[, c("p", "m")], 1, function(x) {
#   p <- x[1]; m <- x[2]
#   Lambda       <- matrix(0, nrow = p, ncol = m)
#   items_per_factor <- p %/% m
#   current_item <- 1
#   for (f in 1:m) {
#     end_item <- if (f == m) p else (current_item + items_per_factor - 1)
#     Lambda[current_item:end_item, f] <- loading
#     current_item <- end_item + 1
#   }
#   Phi          <- matrix(0.3, nrow = m, ncol = m); diag(Phi) <- 1
#   Psi          <- diag(1 - loading^2, p)
#   sigma.target <- Lambda %*% Phi %*% t(Lambda) + Psi
# 
#   margins          <- lapply(1:p, function(x) list(distr = "chisq", df = target_df))
#   margin.variances <- rep(target_df * 2, p)
#   pre              <- diag(sqrt(margin.variances / diag(sigma.target)))
#   vita.target      <- pre %*% sigma.target %*% pre
#   set.seed(1)
#   calibrated.vita  <- vita(margins, vita.target, verbose = TRUE)
#   post             <- diag(1 / diag(pre))
#   list(v = calibrated.vita, post = post)
# })
# 
# saveRDS(calibrated, "calibratedvita.rds")

#####################################################################
# VITA SIMULATION RUNNER
#####################################################################

run_simulation_vita <- function(use_parallel = FALSE, n_cores = NULL,
                                 n_replications = 2000) {
  sim_conditions <- data.frame(
    ratio = c(0.10, 0.10, 0.15, 0.15, 0.20, 0.20, 0.25, 0.30, 0.30),
    p     = c(  12,   18,   15,   21,   18,   24,   20,   18,   24),
    n     = c( 120,  180,  100,  140,   90,  120,   80,   60,   80),
    m     = c(   3,    3,    3,    3,    3,    4,    4,    3,    4)
  )

  tasks <- expand.grid(
    sim_idx     = 1:nrow(sim_conditions),
    replication = 1:n_replications,
    stringsAsFactors = FALSE
  )

  run_task <- function(sim_idx, replication) {
    sim_row <- sim_conditions[sim_idx, ]
    model   <- generate_model_syntax(sim_row$p, sim_row$m)
    cvita   <- calibrated[[sim_idx]][[1]]
    post    <- calibrated[[sim_idx]][[2]]
    data    <- rvinecopulib::rvine(sim_row$n, cvita) %*% post %>% as.data.frame()
    colnames(data) <- paste0("x", 1:sim_row$p)
    stats   <- compute_all_statistics(model, data)
    data.frame(pvals = stats$pvalues, replication = replication, sim_row = sim_row)
  }

  results_list <- list()
  pb <- txtProgressBar(min = 0, max = nrow(tasks), style = 3)
  for (i in 1:nrow(tasks)) {
    results_list[[i]] <- run_task(tasks$sim_idx[i], tasks$replication[i])
    setTxtProgressBar(pb, i)
  }
  close(pb)
  results_list
}

#####################################################################
# MAIN RUNS
#####################################################################

results <- run_simulation(use_parallel = TRUE, dist_conditions = dist_conditions,
                          n_replications = 5)

closeAllConnections()
if (exists("cl")) try(stopCluster(cl), silent = TRUE)
gc()

# PL conditions (calibration run internally per call)
plresults <- run_simulation_plsim(n_replications = 5)

pl_rr <- plresults %>%
  group_by(Condition, Statistic) %>%
  summarise(rr = mean(Pvalue < 0.05), .groups = "drop") %>%
  pivot_wider(names_from = Statistic, values_from = rr) %>%
  select(Condition, all_of(stat_names)) %>%
  mutate(Condition = recode(Condition,
                            "PL1" = "PL s2,k7",
                            "PL2" = "PL s3,k21"))

# VITA run (requires calibrated object from calibration block above)
res <- run_simulation_vita(use_parallel = FALSE, n_replications = 5)

res2 <- lapply(res, function(x) { x$statistic <- stat_names; x })
df   <- do.call(rbind, res2)
df %>% group_by(statistic) %>% summarise(n(), rr = mean(pvals < .05))

