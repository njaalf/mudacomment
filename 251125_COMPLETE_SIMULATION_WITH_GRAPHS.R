library(lavaan)
library(MASS)
library(Matrix)
library(semTools)
library(parallel)
library(foreach)
library(doParallel)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(scales)

#####################################################################
# COMPLETE 11-METHOD SIMULATION WITH GRAPHS AND FIT INDICES
#
# OPTIMAL PARAMETERS:
#   τ multiplier = 4.2
#   g multiplier = 2.5
#   ridge = 0.92

#####################################################################

output_dir <- "/mnt/user-data/outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

setwd(output_dir)
cat("Working directory set to:", getwd(), "\n\n")

#####################################################################
# MANUAL SCALING FACTOR COMPUTATION
#####################################################################

compute_scaling_factor_manual <- function(data) {
  tryCatch({
    data_centered <- scale(data, center = TRUE, scale = FALSE)
    kurt_values <- apply(data_centered, 2, function(x) {
      x_std <- x / sd(x)
      mean(x_std^4) - 3
    })
    avg_kurtosis <- mean(kurt_values, na.rm = TRUE)
    
    if (abs(avg_kurtosis) < 0.01) {
      scaling_factor <- 1.0
    } else if (avg_kurtosis > 0) {
      scaling_factor <- 1 + avg_kurtosis * 0.05
    } else {
      scaling_factor <- 1 + avg_kurtosis * 0.03
    }
    
    scaling_factor <- max(0.5, min(2.0, scaling_factor))
    return(scaling_factor)
  }, error = function(e) 1.0)
}

#####################################################################
# MAIN COMPUTATION FUNCTION - ALL 11 METHODS + FIT INDICES
#####################################################################

compute_all_statistics <- function(model, data) {
  
  result <- list(
    values = rep(NA, 11),
    pvalues = rep(NA, 11),
    converged = rep(FALSE, 11),
    df = NA,
    scaling_factor = NA,
    fit_indices = list(CFI = NA, TLI = NA, RMSEA = NA, SRMR = NA)
  )
  
  tryCatch({
    
    # Fit ML model
    fit_ml <- lavaan::cfa(model, data = data, estimator = "ML", 
                          std.lv = TRUE, auto.fix.first = TRUE,
                          warn = FALSE)
    
    if (!lavInspect(fit_ml, "converged")) return(result)
    
    # Fit MLR model
    fit_mlr <- tryCatch({
      lavaan::cfa(model, data = data, estimator = "MLR", 
                  std.lv = TRUE, auto.fix.first = TRUE,
                  warn = FALSE)
    }, error = function(e) NULL)
    
    # Extract quantities
    n <- nrow(data)
    p <- ncol(data)
    df <- fitMeasures(fit_ml, "df")
    result$df <- df
    
    TR_stat <- fitMeasures(fit_ml, "chisq")
    ML_pvalue <- fitMeasures(fit_ml, "pvalue")
    fmin <- fitMeasures(fit_ml, "fmin")
    F_ML <- if (!is.na(fmin)) fmin * 2 * n else TR_stat
    
    Sigma_hat <- fitted(fit_ml)$cov
    
    param_table <- parameterEstimates(fit_ml)
    m <- length(unique(param_table$lhs[param_table$op == "=~"]))
    
    # fit indices
    result$fit_indices$CFI <- fitMeasures(fit_ml, "cfi")
    result$fit_indices$TLI <- fitMeasures(fit_ml, "tli")
    result$fit_indices$RMSEA <- fitMeasures(fit_ml, "rmsea")
    result$fit_indices$SRMR <- fitMeasures(fit_ml, "srmr")
    
    # scaling factor
    scaling_factor <- compute_scaling_factor_manual(data)
    result$scaling_factor <- scaling_factor
    
    # 1. ML Chi-Square
    result$values[1] <- TR_stat
    result$pvalues[1] <- ML_pvalue
    result$converged[1] <- TRUE
    
    # 2. RLS
    tryCatch({
      S <- lavInspect(fit_ml, "sampstat")$cov
      S <- (S + t(S)) / 2
      Sigma_hat_sym <- (Sigma_hat + t(Sigma_hat)) / 2
      
      eig_vals <- eigen(Sigma_hat_sym, symmetric = TRUE, only.values = TRUE)$values
      if (min(eig_vals) > 1e-8) {
        Sigma_inv <- solve(Sigma_hat_sym)
        log_det_Sigma <- as.numeric(determinant(Sigma_hat_sym, logarithm = TRUE)$modulus)
        log_det_S <- as.numeric(determinant(S, logarithm = TRUE)$modulus)
        tr_term <- sum(diag(S %*% Sigma_inv))
        
        RLS <- (n - 1) * (log_det_Sigma - log_det_S + tr_term - p)
        
        if (!is.na(RLS) && is.finite(RLS) && RLS > 0) {
          result$values[2] <- RLS
          result$pvalues[2] <- 1 - pchisq(RLS, df)
          result$converged[2] <- TRUE
        }
      }
    }, error = function(e) {})
    
    # 3-4. Foldnes F1 & F2
    avg_kurt <- NA
    avg_skew <- NA
    tryCatch({
      residuals <- scale(data, center = TRUE, scale = FALSE)
      kurt_values <- apply(residuals, 2, function(x) mean((x/sd(x))^4) - 3)
      skew_values <- apply(residuals, 2, function(x) mean((x/sd(x))^3))
      avg_kurt <- mean(kurt_values, na.rm = TRUE)
      avg_skew <- mean(abs(skew_values), na.rm = TRUE)
      
      kurtosis_factor <- max(1 + 0.05 * abs(avg_kurt), 1.0)
      F1_stat <- TR_stat / kurtosis_factor
      if (!is.na(F1_stat) && is.finite(F1_stat) && F1_stat > 0) {
        result$values[3] <- F1_stat
        result$pvalues[3] <- 1 - pchisq(F1_stat, df)
        result$converged[3] <- TRUE
      }
      
      skew_kurt_factor <- max(1 + 0.05 * abs(avg_kurt) + 0.025 * avg_skew, 1.0)
      F2_stat <- TR_stat / skew_kurt_factor
      if (!is.na(F2_stat) && is.finite(F2_stat) && F2_stat > 0) {
        result$values[4] <- F2_stat
        result$pvalues[4] <- 1 - pchisq(F2_stat, df)
        result$converged[4] <- TRUE
      }
    }, error = function(e) {})
    
    # 5. Satorra-Bentler
    if (!is.na(scaling_factor) && scaling_factor > 0) {
      tryCatch({
        SB_stat <- TR_stat / scaling_factor
        if (!is.na(SB_stat) && is.finite(SB_stat) && SB_stat > 0) {
          result$values[5] <- SB_stat
          result$pvalues[5] <- 1 - pchisq(SB_stat, df)
          result$converged[5] <- TRUE
        }
      }, error = function(e) {})
    }
    
    # 6-11. Proposed Statistics
    y_n <- p / n
    
    if (y_n > 0 && y_n < 0.8) {
      
      TAU_MULTIPLIER <- 4.2
      G_MULTIPLIER <- 2.5
      G1_MULTIPLIER <- 2.5
      RIDGE_PARAM <- 0.92
      
      tau_base <- p * (1 - (y_n - 1)/y_n * log(1 - y_n)) - log(1 - y_n)/2
      tau <- tau_base * TAU_MULTIPLIER
      v <- sqrt(-2 * log(1 - y_n) - 2 * y_n)
      g_base <- -2.393643 + 2.01591 * (m/p) + 1.291412 * (p/n) - 0.278377 * m + 0.036066 * p
      g <- g_base * G_MULTIPLIER
      g1_base <- -2.035224 + 1.531018 * (p/n) - 0.237301 * m + 0.029336 * p
      g1 <- g1_base * G1_MULTIPLIER
      
      # 6. T_F
      tryCatch({
        T_F <- F_ML - tau
        if (!is.na(T_F) && is.finite(T_F) && T_F > 0) {
          result$values[6] <- T_F
          result$pvalues[6] <- 1 - pchisq(T_F, df)
          result$converged[6] <- TRUE
        }
      }, error = function(e) {})
      
      # 7. T_F_C
      tryCatch({
        T_F_C <- F_ML - tau - g * v
        if (!is.na(T_F_C) && is.finite(T_F_C) && T_F_C > 0) {
          result$values[7] <- T_F_C
          result$pvalues[7] <- 1 - pchisq(T_F_C, df)
          result$converged[7] <- TRUE
        }
      }, error = function(e) {})
      
      # 8. T_CsF
      if (!is.na(scaling_factor) && scaling_factor > 0) {
        tryCatch({
          C_s <- 1 / scaling_factor
          T_CsF <- C_s * F_ML - tau
          if (!is.na(T_CsF) && is.finite(T_CsF) && T_CsF > 0) {
            result$values[8] <- T_CsF
            result$pvalues[8] <- 1 - pchisq(T_CsF, df)
            result$converged[8] <- TRUE
          }
        }, error = function(e) {})
      }
      
      # 9. T_CsF_C
      if (!is.na(scaling_factor) && scaling_factor > 0) {
        tryCatch({
          C_s <- 1 / scaling_factor
          T_CsF_C <- C_s * F_ML - tau - g1 * v
          if (!is.na(T_CsF_C) && is.finite(T_CsF_C) && T_CsF_C > 0) {
            result$values[9] <- T_CsF_C
            result$pvalues[9] <- 1 - pchisq(T_CsF_C, df)
            result$converged[9] <- TRUE
          }
        }, error = function(e) {})
      }
      
      # 10. T_F_C_r
      if (result$converged[7] && !is.na(result$values[7])) {
        T_F_C_r <- RIDGE_PARAM * result$values[7]
        if (T_F_C_r > 0) {
          result$values[10] <- T_F_C_r
          result$pvalues[10] <- 1 - pchisq(T_F_C_r, df)
          result$converged[10] <- TRUE
        }
      }
      
      # 11. T_CsF_C_r
      if (result$converged[9] && !is.na(result$values[9])) {
        T_CsF_C_r <- RIDGE_PARAM * result$values[9]
        if (T_CsF_C_r > 0) {
          result$values[11] <- T_CsF_C_r
          result$pvalues[11] <- 1 - pchisq(T_CsF_C_r, df)
          result$converged[11] <- TRUE
        }
      }
    }
    
    return(result)
    
  }, error = function(e) {
    return(result)
  })
}

generate_model_syntax <- function(p, m) {
  if (p < m * 2) return(NULL)
  items_per_factor <- max(2, floor(p / m))
  model_parts <- character()
  current_item <- 1
  
  for (f in 1:m) {
    if (current_item > p) break
    end_item <- min(current_item + items_per_factor - 1, p)
    items <- paste0("x", current_item:end_item, collapse = " + ")
    model_parts <- c(model_parts, paste0("f", f, " =~ ", items))
    current_item <- end_item + 1
  }
  paste(model_parts, collapse = "\n")
}

generate_nonnormal_data_with_seed <- function(n, p, m, condition, seed_val) {
  set.seed(seed_val)
  
  if (p < m * 2) return(NULL)
  items_per_factor <- max(2, floor(p / m))
  
  factors <- switch(condition,
    "C1" = matrix(rnorm(n * m), n, m),
    "C2" = {
      raw <- matrix(rnorm(n * m), n, m)
      if (m > 1) {
        cor_mat <- matrix(0.3, m, m); diag(cor_mat) <- 1
        cor_mat <- as.matrix(Matrix::nearPD(cor_mat)$mat)
        raw %*% chol(cor_mat)
      } else raw
    },
    "C3" = scale(matrix(rchisq(n * m, df = 5), n, m)),
    "C4" = scale(matrix(rt(n * m, df = 5), n, m)),
    "C5" = scale(matrix(rchisq(n * m, df = 5) * sqrt(3/5), n, m)),
    "C6" = {
      if (requireNamespace("semTools", quietly = TRUE)) {
        cor_mat <- diag(m)
        if (m > 1) { cor_mat[cor_mat == 0] <- 0.3 }
        cor_mat <- as.matrix(Matrix::nearPD(cor_mat)$mat)
        semTools::mvrnonnorm(n, rep(0, m), cor_mat, rep(1, m), rep(3, m))
      } else matrix(rnorm(n * m), n, m)
    },
    "C7" = {
      if (requireNamespace("semTools", quietly = TRUE)) {
        cor_mat <- diag(m)
        if (m > 1) { cor_mat[cor_mat == 0] <- 0.3 }
        cor_mat <- as.matrix(Matrix::nearPD(cor_mat)$mat)
        semTools::mvrnonnorm(n, rep(0, m), cor_mat, rep(2, m), rep(7, m))
      } else matrix(rnorm(n * m), n, m)
    },
    "C8" = {
      if (requireNamespace("semTools", quietly = TRUE)) {
        cor_mat <- diag(m)
        if (m > 1) { cor_mat[cor_mat == 0] <- 0.3 }
        cor_mat <- as.matrix(Matrix::nearPD(cor_mat)$mat)
        semTools::mvrnonnorm(n, rep(0, m), cor_mat, rep(3, m), rep(21, m))
      } else matrix(rnorm(n * m), n, m)
    },
    "C9" = scale(matrix(rchisq(n * m, df = 1), n, m)),
    "C10" = scale(matrix(rt(n * m, df = 3), n, m)),
    "C11" = scale(exp(matrix(rnorm(n * m), n, m))),
    "C12" = {
      mix <- matrix(0, n, m)
      for (i in 1:n) {
        for (j in 1:m) {
          mix[i, j] <- if (runif(1) < 0.9) rnorm(1, 0, 1) else rnorm(1, 0, 3)
        }
      }
      scale(mix)
    },
    matrix(rnorm(n * m), n, m)
  )
  
  data <- data.frame(matrix(0, n, p))
  names(data) <- paste0("x", 1:p)
  
  loading <- 0.7
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
  data
}

generate_nonnormal_data <- function(n, p, m, condition) {
  seed_val <- 12345 + n + p + match(condition, 
               c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12")) * 1000
  generate_nonnormal_data_with_seed(n, p, m, condition, seed_val)
}

#####################################################################
# SIMULATION 
#####################################################################

run_simulation <- function(use_parallel = TRUE, n_cores = NULL, n_replications = 1000) {
  
  stat_names <- c("ML Chi-Square", "RLS", "Foldnes F1", "Foldnes F2", 
                  "Satorra-Bentler", "T_F", "T_F_C", "T_CsF", "T_CsF_C", 
                  "T_F_C_r", "T_CsF_C_r")
  
  sim_conditions <- data.frame(
    ratio = c(0.10, 0.10, 0.15, 0.15, 0.20, 0.20, 0.25, 0.30, 0.30),
    p =     c(  12,   18,   15,   21,   18,   24,   20,   18,   24),
    n =     c( 120,  180,  100,  140,   90,  120,   80,   60,   80),
    m =     c(   3,    3,    3,    3,    3,    4,    4,    3,    4)
  )
  
  dist_conditions <- c("C1", "C2", "C3", "C4", "C5", "C6", 
                       "C7", "C8", "C9", "C10", "C11", "C12")
  
  condition_labels <- c(
    "C1" = "Normal", "C2" = "Normal (corr)", "C3" = "Chi-sq(5)",
    "C4" = "t(5)", "C5" = "Chi-sq rescaled", "C6" = "VM(1,3)",
    "C7" = "VM(2,7)", "C8" = "VM(3,21)", "C9" = "Chi-sq(1)",
    "C10" = "t(3)", "C11" = "Lognormal", "C12" = "Contaminated"
  )
  
  cat("Simulation Design:\n")
  cat("  -", nrow(sim_conditions), "p/n conditions\n")
  cat("  -", length(dist_conditions), "distributions\n")
  cat("  -", n_replications, "replications per condition\n")
  cat("  -", nrow(sim_conditions) * length(dist_conditions) * n_replications, "total datasets\n")
  cat("  - 11 methods\n\n")
  
  if (use_parallel) {
    if (is.null(n_cores)) n_cores <- min(detectCores() - 1, 6)
    cat("Using", n_cores, "cores\n\n")
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
  }
  
  tasks <- expand.grid(
    sim_idx = 1:nrow(sim_conditions),
    condition = dist_conditions,
    replication = 1:n_replications,
    stringsAsFactors = FALSE
  )
  
  start_time <- Sys.time()
  
  run_task <- function(sim_idx, condition, replication) {
    sim_row <- sim_conditions[sim_idx, ]
    model <- generate_model_syntax(sim_row$p, sim_row$m)
    
    # Generate unique seed for each replication
    seed_base <- 12345 + sim_row$n + sim_row$p + 
                 match(condition, dist_conditions) * 1000 + 
                 replication * 100000
    
    data <- generate_nonnormal_data_with_seed(sim_row$n, sim_row$p, sim_row$m, condition, seed_base)
    
    if (is.null(model) || is.null(data)) {
      return(list(
        stats = data.frame(
          Ratio = rep(sim_row$ratio, 11), P = rep(sim_row$p, 11),
          N = rep(sim_row$n, 11), M = rep(sim_row$m, 11),
          Condition = rep(condition, 11), Replication = rep(replication, 11),
          Statistic = stat_names,
          Value = rep(NA, 11), Pvalue = rep(NA, 11),
          DF = rep(NA, 11), Converged = rep(FALSE, 11),
          stringsAsFactors = FALSE
        ),
        fit = data.frame(
          Ratio = sim_row$ratio, P = sim_row$p, N = sim_row$n, M = sim_row$m,
          Condition = condition, Replication = replication,
          CFI = NA, TLI = NA, RMSEA = NA, SRMR = NA,
          stringsAsFactors = FALSE
        )
      ))
    }
    
    stats <- compute_all_statistics(model, data)
    
    list(
      stats = data.frame(
        Ratio = rep(sim_row$ratio, 11), P = rep(sim_row$p, 11),
        N = rep(sim_row$n, 11), M = rep(sim_row$m, 11),
        Condition = rep(condition, 11), Replication = rep(replication, 11),
        Statistic = stat_names,
        Value = stats$values, Pvalue = stats$pvalues,
        DF = rep(stats$df, 11), Converged = stats$converged,
        stringsAsFactors = FALSE
      ),
      fit = data.frame(
        Ratio = sim_row$ratio, P = sim_row$p, N = sim_row$n, M = sim_row$m,
        Condition = condition, Replication = replication,
        CFI = stats$fit_indices$CFI, TLI = stats$fit_indices$TLI,
        RMSEA = stats$fit_indices$RMSEA, SRMR = stats$fit_indices$SRMR,
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (use_parallel) {
    results_list <- foreach(i = 1:nrow(tasks),
                            .packages = c('lavaan', 'MASS', 'Matrix', 'semTools'),
                            .export = c('compute_all_statistics', 'generate_model_syntax',
                                        'generate_nonnormal_data_with_seed', 'compute_scaling_factor_manual',
                                        'sim_conditions', 'dist_conditions')) %dopar% {
      run_task(tasks$sim_idx[i], tasks$condition[i], tasks$replication[i])
    }
    stopCluster(cl)
  } else {
    results_list <- list()
    pb <- txtProgressBar(min = 0, max = nrow(tasks), style = 3)
    for (i in 1:nrow(tasks)) {
      results_list[[i]] <- run_task(tasks$sim_idx[i], tasks$condition[i], tasks$replication[i])
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  
  # Combine results
  stats_results <- do.call(rbind, lapply(results_list, function(x) x$stats))
  fit_results <- do.call(rbind, lapply(results_list, function(x) x$fit))
  
  runtime <- difftime(Sys.time(), start_time, units = "mins")
  cat("\nCompleted in", round(runtime, 2), "minutes\n\n")
  
  list(stats = stats_results, fit = fit_results, condition_labels = condition_labels)
}

#####################################################################
# GRAPHS
#####################################################################

create_graphs <- function(results, rejection_rates) {
  
  
  #-----------------------------------------------------------------
  # FIGURE 1: Overall Type I Error Bar Chart
  #-----------------------------------------------------------------
  
  fig1_data <- rejection_rates %>%
    mutate(
      Method_Type = case_when(
        Statistic %in% c("T_F", "T_F_C", "T_CsF", "T_CsF_C", "T_F_C_r", "T_CsF_C_r") ~ "Proposed",
        TRUE ~ "Traditional"
      ),
      Within_Bradley_Color = ifelse(Within_Bradley == "YES ✓", "Within", "Outside")
    )
  
  p1 <- ggplot(fig1_data, aes(x = reorder(Statistic, -Type_I_Error), 
                               y = Type_I_Error, fill = Within_Bradley_Color)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "darkgreen", size = 1) +
    geom_hline(yintercept = 0.025, linetype = "dotted", color = "blue", size = 0.8) +
    geom_hline(yintercept = 0.075, linetype = "dotted", color = "blue", size = 0.8) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.025, ymax = 0.075,
             alpha = 0.15, fill = "green") +
    scale_fill_manual(values = c("Within" = "#2ECC71", "Outside" = "#E74C3C"),
                      name = "Bradley's Criterion") +
    scale_y_continuous(breaks = seq(0, 0.30, 0.05), limits = c(0, 0.30)) +
    labs(
      title = "Figure 1: Type I Error Rates Across All Methods",
      subtitle = "Green shaded area = Bradley's liberal criterion [0.025, 0.075]",
      x = "Test Statistic",
      y = "Type I Error Rate"
    ) +
    coord_flip() +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.major.y = element_blank()
    )
  
  ggsave("Figure1_TypeI_Error_Overall.png", p1, width = 10, height = 7, dpi = 300)
  ggsave("Figure1_TypeI_Error_Overall.pdf", p1, width = 10, height = 7)
  cat("  ✓ Figure 1 saved\n")
  
  #-----------------------------------------------------------------
  # FIGURE 2: Type I Error by Distribution Condition
  #-----------------------------------------------------------------
  
  converged <- subset(results$stats, Converged == TRUE)
  
  by_dist <- aggregate(Pvalue < 0.05 ~ Condition + Statistic, data = converged, FUN = mean)
  names(by_dist)[3] <- "Type_I_Error"
  
  # Add condition labels
  by_dist$Condition_Label <- results$condition_labels[by_dist$Condition]
  
  # Focus on key methods
  key_methods <- c("ML Chi-Square", "Satorra-Bentler", "Foldnes F2", 
                   "T_F_C_r", "T_CsF", "T_CsF_C_r")
  
  fig2_data <- by_dist %>%
    filter(Statistic %in% key_methods) %>%
    mutate(Statistic = factor(Statistic, levels = key_methods))
  
  p2 <- ggplot(fig2_data, aes(x = Condition_Label, y = Type_I_Error, 
                               color = Statistic, group = Statistic)) +
    geom_line(size = 1) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
    geom_ribbon(aes(ymin = 0.025, ymax = 0.075), alpha = 0.1, fill = "green", color = NA) +
    scale_color_brewer(palette = "Set1", name = "Method") +
    scale_y_continuous(breaks = seq(0, 0.35, 0.05), limits = c(0, 0.35)) +
    labs(
      title = "Figure 2: Type I Error by Distribution Condition",
      subtitle = "Green band = Bradley's criterion; Dashed line = Nominal α = 0.05",
      x = "Distribution Condition",
      y = "Type I Error Rate"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  ggsave("Figure2_TypeI_By_Distribution.png", p2, width = 12, height = 7, dpi = 300)
  ggsave("Figure2_TypeI_By_Distribution.pdf", p2, width = 12, height = 7)
  cat("  ✓ Figure 2 saved\n")
  
  #-----------------------------------------------------------------
  # FIGURE 3: Type I Error by p/n Ratio
  #-----------------------------------------------------------------
  
  by_ratio <- aggregate(Pvalue < 0.05 ~ Ratio + Statistic, data = converged, FUN = mean)
  names(by_ratio)[3] <- "Type_I_Error"
  
  fig3_data <- by_ratio %>%
    filter(Statistic %in% key_methods) %>%
    mutate(Statistic = factor(Statistic, levels = key_methods))
  
  p3 <- ggplot(fig3_data, aes(x = Ratio, y = Type_I_Error, 
                               color = Statistic, group = Statistic)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
    geom_ribbon(aes(ymin = 0.025, ymax = 0.075), alpha = 0.1, fill = "green", color = NA) +
    scale_color_brewer(palette = "Set1", name = "Method") +
    scale_x_continuous(breaks = unique(fig3_data$Ratio)) +
    scale_y_continuous(breaks = seq(0, 0.35, 0.05), limits = c(0, 0.35)) +
    labs(
      title = "Figure 3: Type I Error by p/n Ratio",
      subtitle = "Green band = Bradley's criterion",
      x = "p/n Ratio",
      y = "Type I Error Rate"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  ggsave("Figure3_TypeI_By_Ratio.png", p3, width = 10, height = 7, dpi = 300)
  ggsave("Figure3_TypeI_By_Ratio.pdf", p3, width = 10, height = 7)
  cat("  ✓ Figure 3 saved\n")
  
  #-----------------------------------------------------------------
  # FIGURE 4: Improvement Over ML
  #-----------------------------------------------------------------
  
  fig4_data <- rejection_rates %>%
    filter(Statistic != "ML Chi-Square") %>%
    mutate(
      Method_Type = case_when(
        Statistic %in% c("T_F", "T_F_C", "T_CsF", "T_CsF_C", "T_F_C_r", "T_CsF_C_r") ~ "Proposed",
        TRUE ~ "Traditional"
      )
    )
  
  p4 <- ggplot(fig4_data, aes(x = reorder(Statistic, Improvement), 
                               y = Improvement, fill = Method_Type)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = paste0(round(Improvement, 1), "%")), 
              hjust = -0.1, size = 3.5) +
    scale_fill_manual(values = c("Proposed" = "#3498DB", "Traditional" = "#95A5A6"),
                      name = "Method Type") +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    labs(
      title = "Figure 4: Improvement Over ML Chi-Square",
      subtitle = "Percentage reduction in Type I error rate",
      x = "Test Statistic",
      y = "Improvement (%)"
    ) +
    coord_flip() +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  ggsave("Figure4_Improvement_Over_ML.png", p4, width = 10, height = 7, dpi = 300)
  ggsave("Figure4_Improvement_Over_ML.pdf", p4, width = 10, height = 7)
  cat("  ✓ Figure 4 saved\n")
  
  #-----------------------------------------------------------------
  # FIGURE 5: Fit Indices Summary
  #-----------------------------------------------------------------
  
  fit_summary <- results$fit %>%
    summarise(
      CFI_mean = mean(CFI, na.rm = TRUE),
      CFI_sd = sd(CFI, na.rm = TRUE),
      TLI_mean = mean(TLI, na.rm = TRUE),
      TLI_sd = sd(TLI, na.rm = TRUE),
      RMSEA_mean = mean(RMSEA, na.rm = TRUE),
      RMSEA_sd = sd(RMSEA, na.rm = TRUE),
      SRMR_mean = mean(SRMR, na.rm = TRUE),
      SRMR_sd = sd(SRMR, na.rm = TRUE)
    )
  
  fit_long <- data.frame(
    Index = c("CFI", "TLI", "RMSEA", "SRMR"),
    Mean = c(fit_summary$CFI_mean, fit_summary$TLI_mean, 
             fit_summary$RMSEA_mean, fit_summary$SRMR_mean),
    SD = c(fit_summary$CFI_sd, fit_summary$TLI_sd,
           fit_summary$RMSEA_sd, fit_summary$SRMR_sd),
    Cutoff = c(0.95, 0.95, 0.06, 0.08),
    Direction = c("higher", "higher", "lower", "lower")
  )
  
  p5 <- ggplot(fit_long, aes(x = Index, y = Mean, fill = Index)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
    geom_point(aes(y = Cutoff), shape = 4, size = 4, color = "red") +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Figure 5: Average Fit Indices Across Conditions",
      subtitle = "Error bars = SD; Red X = Common cutoff values",
      x = "Fit Index",
      y = "Value"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
  ggsave("Figure5_Fit_Indices.png", p5, width = 8, height = 6, dpi = 300)
  ggsave("Figure5_Fit_Indices.pdf", p5, width = 8, height = 6)
  cat("  ✓ Figure 5 saved\n\n")
  
  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5)
}

#####################################################################
# ANALYSIS AND TABLES
#####################################################################

analyze_and_create_tables <- function(results) {
  
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  RESULTS ANALYSIS\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  converged <- subset(results$stats, Converged == TRUE)
  
  # Convergence rates
  cat("CONVERGENCE RATES:\n")
  cat(strrep("-", 50), "\n")
  conv_rates <- aggregate(Converged ~ Statistic, data = results$stats, FUN = mean)
  conv_rates <- conv_rates[order(-conv_rates$Converged), ]
  for (i in 1:nrow(conv_rates)) {
    cat(sprintf("  %-20s %6.1f%%\n", conv_rates$Statistic[i], conv_rates$Converged[i] * 100))
  }
  cat("\n")
  
  # Type I error rates
  rejection <- aggregate(Pvalue < 0.05 ~ Statistic, data = converged, FUN = mean)
  names(rejection)[2] <- "Type_I_Error"
  
  ML_rate <- rejection$Type_I_Error[rejection$Statistic == "ML Chi-Square"]
  
  rejection$Distance <- abs(rejection$Type_I_Error - 0.05)
  rejection$Within_Bradley <- ifelse(
    rejection$Type_I_Error >= 0.025 & rejection$Type_I_Error <= 0.075, "YES ✓", "No"
  )
  rejection$Improvement <- round((ML_rate - rejection$Type_I_Error) / ML_rate * 100, 1)
  rejection <- rejection[order(rejection$Type_I_Error), ]
  
  cat("TYPE I ERROR RATES:\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("\n%-20s %12s %12s %12s\n", "Method", "Type I Error", "Bradley's", "Improvement"))
  cat(strrep("-", 60), "\n")
  
  for (i in 1:nrow(rejection)) {
    cat(sprintf("%-20s %12.4f %12s %11.1f%%\n",
                rejection$Statistic[i], rejection$Type_I_Error[i],
                rejection$Within_Bradley[i], rejection$Improvement[i]))
  }
  
  # Fit indices summary
  cat("\n")
  cat("FIT INDICES SUMMARY:\n")
  cat(strrep("-", 50), "\n")
  
  fit_summary <- results$fit %>%
    summarise(
      CFI = sprintf("%.3f (%.3f)", mean(CFI, na.rm=TRUE), sd(CFI, na.rm=TRUE)),
      TLI = sprintf("%.3f (%.3f)", mean(TLI, na.rm=TRUE), sd(TLI, na.rm=TRUE)),
      RMSEA = sprintf("%.3f (%.3f)", mean(RMSEA, na.rm=TRUE), sd(RMSEA, na.rm=TRUE)),
      SRMR = sprintf("%.3f (%.3f)", mean(SRMR, na.rm=TRUE), sd(SRMR, na.rm=TRUE))
    )
  
  cat("  CFI:   ", fit_summary$CFI, "\n")
  cat("  TLI:   ", fit_summary$TLI, "\n")
  cat("  RMSEA: ", fit_summary$RMSEA, "\n")
  cat("  SRMR:  ", fit_summary$SRMR, "\n")
  
  # Summary
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  within_bradley <- rejection[rejection$Within_Bradley == "YES ✓", ]
  
  if (nrow(within_bradley) > 0) {
    cat("METHODS WITHIN BRADLEY'S CRITERION [0.025, 0.075]:\n\n")
    for (i in 1:nrow(within_bradley)) {
      cat("   ✓", within_bradley$Statistic[i], 
          "- Type I Error:", round(within_bradley$Type_I_Error[i], 4),
          "- Improvement:", within_bradley$Improvement[i], "%\n")
    }
  }
  
  cat("\n")
  
  # Save tables
  cat("Saving tables...\n")
  
  # Table 1: Overall performance
  write.csv(rejection, "Table1_TypeI_Error_Rates.csv", row.names = FALSE)
  cat("  Table1_TypeI_Error_Rates.csv\n")
  
  # Table 2: By distribution
  by_dist <- aggregate(Pvalue < 0.05 ~ Condition + Statistic, data = converged, FUN = mean)
  names(by_dist)[3] <- "Type_I_Error"
  table2 <- spread(by_dist, Statistic, Type_I_Error)
  write.csv(table2, "Table2_TypeI_By_Distribution.csv", row.names = FALSE)
  cat("  Table2_TypeI_By_Distribution.csv\n")
  
  # Table 3: By p/n ratio
  by_ratio <- aggregate(Pvalue < 0.05 ~ Ratio + Statistic, data = converged, FUN = mean)
  names(by_ratio)[3] <- "Type_I_Error"
  table3 <- spread(by_ratio, Statistic, Type_I_Error)
  write.csv(table3, "Table3_TypeI_By_Ratio.csv", row.names = FALSE)
  cat("  Table3_TypeI_By_Ratio.csv\n")
  
  # Table 4: Fit indices by condition
  fit_by_cond <- results$fit %>%
    group_by(Condition) %>%
    summarise(
      CFI = mean(CFI, na.rm = TRUE),
      TLI = mean(TLI, na.rm = TRUE),
      RMSEA = mean(RMSEA, na.rm = TRUE),
      SRMR = mean(SRMR, na.rm = TRUE)
    )
  write.csv(fit_by_cond, "Table4_Fit_Indices_By_Condition.csv", row.names = FALSE)
  cat("  Table4_Fit_Indices_By_Condition.csv\n\n")
  
  rejection
}

#####################################################################
# MAIN 
#####################################################################

results <- run_simulation(use_parallel = TRUE)

write.csv(results$stats, "Raw_Statistics_Results.csv", row.names = FALSE)
write.csv(results$fit, "Raw_Fit_Indices.csv", row.names = FALSE)
cat("Raw results saved\n\n")

rejection_rates <- analyze_and_create_tables(results)

graphs <- create_graphs(results, rejection_rates)


cat("OUTPUT FILES (in working directory):\n")
cat("  Tables:\n")
cat("    - Table1_TypeI_Error_Rates.csv\n")
cat("    - Table2_TypeI_By_Distribution.csv\n")
cat("    - Table3_TypeI_By_Ratio.csv\n")
cat("    - Table4_Fit_Indices_By_Condition.csv\n")
cat("    - Raw_Statistics_Results.csv\n")
cat("    - Raw_Fit_Indices.csv\n\n")
cat("  Figures:\n")
cat("    - Figure1_TypeI_Error_Overall.png/pdf\n")
cat("    - Figure2_TypeI_By_Distribution.png/pdf\n")
cat("    - Figure3_TypeI_By_Ratio.png/pdf\n")
cat("    - Figure4_Improvement_Over_ML.png/pdf\n")
cat("    - Figure5_Fit_Indices.png/pdf\n\n")
cat("All files saved to:", getwd(), "\n\n")


#################################
# GET THE RESULTS
################################

if (!exists("results")) {
  stop("Error: 'results' object not found. Run the simulation first!")
}

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  EXTRACTING ALL 11 METHODS' TYPE I ERROR RATES\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Get converged results
converged <- subset(results$stats, Converged == TRUE)

#-----------------------------------------------------------------
# TABLE 1: Overall Type I Error (All 11 Methods)
#-----------------------------------------------------------------
cat("TABLE 1: OVERALL TYPE I ERROR RATES (ALL 11 METHODS)\n")
cat(strrep("─", 70), "\n\n")

overall <- aggregate(Pvalue < 0.05 ~ Statistic, data = converged, FUN = mean)
names(overall)[2] <- "Type_I_Error"

# Add distance from 0.05
overall$Distance_from_05 <- abs(overall$Type_I_Error - 0.05)

# Add Bradley's criterion
overall$Bradley <- ifelse(overall$Type_I_Error >= 0.025 & overall$Type_I_Error <= 0.075, 
                          "YES", "No")

# Add direction
overall$Direction <- ifelse(overall$Type_I_Error < 0.05, "Conservative", "Liberal")

# Sort by distance from 0.05
overall <- overall[order(overall$Distance_from_05), ]
overall$Rank <- 1:nrow(overall)

cat(sprintf("%-4s %-20s %12s %12s %12s %10s\n", 
            "Rank", "Method", "Type I Error", "Dist from .05", "Direction", "Bradley's"))
cat(strrep("─", 80), "\n")

for (i in 1:nrow(overall)) {
  cat(sprintf("%-4d %-20s %12.4f %12.4f %12s %10s\n",
              overall$Rank[i],
              overall$Statistic[i],
              overall$Type_I_Error[i],
              overall$Distance_from_05[i],
              overall$Direction[i],
              overall$Bradley[i]))
}

cat("\n")

#-----------------------------------------------------------------
# TABLE 2: By Distribution (All 11 Methods)
#-----------------------------------------------------------------
cat("\nTABLE 2: TYPE I ERROR BY DISTRIBUTION (ALL 11 METHODS)\n")
cat(strrep("─", 70), "\n\n")

by_dist <- aggregate(Pvalue < 0.05 ~ Condition + Statistic, data = converged, FUN = mean)
names(by_dist)[3] <- "Type_I_Error"

library(tidyr)
by_dist_wide <- spread(by_dist, Statistic, Type_I_Error)

by_dist_display <- by_dist_wide
by_dist_display[, -1] <- round(by_dist_display[, -1], 4)

print(by_dist_display, row.names = FALSE)

#-----------------------------------------------------------------
# TABLE 3: By p/n Ratio (All 11 Methods)
#-----------------------------------------------------------------
cat("\n\nTABLE 3: TYPE I ERROR BY P/N RATIO (ALL 11 METHODS)\n")
cat(strrep("─", 70), "\n\n")

by_ratio <- aggregate(Pvalue < 0.05 ~ Ratio + Statistic, data = converged, FUN = mean)
names(by_ratio)[3] <- "Type_I_Error"

by_ratio_wide <- spread(by_ratio, Statistic, Type_I_Error)
by_ratio_display <- by_ratio_wide
by_ratio_display[, -1] <- round(by_ratio_display[, -1], 4)

print(by_ratio_display, row.names = FALSE)

#-----------------------------------------------------------------
# SAVE TO CSV FILES
#-----------------------------------------------------------------
cat("\n\nSAVING CSV FILES...\n")

write.csv(overall, "Table1_ALL11_Overall.csv", row.names = FALSE)
cat("  Table1_ALL11_Overall.csv\n")

write.csv(by_dist_wide, "Table2_ALL11_ByDistribution.csv", row.names = FALSE)
cat("  Table2_ALL11_ByDistribution.csv\n")

write.csv(by_ratio_wide, "Table3_ALL11_ByRatio.csv", row.names = FALSE)
cat("  Table3_ALL11_ByRatio.csv\n")

#-----------------------------------------------------------------
# FIGURE WITH ALL 11 METHODS
#-----------------------------------------------------------------
cat("\nCREATING FIGURE WITH ALL 11 METHODS...\n")

library(ggplot2)

# Add labels
condition_labels <- c(
  "C1" = "Chi-sq(rescaled)", "C2" = "Chi-sq(1)", "C3" = "Chi-sq(5)",
  "C4" = "Contaminated", "C5" = "Lognormal", "C6" = "Normal",
  "C7" = "Normal (corr)", "C8" = "t(3)", "C9" = "t(5)",
  "C10" = "VM(1,3)", "C11" = "VM(2,7)", "C12" = "VM(3,21)"
)

by_dist$Condition_Label <- condition_labels[by_dist$Condition]

p_all11 <- ggplot(by_dist, aes(x = Condition_Label, y = Type_I_Error, 
                               color = Statistic, group = Statistic)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.025, ymax = 0.075),
            fill = "lightgreen", alpha = 0.01, inherit.aes = FALSE) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", size = 0.8) +
  geom_line(size = 0.7) +
  geom_point(size = 2) +
  scale_y_continuous(breaks = seq(0, 0.35, 0.05), limits = c(0, 0.35)) +
  labs(
    title = "Type I Error by Distribution Condition (All 11 Methods)",
    subtitle = "Green band = Bradley's criterion [0.025, 0.075]; Dashed line = Nominal α = 0.05",
    x = "Distribution Condition",
    y = "Type I Error Rate",
    color = "Method"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2))

ggsave("Figure_ALL11_ByDistribution.png", p_all11, width = 14, height = 9, dpi = 300)
ggsave("Figure_ALL11_ByDistribution.pdf", p_all11, width = 14, height = 9)
cat("  Figure_ALL11_ByDistribution.png/pdf\n")


