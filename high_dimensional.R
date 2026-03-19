library(dplyr)
library(knitr)
library(tidyr)

res.df <- readRDS("highdim.rds")
unique(res.df$dim/res.df$n) %>% sort
# Check available conditions
filter(res.df, dim/n == .125) %>% distinct(dim, n)
filter(res.df, dim/n == .15)  %>% distinct(dim, n)
filter(res.df, dim/n == .25)  %>% distinct(dim, n)

# Check dist values in the data
unique(res.df$dist)


#   T_CsF, T_CsF_C, T_CsF_C_r,   <- Cs/SB-based ridge (stats 8,9,11)
#   T_F, T_F_C, T_F_C_r           <- ML-based ridge    (stats 6,7,10)
tests <- c("std_ml", "std_rls", "sb", "peba4_rls", "pols2_rls",
           "T_CsF", "T_CsF_C", "T_CsF_C_r",
           "T_F", "T_F_C", "T_F_C_r")

# Muda constants
TAU_MULTIPLIER <- 4.2
G_MULTIPLIER   <- 2.5
G1_MULTIPLIER  <- 2.5
RIDGE_PARAM    <- 0.92
# Cap chi-sq when p-value = 0 (qchisq(1, df) = Inf)
LARGE_STAT     <- 1e5

compute_ridge_stats <- function(data, p_val, n_val, df_val, m_val = 5) {

  y_n <- p_val / n_val

  tau_base <- p_val * (1 - (y_n - 1)/y_n * log(1 - y_n)) - log(1 - y_n)/2
  tau <- tau_base * TAU_MULTIPLIER
  v   <- sqrt(-2 * log(1 - y_n) - 2 * y_n)

  # g for ML-based corrected (T_F_C); g1 for Cs-based corrected (T_CsF_C)
  g_base  <- -2.393643 + 2.01591*(m_val/p_val) + 1.291412*(p_val/n_val) -
              0.278377*m_val + 0.036066*p_val
  g       <- g_base * G_MULTIPLIER

  g1_base <- -2.035224 + 1.531018*(p_val/n_val) - 0.237301*m_val + 0.029336*p_val
  g1      <- g1_base * G1_MULTIPLIER

  data %>%
    filter(dim == p_val & n == n_val) %>%
    mutate(
      # Recover chi-sq values from p-values; cap at LARGE_STAT when p-value = 0
      T_ML_stat = pmin(qchisq(1 - std_ml, df = df_val), LARGE_STAT),
      T_SB_stat = pmin(qchisq(1 - sb,     df = df_val), LARGE_STAT),

      # Intermediate chi-sq values (floored at 0 before ridge multiplication)
      # Cs/SB-based (stats 8, 9, 11)
      T_CsF_stat     = pmax(T_SB_stat - tau,           0),
      T_CsF_C_stat   = pmax(T_SB_stat - tau - g1 * v,  0),
      T_CsF_C_r_stat = RIDGE_PARAM * T_CsF_C_stat,
      # ML-based (stats 6, 7, 10)
      T_F_stat       = pmax(T_ML_stat - tau,           0),
      T_F_C_stat     = pmax(T_ML_stat - tau - g  * v,  0),
      T_F_C_r_stat   = RIDGE_PARAM * T_F_C_stat,

      # Convert to p-values
      T_CsF     = 1 - pchisq(T_CsF_stat,     df_val),
      T_CsF_C   = 1 - pchisq(T_CsF_C_stat,   df_val),
      T_CsF_C_r = 1 - pchisq(T_CsF_C_r_stat, df_val),
      T_F       = 1 - pchisq(T_F_stat,        df_val),
      T_F_C     = 1 - pchisq(T_F_C_stat,      df_val),
      T_F_C_r   = 1 - pchisq(T_F_C_r_stat,    df_val)
    ) %>%
    group_by(dist) %>%
    summarise(across(all_of(tests), ~ mean(.x < 0.05, na.rm = TRUE))) %>%
    mutate(across(where(is.numeric), ~ round(., 3)))
}

# Setting 1: p = 100, n = 800 (p/n = 0.125)
res_125 <- compute_ridge_stats(res.df, p_val = 100, n_val = 800, df_val = 4840)
kable(res_125, caption = "Rejection Rates by Dist (p=100, n=800, p/n=0.125)")

# Setting 2: p = 60, n = 400 (p/n = 0.15)
res_150 <- compute_ridge_stats(res.df, p_val = 60, n_val = 400, df_val = 1700)
kable(res_150, caption = "Rejection Rates by Dist (p=60, n=400, p/n=0.15)")

# Setting 3: p = 100, n = 400 (p/n = 0.25)
res_250 <- compute_ridge_stats(res.df, p_val = 100, n_val = 400, df_val = 4840)
kable(res_250, caption = "Rejection Rates by Dist (p=100, n=400, p/n=0.25)")

