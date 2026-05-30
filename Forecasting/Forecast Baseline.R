###############################################################
# Vasicek + CIR FORECATS
# 5-year rolling window, h = 6 months, n_starts_oos = 10
###############################################################


suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(ggplot2)
})


SCRIPT_DIR <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
setwd(SCRIPT_DIR)

set.seed(123)

swap_file <- "~/Desktop/ægte data/EURIBOR swap data.xlsx"
start_date <- as.Date("2008-10-06")

taus <- c(1, 2, 3, 5, 7, 10, 15, 20, 30)
rate_cols <- paste0("Y", taus)

dt <- 1 / 12
pay_freq <- 1

h <- 6                    # 6-month-ahead forecast
rolling_window <- 60      # 5 years of monthly data

n_starts_oos <- 10
optim_maxit_oos <- 1500
jitter_scale <- 0.03

n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1)
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
doParallel::registerDoParallel(cl)
on.exit({
  try(parallel::stopCluster(cl), silent = TRUE)
  try(foreach::registerDoSEQ(), silent = TRUE)
  try(closeAllConnections(), silent = TRUE)
}, add = TRUE)
cat("Using", n_cores, "cores\n")
cat("Registered foreach workers:", foreach::getDoParWorkers(), "\n")

swap <- read_excel(swap_file)
names(swap)[1] <- "Date"

swap1 <- swap %>%
  mutate(Date = as.Date(Date)) %>%
  filter(!is.na(Date)) %>%
  filter(Date >= start_date) %>%
  arrange(Date)

Y_real <- as.matrix(swap1 %>% select(all_of(rate_cols))) / 100
dates  <- swap1$Date

stopifnot(!anyNA(Y_real), all(is.finite(Y_real)))

cat("Data loaded.\n")
cat("Number of observations:", nrow(Y_real), "\n")
cat("Range of swap rates:", range(Y_real), "\n\n")

###############################################################
# HELPERS
###############################################################

chol_pd <- function(S, jitter0 = 1e-10, max_tries = 12) {
  S <- (S + t(S)) / 2
  jitter <- jitter0
  for (i in seq_len(max_tries)) {
    out <- try(chol(S + diag(jitter, nrow(S))), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    jitter <- jitter * 10
  }
  stop("chol failed even after jitter inflation")
}

project_positive <- function(x, eps = 1e-10) pmax(x, eps)

jitter_transformed <- function(p, scale = 0.03) {
  # Works in transformed parameter space.
  # For Vasicek: p = (log kappa, log sigma, theta, log meas_sd)
  # For CIR:     p = (log kappa, log theta, log sigma, log meas_sd)
  out <- p
  out <- out + rnorm(length(p), mean = 0, sd = scale)
  out
}

###############################################################
# VASICEK 
###############################################################

vasicek_B <- function(tau, kappa) {
  (1 - exp(-kappa * tau)) / kappa
}

vasicek_A <- function(tau, kappa, theta, sigma) {
  B <- vasicek_B(tau, kappa)
  (theta - sigma^2 / (2 * kappa^2)) * (B - tau) -
    (sigma^2 / (4 * kappa)) * B^2
}

vasicek_P <- function(tau, r, kappa, theta, sigma) {
  A <- vasicek_A(tau, kappa, theta, sigma)
  B <- vasicek_B(tau, kappa)
  exp(A - B * r)
}

swap_rate_from_r_vas <- function(r, T, kappa, theta, sigma, pay_freq = 1) {
  delta <- 1 / pay_freq
  pay_times <- seq(delta, T, by = delta)
  Pn <- vasicek_P(T, r, kappa, theta, sigma)
  D  <- sum(delta * vasicek_P(pay_times, r, kappa, theta, sigma))
  (1 - Pn) / D
}

swap_rate_dr_vas <- function(r, T, kappa, theta, sigma, pay_freq = 1) {
  delta <- 1 / pay_freq
  pay_times <- seq(delta, T, by = delta)
  Pn <- vasicek_P(T, r, kappa, theta, sigma)
  Bn <- vasicek_B(T, kappa)
  dPn_dr <- -Bn * Pn
  Pk <- vasicek_P(pay_times, r, kappa, theta, sigma)
  Bk <- vasicek_B(pay_times, kappa)
  dPk_dr <- -Bk * Pk
  N <- 1 - Pn
  D <- sum(delta * Pk)
  dN_dr <- -dPn_dr
  dD_dr <- sum(delta * dPk_dr)
  (dN_dr * D - N * dD_dr) / D^2
}

kalman_vasicek_swaps_EKF <- function(
    Y, taus, dt, kappa, theta, sigma, meas_sd,
    pay_freq = 1, m10 = NULL, P10 = NULL, eps = 1e-12,
    return_states = TRUE
) {
  N <- nrow(Y)
  m <- ncol(Y)
  phi <- exp(-kappa * dt)
  a <- theta * (1 - phi)
  Q <- sigma^2 / (2 * kappa) * (1 - exp(-2 * kappa * dt))
  Q <- max(Q, eps)
  meas_sd <- max(meas_sd, 1e-8)
  R <- diag(rep(meas_sd^2, m), m, m)
  if (is.null(m10)) m10 <- theta
  if (is.null(P10)) P10 <- sigma^2 / (2 * kappa)
  mt <- m10
  Pt <- max(P10, eps)
  loglik <- 0
  if (return_states) {
    r_pred <- numeric(N)
    P_pred_store <- numeric(N)
    r_filt <- numeric(N)
    P_filt <- numeric(N)
    Y_pred <- matrix(NA_real_, N, m)
    Y_filt <- matrix(NA_real_, N, m)
    innovations <- matrix(NA_real_, N, m)
    colnames(Y_pred) <- paste0("Y", taus)
    colnames(Y_filt) <- paste0("Y", taus)
    colnames(innovations) <- paste0("Y", taus)
  }
  for (tt in seq_len(N)) {
    m_pred <- a + phi * mt
    P_pred <- phi^2 * Pt + Q
    P_pred <- max(P_pred, eps)
    yhat_pred <- sapply(taus, function(T) swap_rate_from_r_vas(m_pred, T, kappa, theta, sigma, pay_freq))
    Ht <- sapply(taus, function(T) swap_rate_dr_vas(m_pred, T, kappa, theta, sigma, pay_freq))
    yhat_pred <- matrix(yhat_pred, nrow = m, ncol = 1)
    Ht <- matrix(Ht, nrow = m, ncol = 1)
    y_t <- matrix(Y[tt, ], nrow = m, ncol = 1)
    v <- y_t - yhat_pred
    S <- Ht %*% t(Ht) * P_pred + R
    cholS <- chol_pd(S)
    Sinv_v <- backsolve(cholS, forwardsolve(t(cholS), v))
    Sinv_H <- backsolve(cholS, forwardsolve(t(cholS), Ht))
    logdetS <- 2 * sum(log(diag(cholS)))
    loglik <- loglik - 0.5 * (m * log(2 * pi) + logdetS + as.numeric(t(v) %*% Sinv_v))
    K <- as.numeric(P_pred) * t(Sinv_H)
    mt_new <- as.numeric(m_pred + K %*% v)
    KH <- as.numeric(K %*% Ht)
    Pt_new <- (1 - KH)^2 * P_pred + as.numeric(K %*% R %*% t(K))
    Pt_new <- max(Pt_new, eps)
    if (return_states) {
      yhat_filt <- sapply(taus, function(T) swap_rate_from_r_vas(mt_new, T, kappa, theta, sigma, pay_freq))
      r_pred[tt] <- m_pred
      P_pred_store[tt] <- P_pred
      r_filt[tt] <- mt_new
      P_filt[tt] <- Pt_new
      Y_pred[tt, ] <- as.numeric(yhat_pred)
      Y_filt[tt, ] <- as.numeric(yhat_filt)
      innovations[tt, ] <- as.numeric(v)
    }
    mt <- mt_new
    Pt <- Pt_new
  }
  if (return_states) {
    list(loglik = loglik, r_pred = r_pred, P_pred = P_pred_store,
         r_filt = r_filt, P_filt = P_filt, Y_pred = Y_pred,
         Y_filt = Y_filt, innovations = innovations)
  } else {
    list(loglik = loglik)
  }
}

negloglik_vas_trans <- function(p, Y_train) {
  kappa   <- exp(p[1]) + 1e-6
  sigma   <- exp(p[2]) + 1e-6
  theta   <- p[3]
  meas_sd <- exp(p[4]) + 1e-8
  val <- tryCatch({
    -kalman_vasicek_swaps_EKF(Y_train, taus, dt, kappa, theta, sigma, meas_sd,
                              pay_freq, return_states = FALSE)$loglik
  }, error = function(e) 1e12)
  if (!is.finite(val)) val <- 1e12
  val
}

fit_vasicek_window <- function(Y_train, main_start, n_starts = 10, maxit = 1500) {
  starts <- c(list(main_start), replicate(n_starts - 1, jitter_transformed(main_start, jitter_scale), simplify = FALSE))
  fits <- foreach(
    p0 = starts,
    .packages = c("stats"),
    .export = c(
      "negloglik_vas_trans", "kalman_vasicek_swaps_EKF", "chol_pd",
      "swap_rate_from_r_vas", "swap_rate_dr_vas",
      "vasicek_P", "vasicek_A", "vasicek_B",
      "Y_train", "taus", "dt", "pay_freq"
    ),
    .errorhandling = "pass"
  ) %dopar% {
    tryCatch(
      optim(p0, negloglik_vas_trans, Y_train = Y_train,
            method = "BFGS", control = list(maxit = maxit)),
      error = function(e) NULL
    )
  }
  fits <- Filter(function(x) !inherits(x, "error") && !is.null(x), fits)
  if (length(fits) == 0) return(NULL)
  vals <- sapply(fits, function(x) x$value)
  best <- fits[[which.min(vals)]]
  pars <- list(
    kappa = exp(best$par[1]) + 1e-6,
    theta = best$par[3],
    sigma = exp(best$par[2]) + 1e-6,
    meas_sd = exp(best$par[4]) + 1e-8,
    pay_freq = pay_freq
  )
  kf <- tryCatch(
    kalman_vasicek_swaps_EKF(Y_train, taus, dt, pars$kappa, pars$theta, pars$sigma,
                             pars$meas_sd, pay_freq, return_states = TRUE),
    error = function(e) NULL
  )
  if (is.null(kf)) return(NULL)
  list(opt = best, params = pars, filtered = kf, objective = best$value)
}

###############################################################
# CIR 
###############################################################

cir_gamma <- function(kappa, sigma) sqrt(kappa^2 + 2 * sigma^2)

cir_B <- function(tau, kappa, sigma) {
  gamma <- cir_gamma(kappa, sigma)
  numerator <- 2 * (exp(gamma * tau) - 1)
  denominator <- (gamma + kappa) * (exp(gamma * tau) - 1) + 2 * gamma
  numerator / denominator
}

cir_A <- function(tau, kappa, theta, sigma) {
  gamma <- cir_gamma(kappa, sigma)
  denominator <- (gamma + kappa) * (exp(gamma * tau) - 1) + 2 * gamma
  base <- (2 * gamma * exp((kappa + gamma) * tau / 2)) / denominator
  power <- 2 * kappa * theta / sigma^2
  base^power
}

cir_P <- function(tau, r, kappa, theta, sigma) {
  r <- pmax(r, 1e-12)
  cir_A(tau, kappa, theta, sigma) * exp(-cir_B(tau, kappa, sigma) * r)
}

swap_rate_from_r_cir <- function(r, T, kappa, theta, sigma, pay_freq = 1) {
  r <- max(r, 1e-12)
  delta <- 1 / pay_freq
  pay_times <- seq(delta, T, by = delta)
  Pn <- cir_P(T, r, kappa, theta, sigma)
  D <- sum(delta * cir_P(pay_times, r, kappa, theta, sigma))
  (1 - Pn) / D
}

swap_rate_dr_cir <- function(r, T, kappa, theta, sigma, pay_freq = 1) {
  r <- max(r, 1e-12)
  delta <- 1 / pay_freq
  pay_times <- seq(delta, T, by = delta)
  Pn <- cir_P(T, r, kappa, theta, sigma)
  Bn <- cir_B(T, kappa, sigma)
  dPn_dr <- -Bn * Pn
  Pk <- cir_P(pay_times, r, kappa, theta, sigma)
  Bk <- cir_B(pay_times, kappa, sigma)
  dPk_dr <- -Bk * Pk
  N <- 1 - Pn
  D <- sum(delta * Pk)
  dN_dr <- -dPn_dr
  dD_dr <- sum(delta * dPk_dr)
  (dN_dr * D - N * dD_dr) / D^2
}

cir_cond_var <- function(r, kappa, theta, sigma, dt) {
  r <- max(r, 1e-12)
  exp1 <- exp(-kappa * dt)
  exp2 <- exp(-2 * kappa * dt)
  term1 <- r * sigma^2 / kappa * (exp1 - exp2)
  term2 <- theta * sigma^2 / (2 * kappa) * (1 - exp1)^2
  term1 + term2
}

kalman_cir_swaps_EKF <- function(
    Y, taus, dt, kappa, theta, sigma, meas_sd,
    pay_freq = 1, m10 = NULL, P10 = NULL, eps = 1e-12,
    return_states = TRUE
) {
  N <- nrow(Y)
  m <- ncol(Y)
  phi <- exp(-kappa * dt)
  meas_sd <- max(meas_sd, 1e-8)
  R <- diag(rep(meas_sd^2, m), m, m)
  if (is.null(m10)) m10 <- theta
  if (is.null(P10)) P10 <- theta * sigma^2 / (2 * kappa)
  mt <- max(m10, eps)
  Pt <- max(P10, eps)
  loglik <- 0
  if (return_states) {
    r_pred <- numeric(N)
    P_pred_store <- numeric(N)
    r_filt <- numeric(N)
    P_filt <- numeric(N)
    Y_pred <- matrix(NA_real_, N, m)
    Y_filt <- matrix(NA_real_, N, m)
    innovations <- matrix(NA_real_, N, m)
    colnames(Y_pred) <- paste0("Y", taus)
    colnames(Y_filt) <- paste0("Y", taus)
    colnames(innovations) <- paste0("Y", taus)
  }
  for (tt in seq_len(N)) {
    m_pred <- theta + phi * (mt - theta)
    m_pred <- max(m_pred, eps)
    Q_pred <- cir_cond_var(mt, kappa, theta, sigma, dt)
    P_pred <- phi^2 * Pt + Q_pred
    P_pred <- max(P_pred, eps)
    yhat_pred <- sapply(taus, function(T) swap_rate_from_r_cir(m_pred, T, kappa, theta, sigma, pay_freq))
    Ht <- sapply(taus, function(T) swap_rate_dr_cir(m_pred, T, kappa, theta, sigma, pay_freq))
    yhat_pred <- matrix(yhat_pred, nrow = m, ncol = 1)
    Ht <- matrix(Ht, nrow = m, ncol = 1)
    y_t <- matrix(Y[tt, ], nrow = m, ncol = 1)
    v <- y_t - yhat_pred
    S <- Ht %*% t(Ht) * P_pred + R
    cholS <- chol_pd(S)
    Sinv_v <- backsolve(cholS, forwardsolve(t(cholS), v))
    Sinv_H <- backsolve(cholS, forwardsolve(t(cholS), Ht))
    logdetS <- 2 * sum(log(diag(cholS)))
    loglik <- loglik - 0.5 * (m * log(2 * pi) + logdetS + as.numeric(t(v) %*% Sinv_v))
    K <- as.numeric(P_pred) * t(Sinv_H)
    mt_new <- as.numeric(m_pred + K %*% v)
    mt_new <- max(mt_new, eps)
    KH <- as.numeric(K %*% Ht)
    Pt_new <- (1 - KH)^2 * P_pred + as.numeric(K %*% R %*% t(K))
    Pt_new <- max(Pt_new, eps)
    if (return_states) {
      yhat_filt <- sapply(taus, function(T) swap_rate_from_r_cir(mt_new, T, kappa, theta, sigma, pay_freq))
      r_pred[tt] <- m_pred
      P_pred_store[tt] <- P_pred
      r_filt[tt] <- mt_new
      P_filt[tt] <- Pt_new
      Y_pred[tt, ] <- as.numeric(yhat_pred)
      Y_filt[tt, ] <- as.numeric(yhat_filt)
      innovations[tt, ] <- as.numeric(v)
    }
    mt <- mt_new
    Pt <- Pt_new
  }
  if (return_states) {
    list(loglik = loglik, r_pred = r_pred, P_pred = P_pred_store,
         r_filt = r_filt, P_filt = P_filt, Y_pred = Y_pred,
         Y_filt = Y_filt, innovations = innovations)
  } else {
    list(loglik = loglik)
  }
}

negloglik_cir_trans <- function(p, Y_train) {
  kappa   <- exp(p[1]) + 1e-6
  theta   <- exp(p[2]) + 1e-8
  sigma   <- exp(p[3]) + 1e-6
  meas_sd <- exp(p[4]) + 1e-8
  val <- tryCatch({
    -kalman_cir_swaps_EKF(Y_train, taus, dt, kappa, theta, sigma, meas_sd,
                          pay_freq, return_states = FALSE)$loglik
  }, error = function(e) 1e12)
  if (!is.finite(val)) val <- 1e12
  val
}

fit_cir_window <- function(Y_train, main_start, n_starts = 10, maxit = 1500) {
  starts <- c(list(main_start), replicate(n_starts - 1, jitter_transformed(main_start, jitter_scale), simplify = FALSE))
  fits <- foreach(
    p0 = starts,
    .packages = c("stats"),
    .export = c(
      "negloglik_cir_trans", "kalman_cir_swaps_EKF", "chol_pd",
      "swap_rate_from_r_cir", "swap_rate_dr_cir", "cir_P", "cir_A", "cir_B",
      "cir_gamma", "cir_cond_var", "Y_train", "taus", "dt", "pay_freq"
    ),
    .errorhandling = "pass"
  ) %dopar% {
    tryCatch(
      optim(p0, negloglik_cir_trans, Y_train = Y_train,
            method = "BFGS", control = list(maxit = maxit)),
      error = function(e) NULL
    )
  }
  fits <- Filter(function(x) !inherits(x, "error") && !is.null(x), fits)
  if (length(fits) == 0) return(NULL)
  vals <- sapply(fits, function(x) x$value)
  best <- fits[[which.min(vals)]]
  pars <- list(
    kappa = exp(best$par[1]) + 1e-6,
    theta = exp(best$par[2]) + 1e-8,
    sigma = exp(best$par[3]) + 1e-6,
    meas_sd = exp(best$par[4]) + 1e-8,
    pay_freq = pay_freq,
    feller = 2 * (exp(best$par[1]) + 1e-6) * (exp(best$par[2]) + 1e-8) - (exp(best$par[3]) + 1e-6)^2
  )
  kf <- tryCatch(
    kalman_cir_swaps_EKF(Y_train, taus, dt, pars$kappa, pars$theta, pars$sigma,
                         pars$meas_sd, pay_freq, return_states = TRUE),
    error = function(e) NULL
  )
  if (is.null(kf)) return(NULL)
  list(opt = best, params = pars, filtered = kf, objective = best$value)
}

###############################################################
# FORECAST HELPERS
###############################################################

forecast_vas_h_steps <- function(last_r, pars, h) {
  r_h <- last_r
  phi <- exp(-pars$kappa * dt)
  for (i in seq_len(h)) {
    r_h <- pars$theta + phi * (r_h - pars$theta)
  }
  sapply(taus, function(T) swap_rate_from_r_vas(r_h, T, pars$kappa, pars$theta, pars$sigma, pars$pay_freq))
}

forecast_cir_h_steps <- function(last_r, pars, h) {
  r_h <- max(last_r, 1e-12)
  phi <- exp(-pars$kappa * dt)
  for (i in seq_len(h)) {
    r_h <- pars$theta + phi * (r_h - pars$theta)
    r_h <- max(r_h, 1e-12)
  }
  sapply(taus, function(T) swap_rate_from_r_cir(r_h, T, pars$kappa, pars$theta, pars$sigma, pars$pay_freq))
}

run_rolling_forecast <- function(model_label, fit_fun, forecast_fun, main_start, window = 60, h = 6) {
  t_grid <- window:(nrow(Y_real) - h)
  n_total <- length(t_grid)
  cat("\n====================================================\n")
  cat("Running 5-year rolling forecast for", model_label, "\n")
  cat("Windows:", n_total, "| Starts per window:", n_starts_oos, "| Horizon:", h, "months\n")
  cat("====================================================\n")
  res_list <- vector("list", n_total)
  global_start <- Sys.time()
  for (i in seq_along(t_grid)) {
    iter_start <- Sys.time()
    t <- t_grid[i]
    pct_done <- round(100 * (i - 1) / n_total, 1)
    cat("\n----------------------------------------------------\n")
    cat("Model:", model_label, "\n")
    cat("Iteration:", i, "of", n_total, "|", pct_done, "% complete before this window\n")
    cat("Training window:", as.character(dates[t - window + 1]), "to", as.character(dates[t]), "\n")
    cat("Forecast date:", as.character(dates[t + h]), "\n")
    cat("Started:", format(iter_start, "%Y-%m-%d %H:%M:%S"), "\n")
    cat("----------------------------------------------------\n")
    flush.console()
    set.seed(100000 + t + ifelse(model_label == "Vasicek", 1000, 2000))
    train_idx <- (t - window + 1):t
    Y_train <- Y_real[train_idx, , drop = FALSE]
    fit <- fit_fun(Y_train, main_start, n_starts_oos, optim_maxit_oos)
    if (!is.null(fit)) {
      last_r <- fit$filtered$r_filt[length(fit$filtered$r_filt)]
      pred_swap <- forecast_fun(last_r, fit$params, h)
      actual_swap <- Y_real[t + h, ]
      res_list[[i]] <- list(
        t = t,
        date = dates[t + h],
        forecast = pred_swap,
        actual = actual_swap,
        error = actual_swap - pred_swap,
        params = fit$params,
        objective = fit$objective
      )
      status_msg <- "SUCCESS"
    } else {
      res_list[[i]] <- NULL
      status_msg <- "FAILED - skipped"
    }
    iter_end <- Sys.time()
    iter_mins <- round(as.numeric(difftime(iter_end, iter_start, units = "mins")), 2)
    elapsed_mins <- as.numeric(difftime(iter_end, global_start, units = "mins"))
    avg_mins <- elapsed_mins / i
    remaining_mins <- round(avg_mins * (n_total - i), 1)
    pct_after <- round(100 * i / n_total, 1)
    cat("Status:", status_msg, "\n")
    cat("Iteration runtime:", iter_mins, "minutes\n")
    cat("Progress:", pct_after, "% complete\n")
    cat("Estimated remaining time:", remaining_mins, "minutes\n")
    flush.console()
  }
  res_list <- Filter(function(x) !is.null(x), res_list)
  if (length(res_list) == 0) stop("No successful forecasts for ", model_label)
  forecast_mat <- do.call(rbind, lapply(res_list, `[[`, "forecast"))
  actual_mat   <- do.call(rbind, lapply(res_list, `[[`, "actual"))
  error_mat    <- do.call(rbind, lapply(res_list, `[[`, "error"))
  colnames(forecast_mat) <- rate_cols
  colnames(actual_mat) <- rate_cols
  colnames(error_mat) <- rate_cols
  out <- list(
    model = model_label,
    method = "rolling",
    h = h,
    window = window,
    dates = as.Date(sapply(res_list, `[[`, "date"), origin = "1970-01-01"),
    forecast = forecast_mat,
    actual = actual_mat,
    errors = error_mat,
    rmse = sqrt(colMeans(error_mat^2, na.rm = TRUE)),
    mae = colMeans(abs(error_mat), na.rm = TRUE),
    objectives = sapply(res_list, `[[`, "objective"),
    params = lapply(res_list, `[[`, "params")
  )
  cat("\nFinished", model_label, "with", length(res_list), "successful forecasts out of", n_total, "windows.\n")
  out
}

###############################################################
# STARTING VALUES
###############################################################

# Estimated parameters supplied by the user:
# Vasicek: kappa=0.0132, theta=0.2125, sigma=0.0137, sigma_meas=0.0035
# CIR:     kappa=0.0073, theta=0.1608, sigma=0.0491, sigma_meas=0.0049

start_vas <- c(
  log(0.0132),  # log kappa
  log(0.0137),  # log sigma
  0.2125,       # theta, unrestricted for Vasicek
  log(0.0035)   # log measurement sd
)

start_cir <- c(
  log(0.0073),  # log kappa
  log(0.1608),  # log theta
  log(0.0491),  # log sigma
  log(0.0049)   # log measurement sd
)

###############################################################
# FORECASTS
###############################################################

all_results <- list()
all_results$rolling5y_vasicek <- run_rolling_forecast(
  model_label = "Vasicek",
  fit_fun = fit_vasicek_window,
  forecast_fun = forecast_vas_h_steps,
  main_start = start_vas,
  window = rolling_window,
  h = h
)

all_results$rolling5y_cir <- run_rolling_forecast(
  model_label = "CIR",
  fit_fun = fit_cir_window,
  forecast_fun = forecast_cir_h_steps,
  main_start = start_cir,
  window = rolling_window,
  h = h
)

###############################################################
# TABLES AND RESULTS
###############################################################

make_metric_table <- function(results, metric = c("rmse", "mae")) {
  metric <- match.arg(metric)
  tab <- data.frame(maturity = rate_cols)
  for (nm in names(results)) {
    tab[[nm]] <- as.numeric(results[[nm]][[metric]])
  }
  tab
}

rmse_table <- make_metric_table(all_results, "rmse")
mae_table  <- make_metric_table(all_results, "mae")

cat("\n================ RMSE TABLE ================\n")
print(rmse_table)
cat("\n================ MAE TABLE ================\n")
print(mae_table)

write.csv(rmse_table, "baseline_forecast_rmse_vasicek_cir.csv", row.names = FALSE)
write.csv(mae_table,  "baseline_forecast_mae_vasicek_cir.csv", row.names = FALSE)
saveRDS(all_results, "baseline_forecast_results_vasicek_cir.rds")

error_long <- do.call(rbind, lapply(names(all_results), function(nm) {
  res <- all_results[[nm]]
  data.frame(
    result = nm,
    model = res$model,
    method = res$method,
    date = rep(res$dates, times = length(rate_cols)),
    maturity = rep(rate_cols, each = length(res$dates)),
    error = as.vector(res$errors)
  )
}))
write.csv(error_long, "baseline_forecast_errors_long.csv", row.names = FALSE)

cat("\nSaved:\n")
cat("- baseline_forecast_rmse_vasicek_cir.csv\n")
cat("- baseline_forecast_mae_vasicek_cir.csv\n")
cat("- baseline_forecast_errors_long.csv\n")
cat("- baseline_forecast_results_vasicek_cir.rds\n")

###############################################################
# ACTUAL VS FORECAST PLOTS
###############################################################

make_plot_df <- function(results) {
  do.call(rbind, lapply(names(results), function(nm) {
    res <- results[[nm]]
    do.call(rbind, lapply(seq_along(rate_cols), function(j) {
      data.frame(
        result = nm,
        model = res$model,
        Date = res$dates,
        maturity = rate_cols[j],
        actual = res$actual[, j],
        forecast = res$forecast[, j]
      )
    }))
  }))
}

plot_df <- make_plot_df(all_results)
write.csv(plot_df, "baseline_forecast_actual_vs_forecast_long.csv", row.names = FALSE)

maturity_order <- rate_cols

actual_df <- plot_df %>%
  select(Date, maturity, actual) %>%
  distinct() %>%
  mutate(series = "Actual", value = actual)

forecast_df <- plot_df %>%
  transmute(Date, maturity, series = model, value = forecast)

plot_df_three_lines <- bind_rows(actual_df, forecast_df) %>%
  mutate(
    maturity = factor(maturity, levels = maturity_order),
    series = factor(series, levels = c("Actual", "Vasicek", "CIR"))
  )

p_three <- ggplot(
  plot_df_three_lines,
  aes(x = Date, y = value, color = series, linetype = series)
) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ maturity, scales = "free_y", ncol = 3) +
  scale_color_manual(
    values = c("Actual" = "black", "Vasicek" = "#ff69b4", "CIR" = "#1f77b4")
  ) +
  scale_linetype_manual(
    values = c("Actual" = "solid", "Vasicek" = "solid", "CIR" = "solid")
  ) +
  labs(
    title = "Actual vs 6-month-ahead forecasted EURIBOR swap rates",
    subtitle = "Vasicek and CIR baseline models, 5-year rolling estimation window",
    x = NULL,
    y = "Swap rate",
    color = NULL,
    linetype = NULL
  ) +
  theme_light() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_three)
ggsave("baseline_actual_vs_forecast_vasicek_cir.png", p_three, width = 13, height = 9, dpi = 300)

# Error plot
error_plot_df <- error_long %>%
  mutate(
    maturity = factor(maturity, levels = maturity_order),
    model = factor(model, levels = c("Vasicek", "CIR"))
  )

p_errors <- ggplot(error_plot_df, aes(x = date, y = error, color = model)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ maturity, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Vasicek" = "#ff69b4", "CIR" = "#1f77b4")) +
  labs(
    title = "6-month-ahead forecast errors",
    subtitle = "Error = actual - forecast",
    x = NULL,
    y = "Forecast error",
    color = NULL
  ) +
  theme_light() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_errors)
ggsave("baseline_forecast_errors_vasicek_cir.png", p_errors, width = 13, height = 9, dpi = 300)

###############################################################
# DONE
###############################################################

save.image("~/desktop/baselinemodeler_forecasting.RData")


