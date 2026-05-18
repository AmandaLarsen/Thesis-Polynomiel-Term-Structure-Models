###############################################################
# COMBINED SCRIPT:
# Vasicek + CIR 1-factor EKF directly on EURIBOR swap rates
#
# Includes:
# - Real-data estimation for Vasicek and CIR
# - One-step-ahead predictions
# - Filtered reconstructions
# - Real-data residual plots
# - Real-data RMSE plots
# - Combined Vasicek vs CIR comparison plots
# - Simulation / recovery study for Vasicek
# - Simulation / recovery study for CIR
# - Simulation state plots and filtering errors
# - Simulation RMSE plots
# - Jacobian sanity checks
###############################################################

library(readxl)
library(dplyr)
library(tidyr)
library(parallel)

###############################################################
# 0) USER SETTINGS
###############################################################

swap_file <- "~/Desktop/ægte data/EURIBOR swap data.xlsx"
start_date <- as.Date("2008-10-06")

taus <- c(1, 2, 3, 5, 7, 10, 15, 20, 30)
rate_cols <- paste0("Y", taus)

dt <- 1 / 12
pay_freq <- 1

n_starts <- 50
n_cores  <- 7

set.seed(123)

###############################################################
# 1) LOAD DATA
###############################################################

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
# 2) COMMON HELPERS
###############################################################

chol_pd <- function(S, jitter0 = 1e-10, max_tries = 12) {
  S <- (S + t(S)) / 2
  jitter <- jitter0
  
  for (i in 1:max_tries) {
    out <- try(chol(S + diag(jitter, nrow(S))), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    jitter <- jitter * 10
  }
  
  stop("chol failed even after jitter inflation")
}

rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))

plot_swap_grid <- function(x_axis, Y_obs, Y_model, title, model_col = "red",
                           obs_label = "Observed", model_label = "Model") {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  par(mfrow = c(3, 3), mar = c(3, 3, 3, 1), oma = c(0, 0, 3, 0))
  
  for (j in seq_along(taus)) {
    y_min <- min(Y_obs[, j], Y_model[, j], na.rm = TRUE)
    y_max <- max(Y_obs[, j], Y_model[, j], na.rm = TRUE)
    
    plot(
      x_axis, Y_obs[, j],
      type = "l",
      main = paste0(taus[j], "Y"),
      xlab = "",
      ylab = "",
      lwd = 1.5,
      col = "black",
      ylim = c(y_min, y_max)
    )
    
    lines(x_axis, Y_model[, j], lty = 2, lwd = 1.5, col = model_col)
  }
  
  mtext(title, outer = TRUE, cex = 1.4)
  
  legend(
    "bottom",
    inset = -0.02,
    legend = c(obs_label, model_label),
    col = c("black", model_col),
    lty = c(1, 2),
    lwd = c(1.5, 1.5),
    horiz = TRUE,
    bty = "n",
    xpd = TRUE
  )
}

plot_residual_grid <- function(x_axis, residuals, title, resid_col = "black",
                               start_idx = 3) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  idx <- start_idx:nrow(residuals)
  
  par(mfrow = c(3, 3), mar = c(3, 3, 3, 1), oma = c(0, 0, 3, 0))
  
  for (j in seq_along(taus)) {
    plot(
      x_axis[idx], residuals[idx, j],
      type = "l",
      main = paste0(taus[j], "Y"),
      xlab = "",
      ylab = "",
      lwd = 1.2,
      col = resid_col
    )
    abline(h = 0, lty = 2, col = "red", lwd = 1.2)
  }
  
  mtext(title, outer = TRUE, cex = 1.4)
}

plot_short_rate <- function(x_axis, r_filt, r_pred, title) {
  plot(
    x_axis, r_filt,
    type = "l",
    xlab = "",
    ylab = "r_t",
    main = title,
    lwd = 1.5,
    col = "black"
  )
  lines(x_axis, r_pred, lty = 2, lwd = 1.5, col = "red")
  legend(
    "topright",
    legend = c("Filtered r_{t|t}", "Predicted r_{t|t-1}"),
    col = c("black", "red"),
    lty = c(1, 2),
    lwd = 1.5,
    bty = "n"
  )
}

plot_rmse_two_lines <- function(taus, rmse_1, rmse_2, title,
                                label_1 = "One-step-ahead", label_2 = "Filtered",
                                col_1 = "black", col_2 = "red") {
  plot(
    taus, rmse_1,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "Maturity",
    ylab = "RMSE",
    main = title,
    col = col_1,
    ylim = range(rmse_1, rmse_2, na.rm = TRUE)
  )
  lines(taus, rmse_2, type = "b", pch = 19, lwd = 2, col = col_2)
  legend(
    "topright",
    legend = c(label_1, label_2),
    col = c(col_1, col_2),
    lty = 1,
    pch = 19,
    lwd = 2,
    bty = "n"
  )
  grid()
}

###############################################################
# 3) VASICEK FUNCTIONS
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
    Y, taus, dt,
    kappa, theta, sigma, meas_sd,
    pay_freq = 1,
    m10 = NULL, P10 = NULL,
    eps = 1e-12,
    return_states = TRUE
) {
  N <- nrow(Y)
  m <- ncol(Y)
  
  phi <- exp(-kappa * dt)
  a   <- theta * (1 - phi)
  Q   <- sigma^2 / (2 * kappa) * (1 - exp(-2 * kappa * dt))
  Q   <- max(Q, eps)
  
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
  
  for (tt in 1:N) {
    m_pred <- a + phi * mt
    P_pred <- phi^2 * Pt + Q
    P_pred <- max(P_pred, eps)
    
    yhat_pred <- sapply(taus, function(T) {
      swap_rate_from_r_vas(m_pred, T, kappa, theta, sigma, pay_freq)
    })
    
    Ht <- sapply(taus, function(T) {
      swap_rate_dr_vas(m_pred, T, kappa, theta, sigma, pay_freq)
    })
    
    yhat_pred <- matrix(yhat_pred, nrow = m, ncol = 1)
    Ht <- matrix(Ht, nrow = m, ncol = 1)
    y_t <- matrix(Y[tt, ], nrow = m, ncol = 1)
    v <- y_t - yhat_pred
    
    S <- Ht %*% t(Ht) * P_pred + R
    cholS <- chol_pd(S)
    
    Sinv_v <- backsolve(cholS, forwardsolve(t(cholS), v))
    Sinv_H <- backsolve(cholS, forwardsolve(t(cholS), Ht))
    logdetS <- 2 * sum(log(diag(cholS)))
    
    loglik <- loglik - 0.5 * (
      m * log(2 * pi) + logdetS + as.numeric(t(v) %*% Sinv_v)
    )
    
    K <- as.numeric(P_pred) * t(Sinv_H)
    mt_new <- as.numeric(m_pred + K %*% v)
    
    KH <- as.numeric(K %*% Ht)
    Pt_new <- (1 - KH)^2 * P_pred + as.numeric(K %*% R %*% t(K))
    Pt_new <- max(Pt_new, eps)
    
    if (return_states) {
      yhat_filt <- sapply(taus, function(T) {
        swap_rate_from_r_vas(mt_new, T, kappa, theta, sigma, pay_freq)
      })
      
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
    list(
      loglik = loglik,
      r_pred = r_pred,
      P_pred = P_pred_store,
      r_filt = r_filt,
      P_filt = P_filt,
      Y_pred = Y_pred,
      Y_filt = Y_filt,
      innovations = innovations
    )
  } else {
    list(loglik = loglik)
  }
}

make_start_values_vas <- function(n_starts = 50) {
  base <- list(
    c(log(0.30), log(0.010), 0.030, log(0.0010)),
    c(log(0.10), log(0.005), 0.020, log(0.0015)),
    c(log(0.60), log(0.015), 0.040, log(0.0008)),
    c(log(1.00), log(0.020), 0.030, log(0.0020)),
    c(log(0.20), log(0.008), 0.010, log(0.0010))
  )
  
  if (n_starts <= length(base)) return(base[1:n_starts])
  
  extra <- replicate(
    n_starts - length(base),
    c(
      log(runif(1, 0.05, 1.50)),
      log(runif(1, 0.002, 0.050)),
      runif(1, -0.02, 0.08),
      log(runif(1, 0.0002, 0.005))
    ),
    simplify = FALSE
  )
  
  c(base, extra)
}

fit_one_start_vas <- function(p0, Y, taus, dt, pay_freq = 1) {
  negloglik <- function(p) {
    kappa   <- exp(p[1]) + 1e-6
    sigma   <- exp(p[2]) + 1e-6
    theta   <- p[3]
    meas_sd <- exp(p[4]) + 1e-8
    
    val <- tryCatch({
      -kalman_vasicek_swaps_EKF(
        Y, taus, dt,
        kappa, theta, sigma, meas_sd,
        pay_freq,
        return_states = FALSE
      )$loglik
    }, error = function(e) 1e12)
    
    if (!is.finite(val)) val <- 1e12
    val
  }
  
  opt <- optim(p0, negloglik, method = "BFGS", control = list(maxit = 3000))
  
  list(
    opt = opt,
    value = opt$value,
    convergence = opt$convergence,
    params = list(
      kappa = exp(opt$par[1]) + 1e-6,
      theta = opt$par[3],
      sigma = exp(opt$par[2]) + 1e-6,
      meas_sd = exp(opt$par[4]) + 1e-8,
      pay_freq = pay_freq
    )
  )
}

fit_vasicek_EKF_multistart <- function(Y, taus, dt, pay_freq = 1,
                                       n_starts = 50, n_cores = 7) {
  starts <- make_start_values_vas(n_starts)
  n_cores <- min(n_cores, n_starts)
  
  cat("Running", n_starts, "Vasicek MLE starts on", n_cores, "cores...\n")
  
  if (.Platform$OS.type == "windows") {
    cl <- makeCluster(n_cores)
    clusterExport(
      cl,
      varlist = c(
        "Y", "taus", "dt", "pay_freq",
        "fit_one_start_vas", "kalman_vasicek_swaps_EKF", "chol_pd",
        "swap_rate_from_r_vas", "swap_rate_dr_vas",
        "vasicek_P", "vasicek_A", "vasicek_B"
      ),
      envir = environment()
    )
    fits <- parLapply(cl, starts, function(p0) {
      fit_one_start_vas(p0, Y, taus, dt, pay_freq)
    })
    stopCluster(cl)
  } else {
    fits <- mclapply(
      starts,
      function(p0) fit_one_start_vas(p0, Y, taus, dt, pay_freq),
      mc.cores = n_cores
    )
  }
  
  values <- sapply(fits, function(x) x$value)
  best_id <- which.min(values)
  best <- fits[[best_id]]
  
  cat("Best Vasicek start:", best_id, "\n")
  cat("Best Vasicek negative loglik:", best$value, "\n")
  cat("Convergence:", best$convergence, "\n")
  print(best$params)
  
  kf <- kalman_vasicek_swaps_EKF(
    Y = Y,
    taus = taus,
    dt = dt,
    kappa = best$params$kappa,
    theta = best$params$theta,
    sigma = best$params$sigma,
    meas_sd = best$params$meas_sd,
    pay_freq = pay_freq,
    return_states = TRUE
  )
  
  list(
    params = best$params,
    filtered = kf,
    opt = best$opt,
    all_fits = fits,
    all_values = values,
    best_id = best_id,
    starts = starts
  )
}

###############################################################
# 4) CIR FUNCTIONS
###############################################################

cir_gamma <- function(kappa, sigma) {
  sqrt(kappa^2 + 2 * sigma^2)
}

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
  D  <- sum(delta * cir_P(pay_times, r, kappa, theta, sigma))
  
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
    Y, taus, dt,
    kappa, theta, sigma, meas_sd,
    pay_freq = 1,
    m10 = NULL, P10 = NULL,
    eps = 1e-12,
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
  
  for (tt in 1:N) {
    m_pred <- theta + phi * (mt - theta)
    m_pred <- max(m_pred, eps)
    
    Q_pred <- cir_cond_var(mt, kappa, theta, sigma, dt)
    P_pred <- phi^2 * Pt + Q_pred
    P_pred <- max(P_pred, eps)
    
    yhat_pred <- sapply(taus, function(T) {
      swap_rate_from_r_cir(m_pred, T, kappa, theta, sigma, pay_freq)
    })
    
    Ht <- sapply(taus, function(T) {
      swap_rate_dr_cir(m_pred, T, kappa, theta, sigma, pay_freq)
    })
    
    yhat_pred <- matrix(yhat_pred, nrow = m, ncol = 1)
    Ht <- matrix(Ht, nrow = m, ncol = 1)
    y_t <- matrix(Y[tt, ], nrow = m, ncol = 1)
    v <- y_t - yhat_pred
    
    S <- Ht %*% t(Ht) * P_pred + R
    cholS <- chol_pd(S)
    
    Sinv_v <- backsolve(cholS, forwardsolve(t(cholS), v))
    Sinv_H <- backsolve(cholS, forwardsolve(t(cholS), Ht))
    logdetS <- 2 * sum(log(diag(cholS)))
    
    loglik <- loglik - 0.5 * (
      m * log(2 * pi) + logdetS + as.numeric(t(v) %*% Sinv_v)
    )
    
    K <- as.numeric(P_pred) * t(Sinv_H)
    mt_new <- as.numeric(m_pred + K %*% v)
    mt_new <- max(mt_new, eps)
    
    KH <- as.numeric(K %*% Ht)
    Pt_new <- (1 - KH)^2 * P_pred + as.numeric(K %*% R %*% t(K))
    Pt_new <- max(Pt_new, eps)
    
    if (return_states) {
      yhat_filt <- sapply(taus, function(T) {
        swap_rate_from_r_cir(mt_new, T, kappa, theta, sigma, pay_freq)
      })
      
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
    list(
      loglik = loglik,
      r_pred = r_pred,
      P_pred = P_pred_store,
      r_filt = r_filt,
      P_filt = P_filt,
      Y_pred = Y_pred,
      Y_filt = Y_filt,
      innovations = innovations
    )
  } else {
    list(loglik = loglik)
  }
}

make_start_values_cir <- function(n_starts = 50) {
  base <- list(
    c(log(0.30), log(0.030), log(0.080), log(0.0010)),
    c(log(0.10), log(0.020), log(0.050), log(0.0015)),
    c(log(0.60), log(0.040), log(0.100), log(0.0008)),
    c(log(1.00), log(0.030), log(0.120), log(0.0020)),
    c(log(0.20), log(0.010), log(0.060), log(0.0010))
  )
  
  if (n_starts <= length(base)) return(base[1:n_starts])
  
  extra <- replicate(
    n_starts - length(base),
    c(
      log(runif(1, 0.05, 1.50)),
      log(runif(1, 0.005, 0.080)),
      log(runif(1, 0.020, 0.150)),
      log(runif(1, 0.0002, 0.005))
    ),
    simplify = FALSE
  )
  
  c(base, extra)
}

fit_one_start_cir <- function(p0, Y, taus, dt, pay_freq = 1) {
  negloglik <- function(p) {
    kappa   <- exp(p[1]) + 1e-6
    theta   <- exp(p[2]) + 1e-8
    sigma   <- exp(p[3]) + 1e-6
    meas_sd <- exp(p[4]) + 1e-8
    
    val <- tryCatch({
      -kalman_cir_swaps_EKF(
        Y, taus, dt,
        kappa, theta, sigma, meas_sd,
        pay_freq,
        return_states = FALSE
      )$loglik
    }, error = function(e) 1e12)
    
    if (!is.finite(val)) val <- 1e12
    val
  }
  
  opt <- optim(p0, negloglik, method = "BFGS", control = list(maxit = 3000))
  
  kappa_hat   <- exp(opt$par[1]) + 1e-6
  theta_hat   <- exp(opt$par[2]) + 1e-8
  sigma_hat   <- exp(opt$par[3]) + 1e-6
  meas_sd_hat <- exp(opt$par[4]) + 1e-8
  
  list(
    opt = opt,
    value = opt$value,
    convergence = opt$convergence,
    params = list(
      kappa = kappa_hat,
      theta = theta_hat,
      sigma = sigma_hat,
      meas_sd = meas_sd_hat,
      pay_freq = pay_freq,
      feller = 2 * kappa_hat * theta_hat - sigma_hat^2
    )
  )
}

fit_cir_EKF_multistart <- function(Y, taus, dt, pay_freq = 1,
                                   n_starts = 50, n_cores = 7) {
  starts <- make_start_values_cir(n_starts)
  n_cores <- min(n_cores, n_starts)
  
  cat("Running", n_starts, "CIR MLE starts on", n_cores, "cores...\n")
  
  if (.Platform$OS.type == "windows") {
    cl <- makeCluster(n_cores)
    clusterExport(
      cl,
      varlist = c(
        "Y", "taus", "dt", "pay_freq",
        "fit_one_start_cir", "kalman_cir_swaps_EKF", "chol_pd",
        "swap_rate_from_r_cir", "swap_rate_dr_cir",
        "cir_P", "cir_A", "cir_B", "cir_gamma", "cir_cond_var"
      ),
      envir = environment()
    )
    fits <- parLapply(cl, starts, function(p0) {
      fit_one_start_cir(p0, Y, taus, dt, pay_freq)
    })
    stopCluster(cl)
  } else {
    fits <- mclapply(
      starts,
      function(p0) fit_one_start_cir(p0, Y, taus, dt, pay_freq),
      mc.cores = n_cores
    )
  }
  
  values <- sapply(fits, function(x) x$value)
  best_id <- which.min(values)
  best <- fits[[best_id]]
  
  cat("Best CIR start:", best_id, "\n")
  cat("Best CIR negative loglik:", best$value, "\n")
  cat("Convergence:", best$convergence, "\n")
  print(best$params)
  
  kf <- kalman_cir_swaps_EKF(
    Y = Y,
    taus = taus,
    dt = dt,
    kappa = best$params$kappa,
    theta = best$params$theta,
    sigma = best$params$sigma,
    meas_sd = best$params$meas_sd,
    pay_freq = pay_freq,
    return_states = TRUE
  )
  
  list(
    params = best$params,
    filtered = kf,
    opt = best$opt,
    all_fits = fits,
    all_values = values,
    best_id = best_id,
    starts = starts
  )
}

###############################################################
# 5) FIT REAL DATA: VASICEK + CIR
###############################################################

res_vas_real <- fit_vasicek_EKF_multistart(
  Y = Y_real,
  taus = taus,
  dt = dt,
  pay_freq = pay_freq,
  n_starts = n_starts,
  n_cores = n_cores
)

res_cir_real <- fit_cir_EKF_multistart(
  Y = Y_real,
  taus = taus,
  dt = dt,
  pay_freq = pay_freq,
  n_starts = n_starts,
  n_cores = n_cores
)

cat("\n=== REAL DATA: VASICEK PARAMS ===\n")
print(res_vas_real$params)

cat("\n=== REAL DATA: CIR PARAMS ===\n")
print(res_cir_real$params)

cat("\nCIR Feller check 2*kappa*theta - sigma^2:\n")
print(res_cir_real$params$feller)

###############################################################
# 6) EXTRACT REAL DATA RESULTS
###############################################################

Y_pred_vas_real <- res_vas_real$filtered$Y_pred
Y_filt_vas_real <- res_vas_real$filtered$Y_filt
Y_pred_cir_real <- res_cir_real$filtered$Y_pred
Y_filt_cir_real <- res_cir_real$filtered$Y_filt

resid_pred_vas_real <- Y_real - Y_pred_vas_real
resid_filt_vas_real <- Y_real - Y_filt_vas_real
resid_pred_cir_real <- Y_real - Y_pred_cir_real
resid_filt_cir_real <- Y_real - Y_filt_cir_real

cat("\nREAL DATA: Vasicek one-step-ahead residual SD:\n")
print(setNames(apply(resid_pred_vas_real, 2, sd), paste0(taus, "Y")))

cat("\nREAL DATA: CIR one-step-ahead residual SD:\n")
print(setNames(apply(resid_pred_cir_real, 2, sd), paste0(taus, "Y")))

###############################################################
# 7) REAL DATA PLOTS: SEPARATE MODELS
###############################################################

plot_swap_grid(
  dates, Y_real, Y_pred_vas_real,
  "VASICEK REAL DATA: observed vs ONE-STEP-AHEAD predicted swaps",
  model_col = "red", model_label = "Vasicek predicted"
)

plot_swap_grid(
  dates, Y_real, Y_filt_vas_real,
  "VASICEK REAL DATA: observed vs FILTERED reconstructed swaps",
  model_col = "red", model_label = "Vasicek filtered"
)

plot_short_rate(
  dates,
  res_vas_real$filtered$r_filt,
  res_vas_real$filtered$r_pred,
  "VASICEK real data: predicted vs filtered short rate"
)

plot_residual_grid(
  dates, resid_pred_vas_real,
  "VASICEK REAL DATA: one-step-ahead residuals",
  resid_col = "black"
)

plot_residual_grid(
  dates, resid_filt_vas_real,
  "VASICEK REAL DATA: filtered reconstruction residuals",
  resid_col = "black"
)

plot_swap_grid(
  dates, Y_real, Y_pred_cir_real,
  "CIR REAL DATA: observed vs ONE-STEP-AHEAD predicted swaps",
  model_col = "blue", model_label = "CIR predicted"
)

plot_swap_grid(
  dates, Y_real, Y_filt_cir_real,
  "CIR REAL DATA: observed vs FILTERED reconstructed swaps",
  model_col = "blue", model_label = "CIR filtered"
)

plot_short_rate(
  dates,
  res_cir_real$filtered$r_filt,
  res_cir_real$filtered$r_pred,
  "CIR real data: predicted vs filtered short rate"
)

plot_residual_grid(
  dates, resid_pred_cir_real,
  "CIR REAL DATA: one-step-ahead residuals",
  resid_col = "black"
)

plot_residual_grid(
  dates, resid_filt_cir_real,
  "CIR REAL DATA: filtered reconstruction residuals",
  resid_col = "black"
)

###############################################################
# 8) REAL DATA PLOTS: COMBINED VASICEK VS CIR
###############################################################

op <- par(no.readonly = TRUE)
par(mfrow = c(3, 3), mar = c(3, 3, 3, 1), oma = c(0, 0, 3, 0))

for (j in seq_along(taus)) {
  y_min <- min(Y_real[, j], Y_pred_vas_real[, j], Y_pred_cir_real[, j], na.rm = TRUE)
  y_max <- max(Y_real[, j], Y_pred_vas_real[, j], Y_pred_cir_real[, j], na.rm = TRUE)
  
  plot(
    dates, Y_real[, j],
    type = "l",
    main = paste0(taus[j], "Y"),
    xlab = "",
    ylab = "",
    lwd = 1.5,
    col = "black",
    ylim = c(y_min, y_max)
  )
  lines(dates, Y_pred_vas_real[, j], lty = 2, lwd = 1.5, col = "red")
  lines(dates, Y_pred_cir_real[, j], lty = 2, lwd = 1.5, col = "blue")
}

mtext("REAL DATA: observed vs one-step-ahead predicted swaps", outer = TRUE, cex = 1.4)
legend(
  "bottom",
  inset = -0.02,
  legend = c("Observed", "Vasicek", "CIR"),
  col = c("black", "red", "blue"),
  lty = c(1, 2, 2),
  lwd = c(1.5, 1.5, 1.5),
  horiz = TRUE,
  bty = "n",
  xpd = TRUE
)
par(op)

op <- par(no.readonly = TRUE)
par(mfrow = c(3, 3), mar = c(3, 3, 3, 1), oma = c(0, 0, 3, 0))

for (j in seq_along(taus)) {
  y_min <- min(Y_real[, j], Y_filt_vas_real[, j], Y_filt_cir_real[, j], na.rm = TRUE)
  y_max <- max(Y_real[, j], Y_filt_vas_real[, j], Y_filt_cir_real[, j], na.rm = TRUE)
  
  plot(
    dates, Y_real[, j],
    type = "l",
    main = paste0(taus[j], "Y"),
    xlab = "",
    ylab = "",
    lwd = 1.5,
    col = "black",
    ylim = c(y_min, y_max)
  )
  lines(dates, Y_filt_vas_real[, j], lty = 2, lwd = 1.5, col = "red")
  lines(dates, Y_filt_cir_real[, j], lty = 2, lwd = 1.5, col = "blue")
}

#################### NYE plots ########



# Fælles plot-funktion
library(tidyverse)
library(patchwork)

# Fælles plot-funktion
library(tidyverse)
library(patchwork)

plot_fit_terms_gg <- function(selected_terms, plot_title = "") {
  
  term_idx <- match(selected_terms, taus)
  term_idx <- term_idx[!is.na(term_idx)]
  
  plot_df <- bind_rows(lapply(term_idx, function(j) {
    tibble(
      Date = dates,
      Term = paste0(taus[j], "Y"),
      Observed = Y_real[, j],
      Vasicek = Y_filt_vas_real[, j],
      CIR = Y_filt_cir_real[, j]
    )
  })) %>%
    pivot_longer(
      cols = c(Observed, Vasicek, CIR),
      names_to = "Model",
      values_to = "Rate"
    ) %>%
    mutate(
      Term = factor(Term, levels = paste0(selected_terms, "Y")),
      Model = factor(Model, levels = c("Observed", "Vasicek", "CIR"))
    )
  
  ggplot(plot_df, aes(x = Date, y = Rate, color = Model, linetype = Model)) +
    geom_line(linewidth = 0.7) +
    
    facet_wrap(~ Term, ncol = 2, scales = "free_y") +
    
    scale_color_manual(
      values = c(
        "Observed" = "black",
        "Vasicek"  = "#1f77b4",  # blue
        "CIR"      = "#ff7f0e"   # orange
      )
    ) +
    
    scale_linetype_manual(
      values = c(
        "Observed" = "solid",
        "Vasicek"  = "dashed",
        "CIR"      = "dashed"
      )
    ) +
    
    scale_x_date(
      date_breaks = "2 years",
      date_labels = "%Y"
    ) +
    
    labs(
      title = plot_title,
      x = "Date",
      y = "Swap rate",
      color = NULL,
      linetype = NULL
    ) +
    
    theme_classic(base_size = 18) +
    
    theme_bw(base_size = 18) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 22,
        hjust = 0.5
      ),
      
      strip.background = element_blank(),
      
      strip.text = element_text(
        face = "bold",
        size = 10
      ),,
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 12
      ),
      
      axis.text.y = element_text(size = 14),
      
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      
      legend.position = "bottom",
      legend.text = element_text(size = 16),
      
      panel.grid.major = element_line(
        color = "grey85",
        linewidth = 0.4
      ),
      
      panel.grid.minor = element_blank()
    )
}



plot_fit_terms_gg(
  selected_terms = c(1, 3, 5, 10, 20, 30)
)

plot_fit_terms_gg(
  selected_terms = c(1, 2, 3, 5, 7, 10, 20, 30)
)

plot_fit_terms_gg(
  selected_terms = c(1, 2, 3, 5, 7, 10, 15, 20, 30)
)


###################

mtext("REAL DATA: observed vs filtered reconstructed swaps", outer = TRUE, cex = 1.4)
legend(
  "bottom",
  inset = -0.02,
  legend = c("Observed", "Vasicek", "CIR"),
  col = c("black", "red", "blue"),
  lty = c(1, 2, 2),
  lwd = c(1.5, 1.5, 1.5),
  horiz = TRUE,
  bty = "n",
  xpd = TRUE
)
par(op)

start_idx <- 3
op <- par(no.readonly = TRUE)
par(mfrow = c(3, 3), mar = c(3, 3, 3, 1), oma = c(0, 0, 3, 0))

for (j in seq_along(taus)) {
  idx <- start_idx:nrow(Y_real)
  y_min <- min(resid_pred_vas_real[idx, j], resid_pred_cir_real[idx, j], na.rm = TRUE)
  y_max <- max(resid_pred_vas_real[idx, j], resid_pred_cir_real[idx, j], na.rm = TRUE)
  
  plot(
    dates[idx], resid_pred_vas_real[idx, j],
    type = "l",
    main = paste0(taus[j], "Y"),
    xlab = "",
    ylab = "",
    lwd = 1.2,
    col = "red",
    lty = 2,
    ylim = c(y_min, y_max)
  )
  lines(dates[idx], resid_pred_cir_real[idx, j], lwd = 1.2, col = "blue", lty = 2)
  abline(h = 0, lty = 1, col = "black")
}

mtext("REAL DATA: one-step-ahead residuals, Vasicek vs CIR", outer = TRUE, cex = 1.4)
legend(
  "bottom",
  inset = -0.02,
  legend = c("Vasicek residual", "CIR residual", "Zero"),
  col = c("red", "blue", "black"),
  lty = c(2, 2, 1),
  lwd = c(1.2, 1.2, 1),
  horiz = TRUE,
  bty = "n",
  xpd = TRUE
)
par(op)

rmse_idx <- 3:nrow(Y_real)
rmse_pred_vas_real <- sqrt(colMeans(resid_pred_vas_real[rmse_idx, ]^2))
rmse_filt_vas_real <- sqrt(colMeans(resid_filt_vas_real[rmse_idx, ]^2))
rmse_pred_cir_real <- sqrt(colMeans(resid_pred_cir_real[rmse_idx, ]^2))
rmse_filt_cir_real <- sqrt(colMeans(resid_filt_cir_real[rmse_idx, ]^2))

plot_rmse_two_lines(
  taus, rmse_pred_vas_real, rmse_filt_vas_real,
  "VASICEK REAL DATA: RMSE by maturity",
  label_1 = "One-step-ahead", label_2 = "Filtered",
  col_1 = "black", col_2 = "red"
)

plot_rmse_two_lines(
  taus, rmse_pred_cir_real, rmse_filt_cir_real,
  "CIR REAL DATA: RMSE by maturity",
  label_1 = "One-step-ahead", label_2 = "Filtered",
  col_1 = "black", col_2 = "blue"
)

plot_rmse_two_lines(
  taus, rmse_pred_vas_real, rmse_pred_cir_real,
  "REAL DATA: one-step-ahead RMSE, Vasicek vs CIR",
  label_1 = "Vasicek", label_2 = "CIR",
  col_1 = "red", col_2 = "blue"
)

###############################################################
# 9) VASICEK SIMULATION STUDY
###############################################################

simulate_vasicek_exact <- function(r0, kappa, theta, sigma, dt, nSteps, seed = 1) {
  set.seed(seed)
  
  phi <- exp(-kappa * dt)
  a <- theta * (1 - phi)
  sd_w <- sigma * sqrt((1 - exp(-2 * kappa * dt)) / (2 * kappa))
  
  r <- numeric(nSteps + 1)
  r[1] <- r0
  
  for (tt in 1:nSteps) {
    r[tt + 1] <- a + phi * r[tt] + sd_w * rnorm(1)
  }
  
  r
}

simulate_swap_panel_vas <- function(r_path, taus, kappa, theta, sigma,
                                    meas_sd, pay_freq = 1, seed = 2) {
  set.seed(seed)
  
  N <- length(r_path) - 1
  m <- length(taus)
  Y <- matrix(NA_real_, nrow = N, ncol = m)
  colnames(Y) <- paste0("Y", taus)
  
  for (tt in 1:N) {
    r_t <- r_path[tt + 1]
    s_true <- sapply(taus, function(T) {
      swap_rate_from_r_vas(r_t, T, kappa, theta, sigma, pay_freq)
    })
    Y[tt, ] <- s_true + rnorm(m, 0, meas_sd)
  }
  
  Y
}

vas_true <- list(
  kappa = 0.30,
  theta = 0.03,
  sigma = 0.01,
  meas_sd = 0.001,
  r0 = 0.03
)

N_sim <- nrow(Y_real)

r_sim_vas <- simulate_vasicek_exact(
  r0 = vas_true$r0,
  kappa = vas_true$kappa,
  theta = vas_true$theta,
  sigma = vas_true$sigma,
  dt = dt,
  nSteps = N_sim,
  seed = 123
)

Y_sim_vas <- simulate_swap_panel_vas(
  r_path = r_sim_vas,
  taus = taus,
  kappa = vas_true$kappa,
  theta = vas_true$theta,
  sigma = vas_true$sigma,
  meas_sd = vas_true$meas_sd,
  pay_freq = pay_freq,
  seed = 999
)

res_vas_sim <- fit_vasicek_EKF_multistart(
  Y = Y_sim_vas,
  taus = taus,
  dt = dt,
  pay_freq = pay_freq,
  n_starts = n_starts,
  n_cores = n_cores
)

cat("\n=== VASICEK SIMULATION RECOVERY TEST ===\n")
print(cbind(
  true = c(
    kappa = vas_true$kappa,
    theta = vas_true$theta,
    sigma = vas_true$sigma,
    meas_sd = vas_true$meas_sd
  ),
  est = unlist(res_vas_sim$params)[c("kappa", "theta", "sigma", "meas_sd")]
))

Y_pred_vas_sim <- res_vas_sim$filtered$Y_pred
Y_filt_vas_sim <- res_vas_sim$filtered$Y_filt
resid_pred_vas_sim <- Y_sim_vas - Y_pred_vas_sim
resid_filt_vas_sim <- Y_sim_vas - Y_filt_vas_sim

cat("\nVASICEK SIMULATION: residual SD, one-step-ahead:\n")
print(setNames(apply(resid_pred_vas_sim, 2, sd), paste0(taus, "Y")))

cat("\nVASICEK SIMULATION: residual SD, filtered reconstruction:\n")
print(setNames(apply(resid_filt_vas_sim, 2, sd), paste0(taus, "Y")))

plot_swap_grid(
  seq_len(nrow(Y_sim_vas)), Y_sim_vas, Y_pred_vas_sim,
  "VASICEK SIMULATION: observed vs ONE-STEP-AHEAD predicted swaps",
  model_col = "red", model_label = "Predicted"
)

plot_swap_grid(
  seq_len(nrow(Y_sim_vas)), Y_sim_vas, Y_filt_vas_sim,
  "VASICEK SIMULATION: observed vs FILTERED reconstructed swaps",
  model_col = "red", model_label = "Filtered"
)

r_true_vas <- r_sim_vas[-1]
r_filt_vas_est <- res_vas_sim$filtered$r_filt
r_pred_vas_est <- res_vas_sim$filtered$r_pred

plot(
  r_true_vas,
  type = "l",
  xlab = "t",
  ylab = "r_t",
  main = "VASICEK simulation: true vs filtered/predicted r_t",
  lwd = 1.5,
  col = "black"
)
lines(r_filt_vas_est, lty = 2, lwd = 1.5, col = "red")
lines(r_pred_vas_est, lty = 3, lwd = 1.5, col = "blue")
legend(
  "topright",
  legend = c("True", "Filtered", "Predicted"),
  col = c("black", "red", "blue"),
  lty = c(1, 2, 3),
  lwd = 1.5,
  bty = "n"
)

plot(
  r_filt_vas_est - r_true_vas,
  type = "l",
  xlab = "t",
  ylab = "Filtered - true",
  main = "VASICEK simulation: filtering error",
  lwd = 1.5,
  col = "black"
)
abline(h = 0, lty = 2, col = "red", lwd = 1.2)

kf_vas_sim_true_params <- kalman_vasicek_swaps_EKF(
  Y = Y_sim_vas,
  taus = taus,
  dt = dt,
  kappa = vas_true$kappa,
  theta = vas_true$theta,
  sigma = vas_true$sigma,
  meas_sd = vas_true$meas_sd,
  pay_freq = pay_freq,
  return_states = TRUE
)

plot(
  r_true_vas,
  type = "l",
  lwd = 2,
  xlab = "t",
  ylab = expression(r[t]),
  main = "VASICEK simulation: true state vs EKF filtered states",
  col = "black"
)
lines(res_vas_sim$filtered$r_filt, lty = 2, lwd = 2, col = "red")
lines(kf_vas_sim_true_params$r_filt, lty = 3, lwd = 2, col = "blue")
legend(
  "topright",
  legend = c("True simulated state", "Filtered, estimated params", "Filtered, true params"),
  col = c("black", "red", "blue"),
  lty = c(1, 2, 3),
  lwd = 2,
  bty = "n"
)

err_vas_est_params <- res_vas_sim$filtered$r_filt - r_true_vas
err_vas_true_params <- kf_vas_sim_true_params$r_filt - r_true_vas

plot(
  err_vas_est_params,
  type = "l",
  lwd = 2,
  xlab = "t",
  ylab = "Filtered - true",
  main = "VASICEK filtering error: estimated params vs true params",
  col = "black"
)
lines(err_vas_true_params, lty = 2, lwd = 2, col = "red")
abline(h = 0, lty = 3)
legend(
  "topright",
  legend = c("Estimated parameters", "True parameters"),
  col = c("black", "red"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

rmse_pred_vas_sim <- sqrt(colMeans(resid_pred_vas_sim^2))
rmse_filt_vas_sim <- sqrt(colMeans(resid_filt_vas_sim^2))

plot_rmse_two_lines(
  taus, rmse_pred_vas_sim, rmse_filt_vas_sim,
  "VASICEK SIMULATION: RMSE by maturity",
  label_1 = "One-step-ahead", label_2 = "Filtered",
  col_1 = "black", col_2 = "red"
)

###############################################################
# 10) CIR SIMULATION STUDY
###############################################################

simulate_cir_exact <- function(r0, kappa, theta, sigma, dt, nSteps, seed = 1) {
  set.seed(seed)
  
  r <- numeric(nSteps + 1)
  r[1] <- max(r0, 1e-12)
  
  for (tt in 1:nSteps) {
    c_scale <- sigma^2 * (1 - exp(-kappa * dt)) / (4 * kappa)
    df <- 4 * kappa * theta / sigma^2
    ncp <- 4 * kappa * exp(-kappa * dt) * r[tt] /
      (sigma^2 * (1 - exp(-kappa * dt)))
    
    r[tt + 1] <- c_scale * rchisq(1, df = df, ncp = ncp)
  }
  
  r
}

simulate_swap_panel_cir <- function(r_path, taus, kappa, theta, sigma,
                                    meas_sd, pay_freq = 1, seed = 2) {
  set.seed(seed)
  
  N <- length(r_path) - 1
  m <- length(taus)
  Y <- matrix(NA_real_, nrow = N, ncol = m)
  colnames(Y) <- paste0("Y", taus)
  
  for (tt in 1:N) {
    r_t <- r_path[tt + 1]
    s_true <- sapply(taus, function(T) {
      swap_rate_from_r_cir(r_t, T, kappa, theta, sigma, pay_freq)
    })
    Y[tt, ] <- s_true + rnorm(m, 0, meas_sd)
  }
  
  Y
}

cir_true <- list(
  kappa = 0.30,
  theta = 0.03,
  sigma = 0.08,
  meas_sd = 0.001,
  r0 = 0.03
)

cat("\nCIR TRUE Feller check 2*kappa*theta - sigma^2 =",
    2 * cir_true$kappa * cir_true$theta - cir_true$sigma^2, "\n")

r_sim_cir <- simulate_cir_exact(
  r0 = cir_true$r0,
  kappa = cir_true$kappa,
  theta = cir_true$theta,
  sigma = cir_true$sigma,
  dt = dt,
  nSteps = N_sim,
  seed = 123
)

Y_sim_cir <- simulate_swap_panel_cir(
  r_path = r_sim_cir,
  taus = taus,
  kappa = cir_true$kappa,
  theta = cir_true$theta,
  sigma = cir_true$sigma,
  meas_sd = cir_true$meas_sd,
  pay_freq = pay_freq,
  seed = 999
)

res_cir_sim <- fit_cir_EKF_multistart(
  Y = Y_sim_cir,
  taus = taus,
  dt = dt,
  pay_freq = pay_freq,
  n_starts = n_starts,
  n_cores = n_cores
)

cat("\n=== CIR SIMULATION RECOVERY TEST ===\n")
print(cbind(
  true = c(
    kappa = cir_true$kappa,
    theta = cir_true$theta,
    sigma = cir_true$sigma,
    meas_sd = cir_true$meas_sd
  ),
  est = unlist(res_cir_sim$params)[c("kappa", "theta", "sigma", "meas_sd")]
))

cat("\nCIR estimated simulation Feller check:\n")
print(res_cir_sim$params$feller)

Y_pred_cir_sim <- res_cir_sim$filtered$Y_pred
Y_filt_cir_sim <- res_cir_sim$filtered$Y_filt
resid_pred_cir_sim <- Y_sim_cir - Y_pred_cir_sim
resid_filt_cir_sim <- Y_sim_cir - Y_filt_cir_sim

cat("\nCIR SIMULATION: residual SD, one-step-ahead:\n")
print(setNames(apply(resid_pred_cir_sim, 2, sd), paste0(taus, "Y")))

cat("\nCIR SIMULATION: residual SD, filtered reconstruction:\n")
print(setNames(apply(resid_filt_cir_sim, 2, sd), paste0(taus, "Y")))

plot_swap_grid(
  seq_len(nrow(Y_sim_cir)), Y_sim_cir, Y_pred_cir_sim,
  "CIR SIMULATION: observed vs ONE-STEP-AHEAD predicted swaps",
  model_col = "blue", model_label = "Predicted"
)

plot_swap_grid(
  seq_len(nrow(Y_sim_cir)), Y_sim_cir, Y_filt_cir_sim,
  "CIR SIMULATION: observed vs FILTERED reconstructed swaps",
  model_col = "blue", model_label = "Filtered"
)

r_true_cir <- r_sim_cir[-1]
r_filt_cir_est <- res_cir_sim$filtered$r_filt
r_pred_cir_est <- res_cir_sim$filtered$r_pred

plot(
  r_true_cir,
  type = "l",
  xlab = "t",
  ylab = "r_t",
  main = "CIR simulation: true vs filtered/predicted r_t",
  lwd = 1.5,
  col = "black"
)
lines(r_filt_cir_est, lty = 2, lwd = 1.5, col = "red")
lines(r_pred_cir_est, lty = 3, lwd = 1.5, col = "blue")
legend(
  "topright",
  legend = c("True", "Filtered", "Predicted"),
  col = c("black", "red", "blue"),
  lty = c(1, 2, 3),
  lwd = 1.5,
  bty = "n"
)

plot(
  r_filt_cir_est - r_true_cir,
  type = "l",
  xlab = "t",
  ylab = "Filtered - true",
  main = "CIR simulation: filtering error",
  lwd = 1.5,
  col = "black"
)
abline(h = 0, lty = 2, col = "red", lwd = 1.2)

kf_cir_sim_true_params <- kalman_cir_swaps_EKF(
  Y = Y_sim_cir,
  taus = taus,
  dt = dt,
  kappa = cir_true$kappa,
  theta = cir_true$theta,
  sigma = cir_true$sigma,
  meas_sd = cir_true$meas_sd,
  pay_freq = pay_freq,
  return_states = TRUE
)

plot(
  r_true_cir,
  type = "l",
  lwd = 2,
  xlab = "t",
  ylab = expression(r[t]),
  main = "CIR simulation: true state vs EKF filtered states",
  col = "black"
)
lines(res_cir_sim$filtered$r_filt, lty = 2, lwd = 2, col = "red")
lines(kf_cir_sim_true_params$r_filt, lty = 3, lwd = 2, col = "blue")
legend(
  "topright",
  legend = c("True simulated state", "Filtered, estimated params", "Filtered, true params"),
  col = c("black", "red", "blue"),
  lty = c(1, 2, 3),
  lwd = 2,
  bty = "n"
)

err_cir_est_params <- res_cir_sim$filtered$r_filt - r_true_cir
err_cir_true_params <- kf_cir_sim_true_params$r_filt - r_true_cir

plot(
  err_cir_est_params,
  type = "l",
  lwd = 2,
  xlab = "t",
  ylab = "Filtered - true",
  main = "CIR filtering error: estimated params vs true params",
  col = "black"
)
lines(err_cir_true_params, lty = 2, lwd = 2, col = "red")
abline(h = 0, lty = 3)
legend(
  "topright",
  legend = c("Estimated parameters", "True parameters"),
  col = c("black", "red"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

rmse_pred_cir_sim <- sqrt(colMeans(resid_pred_cir_sim^2))
rmse_filt_cir_sim <- sqrt(colMeans(resid_filt_cir_sim^2))

plot_rmse_two_lines(
  taus, rmse_pred_cir_sim, rmse_filt_cir_sim,
  "CIR SIMULATION: RMSE by maturity",
  label_1 = "One-step-ahead", label_2 = "Filtered",
  col_1 = "black", col_2 = "blue"
)

###############################################################
# 11) SIMULATION COMPARISON: VASICEK VS CIR
###############################################################

plot_rmse_two_lines(
  taus, rmse_pred_vas_sim, rmse_pred_cir_sim,
  "SIMULATION: one-step-ahead RMSE, Vasicek vs CIR",
  label_1 = "Vasicek simulation", label_2 = "CIR simulation",
  col_1 = "red", col_2 = "blue"
)

###############################################################
# 12) JACOBIAN SANITY CHECKS
###############################################################

r0 <- 0.03
T0 <- 10
h  <- 1e-6

num_vas <- (
  swap_rate_from_r_vas(r0 + h, T0, vas_true$kappa, vas_true$theta, vas_true$sigma, pay_freq) -
    swap_rate_from_r_vas(r0 - h, T0, vas_true$kappa, vas_true$theta, vas_true$sigma, pay_freq)
) / (2 * h)

ana_vas <- swap_rate_dr_vas(r0, T0, vas_true$kappa, vas_true$theta, vas_true$sigma, pay_freq)

cat("\n=== VASICEK Jacobian sanity check ===\n")
print(c(numerical = num_vas, analytical = ana_vas))

num_cir <- (
  swap_rate_from_r_cir(r0 + h, T0, cir_true$kappa, cir_true$theta, cir_true$sigma, pay_freq) -
    swap_rate_from_r_cir(r0 - h, T0, cir_true$kappa, cir_true$theta, cir_true$sigma, pay_freq)
) / (2 * h)

ana_cir <- swap_rate_dr_cir(r0, T0, cir_true$kappa, cir_true$theta, cir_true$sigma, pay_freq)

cat("\n=== CIR Jacobian sanity check ===\n")
print(c(numerical = num_cir, analytical = ana_cir))


###############################################################
# DONE
###############################################################


vas_params <- unlist(res_vas_real$params)
cir_params <- unlist(res_cir_real$params)

param_table <- rbind(
  Vasicek = vas_params[c("kappa", "theta", "sigma", "meas_sd")],
  CIR     = cir_params[c("kappa", "theta", "sigma", "meas_sd")]
)

print(param_table)



# Optional:
# save.image("~/Desktop/combined_vasicek_cir_ekf_simulation_workspace.RData")
# load("~/Desktop/combined_vasicek_cir_ekf_simulation_workspace.RData")
