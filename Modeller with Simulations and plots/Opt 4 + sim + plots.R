# ============================================================
# OPT 4
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(parallel)
library(doParallel)
library(foreach)
library(ggplot2)

dt_monthly <- 1 / 12

n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1)

n_starts     <- 50
n_starts_sim <- 50

set.seed(123)

cl <- parallel::makeCluster(n_cores, type = "PSOCK")
doParallel::registerDoParallel(cl)

on.exit({
  try(parallel::stopCluster(cl), silent = TRUE)
}, add = TRUE)

cat("Using", n_cores, "cores\n")
cat("Registered foreach workers:", foreach::getDoParWorkers(), "\n")


swap <- read_excel("EURIBOR swap data.xlsx")
ois  <- read_excel("OIS zero bond.xlsx")

names(swap)[1] <- "Date"
names(ois)[1]  <- "Date"

swap <- swap %>%
  mutate(Date = as.Date(Date)) %>%
  filter(!is.na(Date)) %>%
  arrange(Date)

ois <- ois %>%
  mutate(Date = as.Date(Date)) %>%
  filter(!is.na(Date)) %>%
  arrange(Date)

swap_maturities <- c(1, 2, 3, 5, 7, 10, 15, 20, 30)
swap_cols       <- c("Y1", "Y2", "Y3", "Y5", "Y7", "Y10", "Y15", "Y20", "Y30")

p_name <- function(T) paste0("P_T_", sprintf("%.1f", T))

# ============================================================
# MERGE DATA + BOOTSTRAP A-BLOCKS
# ============================================================

first_ois_date <- min(ois$Date, na.rm = TRUE)

data_merged <- full_join(
  swap %>% select(Date, all_of(swap_cols)),
  ois,
  by = "Date"
) %>%
  arrange(Date) %>%
  fill(starts_with("P_T_"), .direction = "down") %>%
  filter(Date >= first_ois_date) %>%
  filter(!is.na(Y1))

fixed_delta <- 0.5

coupon_grid <- function(n) seq(0.5, n, by = 0.5)

needed_ois_cols <- unique(unlist(lapply(
  swap_maturities,
  function(n) p_name(coupon_grid(n))
)))

missing_cols <- setdiff(needed_ois_cols, names(data_merged))

if (length(missing_cols) > 0) {
  stop("Missing OIS columns: ", paste(missing_cols, collapse = ", "))
}

compute_B_row <- function(row_df) {
  out <- numeric(length(swap_maturities))
  
  for (j in seq_along(swap_maturities)) {
    n <- swap_maturities[j]
    S <- as.numeric(row_df[[swap_cols[j]]]) / 100
    
    pay_dates <- coupon_grid(n)
    p_cols <- p_name(pay_dates)
    
    P_vals <- sapply(p_cols, function(col) as.numeric(row_df[[col]]))
    P_sum  <- sum(P_vals, na.rm = TRUE)
    Pn     <- as.numeric(row_df[[p_name(n)]])
    
    out[j] <- fixed_delta * S * P_sum - 1 + Pn
  }
  
  names(out) <- paste0("B_", swap_maturities, "Y")
  out
}

B_list <- lapply(seq_len(nrow(data_merged)), function(i) {
  compute_B_row(data_merged[i, , drop = FALSE])
})

B_mat <- do.call(rbind, B_list)
B_df  <- bind_cols(data_merged["Date"], as.data.frame(B_mat))

A_df <- B_df %>%
  mutate(
    A_0_1Y   = B_1Y,
    A_1_2Y   = B_2Y  - B_1Y,
    A_2_3Y   = B_3Y  - B_2Y,
    A_3_5Y   = B_5Y  - B_3Y,
    A_5_7Y   = B_7Y  - B_5Y,
    A_7_10Y  = B_10Y - B_7Y,
    A_10_15Y = B_15Y - B_10Y,
    A_15_20Y = B_20Y - B_15Y,
    A_20_30Y = B_30Y - B_20Y
  )

ois_meas_maturities <- c(0.5, 1, 2, 3, 5, 7, 10, 15, 20, 30)
ois_meas_cols <- p_name(ois_meas_maturities)

A_meas_cols <- c(
  "A_0_1Y", "A_1_2Y", "A_2_3Y", "A_3_5Y", "A_5_7Y",
  "A_7_10Y", "A_10_15Y", "A_15_20Y", "A_20_30Y"
)

model_df <- data_merged %>%
  select(Date, all_of(ois_meas_cols)) %>%
  left_join(
    A_df %>% select(Date, all_of(A_meas_cols)),
    by = "Date"
  ) %>%
  drop_na()

Z_obs <- as.matrix(model_df %>% select(-Date))
dates <- model_df$Date

colnames(Z_obs) <- c(ois_meas_cols, A_meas_cols)

cat("Observations:", nrow(model_df), "\n")
cat("Z_obs dimension:", dim(Z_obs)[1], "x", dim(Z_obs)[2], "\n")

# ============================================================
# HELPER 
# ============================================================

vasicek_mean <- function(x, tau, kappa, theta) {
  theta + exp(-kappa * tau) * (x - theta)
}

vasicek_var <- function(tau, kappa, sigma) {
  if (kappa < 1e-8) return(sigma^2 * tau)
  sigma^2 / (2 * kappa) * (1 - exp(-2 * kappa * tau))
}

vasicek_second_moment <- function(x, tau, kappa, theta, sigma) {
  m <- vasicek_mean(x, tau, kappa, theta)
  v <- vasicek_var(tau, kappa, sigma)
  v + m^2
}

state_transition <- function(x_prev, pars, dt = dt_monthly) {
  kx  <- pars["kappa_x"]
  thx <- pars["theta_x"]
  ky  <- pars["kappa_y"]
  thy <- pars["theta_y"]
  
  c(
    thx + exp(-kx * dt) * (x_prev[1] - thx),
    thy + exp(-ky * dt) * (x_prev[2] - thy)
  )
}

state_Q <- function(pars, dt = dt_monthly) {
  qx  <- vasicek_var(dt, pars["kappa_x"], pars["sigma_x"])
  qy  <- vasicek_var(dt, pars["kappa_y"], pars["sigma_y"])
  rho <- pars["rho_xy"]
  
  qxy <- rho * sqrt(qx * qy)
  
  matrix(
    c(qx,  qxy,
      qxy, qy),
    nrow = 2,
    ncol = 2,
    byrow = TRUE
  )
}

# ============================================================
# MEASUREMENT FUNCTION
# ============================================================

measurement_function <- function(state, pars) {
  x <- state[1]
  y <- state[2]
  
  kx  <- pars["kappa_x"]
  thx <- pars["theta_x"]
  sx  <- pars["sigma_x"]
  
  ky  <- pars["kappa_y"]
  thy <- pars["theta_y"]
  sy  <- pars["sigma_y"]
  
  alpha <- pars["alpha"]
  eta   <- pars["eta"]
  a     <- pars["a"]
  
  denom <- a + x + x^2
  
  if (!is.finite(denom) || abs(denom) < 1e-10) {
    return(rep(NA_real_, length(ois_meas_maturities) + length(A_meas_cols)))
  }
  
  P_model <- sapply(ois_meas_maturities, function(Tau) {
    EX  <- vasicek_mean(x, Tau, kx, thx)
    EX2 <- vasicek_second_moment(x, Tau, kx, thx, sx)
    
    exp(-alpha * Tau) * (a + EX + EX2) / denom
  })
  
  A_period <- function(T_start, delta) {
    EY  <- vasicek_mean(y, T_start, ky, thy)
    EY2 <- vasicek_second_moment(y, T_start, ky, thy, sy)
    
    delta * exp(-alpha * T_start) * (eta + EY + EY2) / denom
  }
  
  A_block_model <- function(start, end, delta) {
    Ts <- seq(start, end - delta, by = delta)
    sum(sapply(Ts, function(Ts_i) A_period(Ts_i, delta)))
  }
  
  A_model <- c(
    A_block_model(0, 1, 0.25),
    A_block_model(1, 2, 0.25),
    A_block_model(2, 3, 0.5),
    A_block_model(3, 5, 0.5),
    A_block_model(5, 7, 0.5),
    A_block_model(7, 10, 0.5),
    A_block_model(10, 15, 0.5),
    A_block_model(15, 20, 0.5),
    A_block_model(20, 30, 0.5)
  )
  
  c(P_model, A_model)
}

measurement_R <- function(pars, nP, nA) {
  diag(c(
    rep(pars["sigma_meas_P"]^2, nP),
    rep(pars["sigma_meas_A"]^2, nA)
  ))
}

num_jacobian <- function(f, x, eps = 1e-6) {
  n <- length(x)
  fx <- f(x)
  m <- length(fx)
  J <- matrix(0, m, n)
  
  for (i in 1:n) {
    xp <- x
    xm <- x
    
    xp[i] <- xp[i] + eps
    xm[i] <- xm[i] - eps
    
    J[, i] <- (f(xp) - f(xm)) / (2 * eps)
  }
  
  J
}

# ============================================================
# EKF 
# ============================================================

safe_chol <- function(S, jitter = 1e-10, max_tries = 8) {
  S <- (S + t(S)) / 2
  
  for (i in 0:max_tries) {
    out <- tryCatch(
      chol(S + diag(jitter * 10^i, nrow(S))),
      error = function(e) NULL
    )
    
    if (!is.null(out)) return(out)
  }
  
  stop("failed")
}

ekf_filter <- function(theta_vec, Z_obs) {
  theta_vec <- unname(theta_vec)
  
  pars <- c(
    kappa_x      = theta_vec[1],
    theta_x      = theta_vec[2],
    sigma_x      = theta_vec[3],
    kappa_y      = theta_vec[4],
    theta_y      = theta_vec[5],
    sigma_y      = theta_vec[6],
    alpha        = theta_vec[7],
    eta          = theta_vec[8],
    sigma_meas_P = theta_vec[9],
    sigma_meas_A = theta_vec[10],
    rho_xy       = theta_vec[11],
    a            = theta_vec[12]
  )
  
  n <- nrow(Z_obs)
  m <- ncol(Z_obs)
  
  nP <- length(ois_meas_maturities)
  nA <- length(A_meas_cols)
  
  x_pred <- matrix(NA_real_, n, 2)
  x_filt <- matrix(NA_real_, n, 2)
  
  P_pred <- array(NA_real_, c(2, 2, n))
  P_filt <- array(NA_real_, c(2, 2, n))
  
  Z_pred <- matrix(NA_real_, n, m)
  Z_filt <- matrix(NA_real_, n, m)
  innovations <- matrix(NA_real_, n, m)
  
  colnames(x_pred) <- c("X_pred", "Y_pred")
  colnames(x_filt) <- c("X_filt", "Y_filt")
  colnames(Z_pred) <- colnames(Z_obs)
  colnames(Z_filt) <- colnames(Z_obs)
  colnames(innovations) <- colnames(Z_obs)
  
  x_tt <- c(pars["theta_x"], pars["theta_y"])
  P_tt <- diag(c(0.01, 0.01), 2)
  
  loglik <- 0
  
  for (t in 1:n) {
    x_t_pred <- state_transition(x_tt, pars)
    
    Fmat <- matrix(
      c(
        exp(-pars["kappa_x"] * dt_monthly), 0,
        0, exp(-pars["kappa_y"] * dt_monthly)
      ),
      2, 2,
      byrow = TRUE
    )
    
    Qmat <- state_Q(pars)
    P_t_pred <- Fmat %*% P_tt %*% t(Fmat) + Qmat
    P_t_pred <- (P_t_pred + t(P_t_pred)) / 2
    
    hfun <- function(xx) measurement_function(xx, pars)
    
    z_t_pred <- hfun(x_t_pred)
    
    if (any(!is.finite(z_t_pred))) {
      stop("Non-finite measurement prediction")
    }
    
    Hmat <- num_jacobian(hfun, x_t_pred)
    
    if (any(!is.finite(Hmat))) {
      stop("Non-finite Jacobian")
    }
    
    Rmat <- measurement_R(pars, nP, nA)
    
    innov <- as.numeric(Z_obs[t, ] - z_t_pred)
    
    S_t <- Hmat %*% P_t_pred %*% t(Hmat) + Rmat
    S_t <- (S_t + t(S_t)) / 2
    
    U <- safe_chol(S_t)
    
    PHt <- P_t_pred %*% t(Hmat)
    Y_tmp <- forwardsolve(t(U), t(PHt), upper.tri = FALSE)
    X_tmp <- backsolve(U, Y_tmp, upper.tri = TRUE)
    K_t <- t(X_tmp)
    
    x_t_filt <- as.numeric(x_t_pred + K_t %*% innov)
    
    I2 <- diag(2)
    KH <- K_t %*% Hmat
    
    P_t_filt <- (I2 - KH) %*% P_t_pred %*% t(I2 - KH) +
      K_t %*% Rmat %*% t(K_t)
    
    P_t_filt <- (P_t_filt + t(P_t_filt)) / 2
    
    z_t_filt <- hfun(x_t_filt)
    
    innov_mat <- matrix(innov, ncol = 1)
    u <- forwardsolve(t(U), innov_mat, upper.tri = FALSE)
    mahal <- sum(u^2)
    logdet <- 2 * sum(log(diag(U)))
    
    ll_t <- -0.5 * (
      m * log(2 * pi) +
        logdet +
        mahal
    )
    
    loglik <- loglik + ll_t
    
    x_pred[t, ] <- x_t_pred
    x_filt[t, ] <- x_t_filt
    
    P_pred[, , t] <- P_t_pred
    P_filt[, , t] <- P_t_filt
    
    Z_pred[t, ] <- z_t_pred
    Z_filt[t, ] <- z_t_filt
    innovations[t, ] <- innov
    
    x_tt <- x_t_filt
    P_tt <- P_t_filt
  }
  
  list(
    loglik = as.numeric(loglik),
    x_pred = x_pred,
    x_filt = x_filt,
    P_pred = P_pred,
    P_filt = P_filt,
    Z_pred = Z_pred,
    Z_filt = Z_filt,
    innovations = innovations,
    pars = pars
  )
}

neg_loglik <- function(theta_vec, Z_obs) {
  out <- tryCatch(
    ekf_filter(theta_vec, Z_obs),
    error = function(e) NULL
  )
  
  if (is.null(out) || !is.finite(out$loglik)) return(1e12)
  -out$loglik
}

# ============================================================
# REAL DATA MLE OPT 4
# ============================================================

param_names <- c(
  "kappa_x", "theta_x", "sigma_x",
  "kappa_y", "theta_y", "sigma_y",
  "alpha", "eta",
  "sigma_meas_P", "sigma_meas_A",
  "rho_xy", "a"
)

lower <- c(
  0.001, -0.1, 1e-4,
  0.001, -0.05, 1e-4,
  0.001, -0.05,
  1e-6, 1e-6,
  -0.999,
  0.2501
)

upper <- c(
  5.00, 0.3, 0.50,
  5.00, 0.7, 0.50,
  0.15, 0.05,
  0.05, 0.05,
  0.999,
  5
)

names(lower) <- param_names
names(upper) <- param_names

random_start <- function() {
  out <- runif(length(lower), lower, upper)
  names(out) <- param_names
  out
}

starts <- replicate(n_starts, random_start(), simplify = FALSE)

export_objs <- c(
  "neg_loglik", "ekf_filter", "state_transition", "state_Q",
  "vasicek_mean", "vasicek_var", "vasicek_second_moment",
  "measurement_function", "measurement_R", "num_jacobian",
  "safe_chol", "dt_monthly", "ois_meas_maturities",
  "A_meas_cols", "ois_meas_cols", "lower", "upper",
  "param_names", "Z_obs"
)

fit_list <- foreach(
  st = starts,
  .packages = c("stats"),
  .export = export_objs,
  .errorhandling = "pass"
) %dopar% {
  optim(
    par = st,
    fn = neg_loglik,
    Z_obs = Z_obs,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 500)
  )
}

fit_list <- Filter(function(x) !inherits(x, "error"), fit_list)

if (length(fit_list) == 0) {
  stop("All real-data MLE optimizations failed.")
}

obj_vals <- sapply(fit_list, `[[`, "value")
best_idx <- which.min(obj_vals)
best_fit <- fit_list[[best_idx]]

theta_hat <- best_fit$par
ekf_out_real <- ekf_filter(theta_hat, Z_obs)

cat("\nEstimated parameters with rho_xy and a:\n")
print(ekf_out_real$pars)

cat("\nBest start:", best_idx, "out of", length(fit_list), "successful starts\n")
cat("Best objective value:", best_fit$value, "\n")

param_real <- data.frame(
  parameter = names(ekf_out_real$pars),
  estimate  = as.numeric(ekf_out_real$pars)
)

print(param_real)

cat("\nEstimated model correlation rho_xy:\n")
print(ekf_out_real$pars["rho_xy"])

cat("\nEstimated a:\n")
print(ekf_out_real$pars["a"])

cat("\nEmpirical correlation between filtered states X and Y:\n")
print(cor(ekf_out_real$x_filt[, 1], ekf_out_real$x_filt[, 2], use = "complete.obs"))

cat("\nEmpirical correlation between predicted states X and Y:\n")
print(cor(ekf_out_real$x_pred[, 1], ekf_out_real$x_pred[, 2], use = "complete.obs"))

all_param_real <- do.call(rbind, lapply(seq_along(fit_list), function(i) {
  c(start = i, objective = fit_list[[i]]$value, fit_list[[i]]$par)
}))
all_param_real <- as.data.frame(all_param_real)
all_param_real <- all_param_real[order(all_param_real$objective), ]

cat("\nAll successful starts sorted by objective value:\n")
print(all_param_real)

# ============================================================
# STATE PLOTS
# ============================================================

plot_df_real <- data.frame(
  Date = dates,
  X_filtered = ekf_out_real$x_filt[, 1],
  Y_filtered = ekf_out_real$x_filt[, 2]
)

print(
  ggplot(plot_df_real, aes(Date)) +
    geom_line(aes(y = X_filtered), color = "hotpink", linewidth = 1) +
    ggtitle("Filtered X state") +
    labs(y = NULL) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
)

print(
  ggplot(plot_df_real, aes(Date)) +
    geom_line(aes(y = Y_filtered), color = "hotpink", linewidth = 1) +
    ggtitle("Filtered Y state") +
    labs(y = NULL) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
)

# ============================================================
# SWAP RATE 
# ============================================================

model_swap_rates <- function(state, pars) {
  x <- state[1]
  y <- state[2]
  
  kx  <- pars["kappa_x"]
  thx <- pars["theta_x"]
  sx  <- pars["sigma_x"]
  
  ky  <- pars["kappa_y"]
  thy <- pars["theta_y"]
  sy  <- pars["sigma_y"]
  
  alpha <- pars["alpha"]
  eta   <- pars["eta"]
  a     <- pars["a"]
  
  denom <- a + x + x^2
  
  P_fun <- function(Tau) {
    EX  <- vasicek_mean(x, Tau, kx, thx)
    EX2 <- vasicek_second_moment(x, Tau, kx, thx, sx)
    
    exp(-alpha * Tau) * (a + EX + EX2) / denom
  }
  
  A_period <- function(T_start, delta) {
    EY  <- vasicek_mean(y, T_start, ky, thy)
    EY2 <- vasicek_second_moment(y, T_start, ky, thy, sy)
    
    delta * exp(-alpha * T_start) *
      (eta + EY + EY2) / denom
  }
  
  A_sum <- function(Tn) {
    if (Tn <= 2) {
      starts <- seq(0, Tn - 0.25, by = 0.25)
      deltas <- rep(0.25, length(starts))
    } else {
      starts_3m <- seq(0, 2 - 0.25, by = 0.25)
      starts_6m <- seq(2, Tn - 0.5, by = 0.5)
      
      starts <- c(starts_3m, starts_6m)
      deltas <- c(
        rep(0.25, length(starts_3m)),
        rep(0.5, length(starts_6m))
      )
    }
    
    sum(mapply(A_period, starts, deltas))
  }
  
  sapply(swap_maturities, function(Tn) {
    pay_dates <- seq(0.5, Tn, by = 0.5)
    P_vals <- sapply(pay_dates, P_fun)
    
    (1 - P_fun(Tn) + A_sum(Tn)) /
      (fixed_delta * sum(P_vals))
  })
}

# ============================================================
# OBSERVED VS PREDICTED SWAPS
# ============================================================

swap_obs_real <- model_df %>%
  select(Date) %>%
  left_join(
    data_merged %>% select(Date, all_of(swap_cols)),
    by = "Date"
  ) %>%
  select(all_of(swap_cols)) %>%
  as.matrix() / 100

swap_pred_real <- t(sapply(seq_len(nrow(ekf_out_real$x_pred)), function(i) {
  model_swap_rates(
    state = ekf_out_real$x_pred[i, ],
    pars  = ekf_out_real$pars
  )
}))

swap_filt_real <- t(sapply(seq_len(nrow(ekf_out_real$x_filt)), function(i) {
  model_swap_rates(
    state = ekf_out_real$x_filt[i, ],
    pars  = ekf_out_real$pars
  )
}))

colnames(swap_pred_real) <- swap_cols
colnames(swap_filt_real) <- swap_cols

resid_swap_pred_real <- swap_obs_real - swap_pred_real
resid_swap_filt_real <- swap_obs_real - swap_filt_real

cat("\nSwap plot dimensions:\n")
print(dim(swap_obs_real))
print(dim(swap_pred_real))

cat("\nREAL DATA: swap residual mean:\n")
print(setNames(colMeans(resid_swap_pred_real), swap_cols))

cat("\nREAL DATA: swap residual SD:\n")
print(setNames(apply(resid_swap_pred_real, 2, sd), swap_cols))

cat("\nREAL DATA: filtered reconstruction swap residual SD:\n")
print(setNames(apply(resid_swap_filt_real, 2, sd), swap_cols))

op <- par(no.readonly = TRUE)
par(mfrow = c(3, 3), mar = c(3, 3, 3, 1), oma = c(0, 0, 3, 0))

idx <- 3:length(dates)

for (j in seq_along(swap_cols)) {
  plot(
    dates[idx],
    swap_obs_real[idx, j],
    type = "l",
    lwd = 1,
    xlab = "",
    ylab = "",
    main = paste0(swap_maturities[j], "Y"),
    ylim = c(-0.015, 0.05)
  )
  
  lines(
    dates[idx],
    swap_pred_real[idx, j],
    lty = 2,
    lwd = 1,
    col = "red"
  )
}

mtext("REAL DATA: observed vs predicted EURIBOR swap rates", outer = TRUE, cex = 1.2)

legend(
  "bottom",
  legend = c("Observed", "Predicted"),
  lty = c(1, 2),
  lwd = c(1, 1),
  horiz = TRUE,
  bty = "n"
)

par(op)

# ============================================================
# PREDICTED RESIDUAL PLOTS
# ============================================================

op <- par(no.readonly = TRUE)
par(mfrow = c(3, 3), mar = c(3, 3, 3, 1), oma = c(0, 0, 3, 0))

for (j in seq_along(swap_cols)) {
  plot(
    dates[idx],
    resid_swap_pred_real[idx, j],
    type = "l",
    lwd = 1,
    xlab = "",
    ylab = "",
    main = paste0(" ", swap_maturities[j], "Y"),
    col = "red",
    ylim = c(-0.03, 0.015)
  )
  
  abline(h = 0, lty = 2)
}

mtext("REAL DATA: swap residuals", outer = TRUE, cex = 1.2)
par(op)

# ============================================================
# SIMULATION
# ============================================================

rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))

simulate_current_model <- function(Tn, true_pars, seed = 123) {
  set.seed(seed)
  
  X <- numeric(Tn)
  Y <- numeric(Tn)
  
  X[1] <- true_pars["theta_x"]
  Y[1] <- true_pars["theta_y"]
  
  phi_x <- exp(-true_pars["kappa_x"] * dt_monthly)
  phi_y <- exp(-true_pars["kappa_y"] * dt_monthly)
  
  Q <- state_Q(true_pars)
  L <- t(chol(Q))
  
  for (t in 2:Tn) {
    eps <- as.numeric(L %*% rnorm(2))
    
    X[t] <- true_pars["theta_x"] +
      phi_x * (X[t - 1] - true_pars["theta_x"]) +
      eps[1]
    
    Y[t] <- true_pars["theta_y"] +
      phi_y * (Y[t - 1] - true_pars["theta_y"]) +
      eps[2]
  }
  
  m <- length(ois_meas_cols) + length(A_meas_cols)
  
  Z_true <- matrix(NA_real_, Tn, m)
  Z_obs_sim <- matrix(NA_real_, Tn, m)
  
  colnames(Z_true) <- colnames(Z_obs)
  colnames(Z_obs_sim) <- colnames(Z_obs)
  
  for (t in 1:Tn) {
    Z_true[t, ] <- measurement_function(c(X[t], Y[t]), true_pars)
    
    noise <- c(
      rnorm(length(ois_meas_cols), 0, true_pars["sigma_meas_P"]),
      rnorm(length(A_meas_cols),  0, true_pars["sigma_meas_A"])
    )
    
    Z_obs_sim[t, ] <- Z_true[t, ] + noise
  }
  
  list(
    X_true = X,
    Y_true = Y,
    Z_true = Z_true,
    Z_obs = Z_obs_sim
  )
}

estimate_current_model_once <- function(Z_sim, n_starts_sim = 10, seed = 123) {
  set.seed(seed)
  
  smart_start <- function() {
    c(
      kappa_x      = runif(1, 0.05, 0.70),
      theta_x      = runif(1, -0.01, 0.06),
      sigma_x      = runif(1, 0.002, 0.05),
      kappa_y      = runif(1, 0.05, 0.90),
      theta_y      = runif(1, -0.005, 0.02),
      sigma_y      = runif(1, 0.002, 0.04),
      alpha        = runif(1, 0.01, 0.10),
      eta          = runif(1, -0.005, 0.01),
      sigma_meas_P = runif(1, 0.0001, 0.005),
      sigma_meas_A = runif(1, 0.0001, 0.005),
      rho_xy       = runif(1, -0.80, 0.80),
      a            = runif(1, 0.3, 3)
    )
  }
  
  p0 <- c(
    kappa_x      = 0.25,
    theta_x      = 0.03,
    sigma_x      = 0.02,
    kappa_y      = 0.35,
    theta_y      = 0.002,
    sigma_y      = 0.01,
    alpha        = 0.05,
    eta          = 0.001,
    sigma_meas_P = 0.001,
    sigma_meas_A = 0.001,
    rho_xy       = 0.50,
    a            = 1
  )
  
  starts <- c(
    list(p0),
    replicate(n_starts_sim - 1, smart_start(), simplify = FALSE)
  )
  
  export_objs_sim <- c(
    "neg_loglik", "ekf_filter", "state_transition", "state_Q",
    "vasicek_mean", "vasicek_var", "vasicek_second_moment",
    "measurement_function", "measurement_R", "num_jacobian",
    "safe_chol", "dt_monthly", "ois_meas_maturities",
    "A_meas_cols", "ois_meas_cols", "lower", "upper",
    "param_names"
  )
  
  fits <- foreach(
    st = starts,
    .packages = c("stats"),
    .export = export_objs_sim,
    .errorhandling = "pass"
  ) %dopar% {
    optim(
      par = st,
      fn = neg_loglik,
      Z_obs = Z_sim,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = 400)
    )
  }
  
  fits <- Filter(function(x) !inherits(x, "error"), fits)
  
  if (length(fits) == 0) {
    stop("All simulated-data MLE optimizations failed.")
  }
  
  vals <- sapply(fits, `[[`, "value")
  best <- fits[[which.min(vals)]]
  ekf <- ekf_filter(best$par, Z_sim)
  
  list(
    theta_hat = best$par,
    ekf = ekf,
    opt_value = best$value,
    all_values = vals
  )
}

run_simulation_study_current_model <- function(
    Tn = nrow(Z_obs),
    n_starts_sim = 10,
    seed = 123
) {
  true_pars <- c(
    kappa_x      = 0.25,
    theta_x      = 0.03,
    sigma_x      = 0.02,
    kappa_y      = 0.35,
    theta_y      = 0.002,
    sigma_y      = 0.01,
    alpha        = 0.05,
    eta          = 0.001,
    sigma_meas_P = 0.001,
    sigma_meas_A = 0.001,
    rho_xy       = 0.50,
    a            = 1
  )
  
  cat("\n=SIMULATE DATA =\n")
  cat("True rho_xy:", true_pars["rho_xy"], "\n")
  cat("True a:", true_pars["a"], "\n")
  
  sim <- simulate_current_model(
    Tn = Tn,
    true_pars = true_pars,
    seed = seed
  )
  
  cat("Empirical correlation of simulated X and Y states:\n")
  print(cor(sim$X_true, sim$Y_true))
  
  cat("\n= FILTER WITH TRUE PARAMETERS =\n")
  
  ekf_true <- ekf_filter(true_pars, sim$Z_obs)
  
  cat("RMSE X, true-param filter:", rmse(sim$X_true, ekf_true$x_filt[, 1]), "\n")
  cat("RMSE Y, true-param filter:", rmse(sim$Y_true, ekf_true$x_filt[, 2]), "\n")
  
  cat("\nObservation RMSE, true parameters, Z_pred vs Z_true:\n")
  print(data.frame(
    series = colnames(Z_obs),
    rmse = sqrt(colMeans((sim$Z_true - ekf_true$Z_pred)^2))
  ))
  
  cat("\n= ESTIMATE PARAMETERS ON SIMULATED DATA =\n")
  cat("Number of starts:", n_starts_sim, "\n")
  
  fit_sim <- estimate_current_model_once(
    Z_sim = sim$Z_obs,
    n_starts_sim = n_starts_sim,
    seed = seed
  )
  
  true_vs_est <- data.frame(
    parameter = names(true_pars),
    true = as.numeric(true_pars),
    estimated = as.numeric(fit_sim$ekf$pars)
  )
  
  cat("\n=TRUE VS ESTIMATED PARAMETERS =\n")
  print(true_vs_est)
  
  cat("\n=STATE RMSE, ESTIMATED PARAMETERS =\n")
  cat("RMSE X:", rmse(sim$X_true, fit_sim$ekf$x_filt[, 1]), "\n")
  cat("RMSE Y:", rmse(sim$Y_true, fit_sim$ekf$x_filt[, 2]), "\n")
  
  cat("\nObservation RMSE, estimated parameters, Z_pred vs Z_true:\n")
  print(data.frame(
    series = colnames(Z_obs),
    rmse = sqrt(colMeans((sim$Z_true - fit_sim$ekf$Z_pred)^2))
  ))
  
  list(
    true_pars = true_pars,
    sim = sim,
    ekf_true = ekf_true,
    fit_sim = fit_sim,
    true_vs_est = true_vs_est
  )
}

res_simstudy <- run_simulation_study_current_model(
  Tn = nrow(Z_obs),
  n_starts_sim = n_starts_sim,
  seed = 123
)

# ============================================================
# SIMULATION STATE PLOTS
# ============================================================

plot_df_sim <- data.frame(
  t = seq_along(res_simstudy$sim$X_true),
  X_true = res_simstudy$sim$X_true,
  X_pred_truepar = res_simstudy$ekf_true$x_pred[, 1],
  X_filt_truepar = res_simstudy$ekf_true$x_filt[, 1],
  X_pred_estpar = res_simstudy$fit_sim$ekf$x_pred[, 1],
  X_filt_estpar = res_simstudy$fit_sim$ekf$x_filt[, 1],
  Y_true = res_simstudy$sim$Y_true,
  Y_pred_truepar = res_simstudy$ekf_true$x_pred[, 2],
  Y_filt_truepar = res_simstudy$ekf_true$x_filt[, 2],
  Y_pred_estpar = res_simstudy$fit_sim$ekf$x_pred[, 2],
  Y_filt_estpar = res_simstudy$fit_sim$ekf$x_filt[, 2]
)

print(
  ggplot(plot_df_sim, aes(t)) +
    geom_line(aes(y = X_true, linetype = "True"), color = "hotpink", linewidth = 1) +
    geom_line(aes(y = X_pred_estpar, linetype = "Predicted")) +
    ggtitle("Simulation: X state") +
    labs(y = NULL, linetype = NULL) +
    theme_light()
)

print(
  ggplot(plot_df_sim, aes(t)) +
    geom_line(aes(y = Y_true, linetype = "True"), color = "hotpink", linewidth = 1) +
    geom_line(aes(y = Y_pred_estpar, linetype = "Predicted")) +
    ggtitle("Simulation: Y state") +
    labs(y = NULL, linetype = NULL) +
    theme_light()
)

# ============================================================
# FILTERING ERROR PLOTS
# ============================================================

plot_df_sim$error_X_pred_estpar <- plot_df_sim$X_pred_estpar - plot_df_sim$X_true
plot_df_sim$error_X_filt_estpar <- plot_df_sim$X_filt_estpar - plot_df_sim$X_true
plot_df_sim$error_Y_pred_estpar <- plot_df_sim$Y_pred_estpar - plot_df_sim$Y_true
plot_df_sim$error_Y_filt_estpar <- plot_df_sim$Y_filt_estpar - plot_df_sim$Y_true

print(
  ggplot(plot_df_sim, aes(t)) +
    geom_line(aes(y = error_X_pred_estpar, linetype = "Predicted - true")) +
    geom_line(aes(y = error_X_filt_estpar, linetype = "Filtered - true")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ggtitle("Simulation: X_t prediction/filtering error") +
    labs(y = "Error", linetype = NULL) +
    theme_light()
)

print(
  ggplot(plot_df_sim, aes(t)) +
    geom_line(aes(y = error_Y_pred_estpar, linetype = "Predicted - true")) +
    geom_line(aes(y = error_Y_filt_estpar, linetype = "Filtered - true")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ggtitle("Simulation: Y_t prediction/filtering error") +
    labs(y = "Error", linetype = NULL) +
    theme_light()
)

# ============================================================
# SIMULATED DATA: TRUE VS PREDICTED SWAP RATES
# ============================================================

swap_true_sim <- t(sapply(seq_len(length(res_simstudy$sim$X_true)), function(i) {
  model_swap_rates(
    state = c(res_simstudy$sim$X_true[i], res_simstudy$sim$Y_true[i]),
    pars  = res_simstudy$true_pars
  )
}))

swap_pred_sim <- t(sapply(seq_len(nrow(res_simstudy$fit_sim$ekf$x_pred)), function(i) {
  model_swap_rates(
    state = res_simstudy$fit_sim$ekf$x_pred[i, ],
    pars  = res_simstudy$fit_sim$ekf$pars
  )
}))

swap_filt_sim <- t(sapply(seq_len(nrow(res_simstudy$fit_sim$ekf$x_filt)), function(i) {
  model_swap_rates(
    state = res_simstudy$fit_sim$ekf$x_filt[i, ],
    pars  = res_simstudy$fit_sim$ekf$pars
  )
}))

colnames(swap_true_sim) <- swap_cols
colnames(swap_pred_sim) <- swap_cols
colnames(swap_filt_sim) <- swap_cols

resid_swap_pred_sim <- swap_true_sim - swap_pred_sim
resid_swap_filt_sim <- swap_true_sim - swap_filt_sim

cat("\nSIMULATED DATA: swap residual mean:\n")
print(setNames(colMeans(resid_swap_pred_sim), swap_cols))

cat("\nSIMULATED DATA: swap residual SD:\n")
print(setNames(apply(resid_swap_pred_sim, 2, sd), swap_cols))

cat("\nSIMULATED DATA: filtered reconstruction swap residual SD:\n")
print(setNames(apply(resid_swap_filt_sim, 2, sd), swap_cols))

op <- par(no.readonly = TRUE)
par(mfrow = c(3, 3), mar = c(3, 3, 3, 1), oma = c(0, 0, 3, 0))

for (j in seq_along(swap_cols)) {
  plot(
    seq_len(nrow(swap_true_sim)),
    swap_true_sim[, j],
    type = "l",
    lwd = 1,
    xlab = "",
    ylab = "",
    main = paste0(swap_maturities[j], "Y")
  )
  
  lines(
    seq_len(nrow(swap_pred_sim)),
    swap_pred_sim[, j],
    lty = 2,
    lwd = 1,
    col = "hotpink"
  )
}

mtext("SIMULATED DATA: true vs predicted EURIBOR swap rates", outer = TRUE, cex = 1.2)

legend(
  "bottom",
  legend = c("True simulated", "Predicted"),
  lty = c(1, 2),
  lwd = c(1, 1),
  horiz = TRUE,
  bty = "n"
)

par(op)

# ============================================================
# SIMULATED DATA: RESIDUAL PLOTS
# ============================================================

op <- par(no.readonly = TRUE)
par(mfrow = c(3, 3), mar = c(3, 3, 3, 1), oma = c(0, 0, 3, 0))

for (j in seq_along(swap_cols)) {
  plot(
    seq_len(nrow(resid_swap_pred_sim)),
    resid_swap_pred_sim[, j],
    type = "l",
    lwd = 1,
    xlab = "",
    ylab = "",
    main = paste0("Residual ", swap_maturities[j], "Y")
  )
  
  abline(h = 0, lty = 2)
}

mtext("SIMULATED DATA: swap residuals", outer = TRUE, cex = 1.2)

par(op)

# ============================================================
# RMSE
# ============================================================

rmse_idx <- seq(3, nrow(swap_obs_real))

rmse_pred <- sqrt(colMeans(resid_swap_pred_real[rmse_idx, ]^2))
rmse_filt <- sqrt(colMeans(resid_swap_filt_real[rmse_idx, ]^2))

cat("\nRMSE:\n")
print(setNames(rmse_pred, swap_cols))

cat("\nRMSE filtered:\n")
print(setNames(rmse_filt, swap_cols))

par(mfrow = c(1, 1))
plot(
  swap_maturities,
  rmse_pred,
  type = "b",
  pch = 19,
  lwd = 2,
  col = "hotpink",
  lty = 2,
  xlab = "Maturity (Years)",
  ylab = "RMSE",
  ylim = range(rmse_pred, rmse_filt)
)

lines(
  swap_maturities,
  rmse_filt,
  type = "b",
  pch = 19,
  lwd = 2,
  col = "blue",
  lty = 1
)

legend(
  "topright",
  legend = c("Predicted (EKF)", "Filtered (EKF)"),
  col = c("hotpink", "blue"),
  lty = c(2, 1),
  pch = 19,
  lwd = 2,
  bty = "n"
)






###############################################################
# INITIAL MODEL: SIMULATION TRUE VS PREDICTED
###############################################################

idx_sim <- 1:nrow(swap_true_sim)

cols <- c(
  "black",
  "#ff69b4"
)

op <- par(no.readonly = TRUE)

par(
  mfrow = c(3, 3),
  mar = c(1.8, 1.8, 1.2, 0.2),
  oma = c(7, 3.5, 1, 0),
  mgp = c(1.4, 0.4, 0),
  tcl = -0.25
)

for (j in seq_along(swap_maturities)) {
  
  y_min <- min(
    swap_true_sim[idx_sim, j],
    swap_pred_sim[idx_sim, j],
    na.rm = TRUE
  )
  
  y_max <- max(
    swap_true_sim[idx_sim, j],
    swap_pred_sim[idx_sim, j],
    na.rm = TRUE
  )
  
  plot(
    idx_sim,
    swap_true_sim[idx_sim, j],
    type = "l",
    main = paste0(swap_maturities[j], "Y"),
    xlab = "",
    ylab = "",
    lwd = 1.8,
    col = cols[1],
    ylim = c(y_min, y_max),
    cex.main = 1.45,
    cex.axis = 1.4
  )
  
  abline(h = axTicks(2), col = "grey85", lwd = 0.8)
  
  lines(idx_sim, swap_true_sim[idx_sim, j], lwd = 1.8, col = cols[1])
  lines(idx_sim, swap_pred_sim[idx_sim, j], lwd = 1.6, col = cols[2])
}

mtext("Time", side = 1, outer = TRUE, line = 1.1, cex = 1.2)
mtext("Swap rate", side = 2, outer = TRUE, line = 0.8, cex = 1.2)

par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE)
plot.new()

legend(
  x = "bottom",
  inset = 0.02,
  horiz = TRUE,
  bty = "n",
  cex = 1.25,
  xpd = NA,
  legend = c("Simulated observed", "Optimization 4"),
  col = cols,
  lty = 1,
  lwd = c(1.8, 1.6)
)

par(op)


# ============================================================
# DONE
# ============================================================

# Optional Windows-friendly save path:
# save.image(file.path(getwd(), "ekf_workspace_tilfoej_a_uden_udevidet.RData"))
# load(file.path(getwd(), "ekf_workspace_tilfoej_a_uden_udevidet.RData"))
