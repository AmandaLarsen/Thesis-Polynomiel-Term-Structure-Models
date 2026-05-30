# ============================================================
# INITIAL MODEL AND OPT. 1
# EKF, 5-year rolling window only, h = 6 months
# ============================================================

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

dt_monthly <- 1 / 12
fixed_delta <- 0.5


h <- 6                       # 6-month-ahead forecast
rolling_window <- 60         # 5 years monthly data


n_starts_oos <- 10
optim_maxit_oos <- 250


swap_file <- "~/Desktop/ægte data/EURIBOR swap data.xlsx"
ois_file  <- "~/Desktop/ægte data/OIS zero bond.xlsx"


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
ois  <- read_excel(ois_file)

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
coupon_grid <- function(n) seq(0.5, n, by = 0.5)

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

needed_ois_cols <- unique(unlist(lapply(swap_maturities, function(n) p_name(coupon_grid(n)))))
missing_cols <- setdiff(needed_ois_cols, names(data_merged))
if (length(missing_cols) > 0) stop("Missing OIS columns: ", paste(missing_cols, collapse = ", "))

compute_B_matrix <- function(data_merged) {
  out <- matrix(NA_real_, nrow(data_merged), length(swap_maturities))
  colnames(out) <- paste0("B_", swap_maturities, "Y")
  for (j in seq_along(swap_maturities)) {
    n <- swap_maturities[j]
    S <- as.numeric(data_merged[[swap_cols[j]]]) / 100
    pay_dates <- coupon_grid(n)
    p_cols <- p_name(pay_dates)
    P_sum <- rowSums(as.matrix(data_merged[, p_cols]), na.rm = TRUE)
    Pn <- as.numeric(data_merged[[p_name(n)]])
    out[, j] <- fixed_delta * S * P_sum - 1 + Pn
  }
  out
}

B_mat <- compute_B_matrix(data_merged)
B_df <- bind_cols(data_merged["Date"], as.data.frame(B_mat))

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
  left_join(A_df %>% select(Date, all_of(A_meas_cols)), by = "Date") %>%
  drop_na()

Z_obs <- as.matrix(model_df %>% select(-Date))
dates <- model_df$Date
colnames(Z_obs) <- c(ois_meas_cols, A_meas_cols)

swap_obs_real <- model_df %>%
  select(Date) %>%
  left_join(data_merged %>% select(Date, all_of(swap_cols)), by = "Date") %>%
  select(all_of(swap_cols)) %>%
  as.matrix() / 100

cat("Observations:", nrow(Z_obs), "\n")
cat("Z_obs dimension:", dim(Z_obs)[1], "x", dim(Z_obs)[2], "\n")

A_blocks <- list(
  list(start = 0,  end = 1,  delta = 0.25),
  list(start = 1,  end = 2,  delta = 0.25),
  list(start = 2,  end = 3,  delta = 0.5),
  list(start = 3,  end = 5,  delta = 0.5),
  list(start = 5,  end = 7,  delta = 0.5),
  list(start = 7,  end = 10, delta = 0.5),
  list(start = 10, end = 15, delta = 0.5),
  list(start = 15, end = 20, delta = 0.5),
  list(start = 20, end = 30, delta = 0.5)
)
A_block_grids <- lapply(A_blocks, function(b) seq(b$start, b$end - b$delta, by = b$delta))
A_block_deltas <- sapply(A_blocks, function(b) b$delta)

param_names_base <- c(
  "kappa_x", "theta_x", "sigma_x",
  "kappa_y", "theta_y", "sigma_y",
  "alpha", "eta",
  "sigma_meas_P", "sigma_meas_A"
)

param_names_opt1 <- c(param_names_base, "rho_xy")

bounds_initial <- list(
  lower = c(0.001, -0.050, 1e-4, 0.001, -0.050, 1e-4, 0.001, -0.050, 1e-6, 1e-6),
  upper = c(5.000,  0.150, 0.500, 5.000,  0.150, 0.500, 0.150,  0.050, 0.050, 0.050)
)
names(bounds_initial$lower) <- param_names_base
names(bounds_initial$upper) <- param_names_base

bounds_opt1 <- list(
  lower = c(0.001, -0.050, 1e-4, 0.001, -0.050, 1e-4, 0.001, -0.050, 1e-6, 1e-6, -0.999),
  upper = c(5.000,  0.150, 0.500, 5.000,  0.150, 0.500, 0.150,  0.050, 0.050, 0.050,  0.999)
)
names(bounds_opt1$lower) <- param_names_opt1
names(bounds_opt1$upper) <- param_names_opt1

start_initial <- c(
  kappa_x      = 0.0210,
  theta_x      = 0.1500,
  sigma_x      = 0.1837,
  kappa_y      = 0.1278,
  theta_y      = -0.0073,
  sigma_y      = 0.0032,
  alpha        = 0.0246,
  eta          = 0.0074,
  sigma_meas_P = 0.0373,
  sigma_meas_A = 0.0237
)

start_opt1 <- c(
  kappa_x      = 0.0216,
  theta_x      = 0.1500,
  sigma_x      = 0.1796,
  kappa_y      = 4.9366,
  theta_y      = -0.0127,
  sigma_y      = 0.0103,
  alpha        = 0.0241,
  eta          = 0.0146,
  sigma_meas_P = 0.0376,
  sigma_meas_A = 0.0244,
  rho_xy       = 0.9767
)

# ============================================================
# MODEL FUNCTIONS
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

expand_pars <- function(theta_vec, model_type = c("initial", "opt1")) {
  model_type <- match.arg(model_type)
  theta_vec <- unname(theta_vec)

  if (model_type == "initial") {
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
      rho_xy       = 0
    )
  } else {
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
      rho_xy       = theta_vec[11]
    )
  }

  pars
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
  matrix(c(qx, qxy, qxy, qy), nrow = 2, byrow = TRUE)
}

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

  # Initial/Opt.1 specification: no a parameter.
  denom <- x + x^2
  if (!is.finite(denom) || abs(denom) < 1e-10) {
    return(rep(NA_real_, length(ois_meas_maturities) + length(A_meas_cols)))
  }

  EX  <- vasicek_mean(x, ois_meas_maturities, kx, thx)
  EX2 <- vasicek_second_moment(x, ois_meas_maturities, kx, thx, sx)
  P_model <- exp(-alpha * ois_meas_maturities) * (EX + EX2) / denom

  A_period <- function(T_start, delta) {
    EY  <- vasicek_mean(y, T_start, ky, thy)
    EY2 <- vasicek_second_moment(y, T_start, ky, thy, sy)
    delta * exp(-alpha * T_start) * (eta + EY + EY2) / denom
  }

  A_model <- numeric(length(A_block_grids))
  for (i in seq_along(A_block_grids)) {
    Ts <- A_block_grids[[i]]
    delta_i <- A_block_deltas[i]
    A_model[i] <- sum(A_period(Ts, delta_i))
  }

  c(P_model, A_model)
}

measurement_R <- function(pars, nP, nA) {
  diag(c(rep(pars["sigma_meas_P"]^2, nP), rep(pars["sigma_meas_A"]^2, nA)))
}

num_jacobian <- function(f, x, eps = 1e-6) {
  n <- length(x)
  fx <- f(x)
  m <- length(fx)
  J <- matrix(0, m, n)
  for (i in seq_len(n)) {
    xp <- x; xm <- x
    xp[i] <- xp[i] + eps
    xm[i] <- xm[i] - eps
    J[, i] <- (f(xp) - f(xm)) / (2 * eps)
  }
  J
}

safe_chol <- function(S, jitter = 1e-10, max_tries = 8) {
  S <- (S + t(S)) / 2
  for (i in 0:max_tries) {
    out <- tryCatch(chol(S + diag(jitter * 10^i, nrow(S))), error = function(e) NULL)
    if (!is.null(out)) return(out)
  }
  stop("failed")
}

ekf_filter <- function(theta_vec, Z_obs, model_type = c("initial", "opt1")) {
  model_type <- match.arg(model_type)
  pars <- expand_pars(theta_vec, model_type)

  n <- nrow(Z_obs)
  m <- ncol(Z_obs)
  nP <- length(ois_meas_maturities)
  nA <- length(A_meas_cols)

  x_pred <- matrix(NA_real_, n, 2)
  x_filt <- matrix(NA_real_, n, 2)
  Z_pred <- matrix(NA_real_, n, m)

  colnames(x_pred) <- c("X_pred", "Y_pred")
  colnames(x_filt) <- c("X_filt", "Y_filt")
  colnames(Z_pred) <- colnames(Z_obs)

  x_tt <- c(pars["theta_x"], pars["theta_y"])
  P_tt <- diag(c(0.01, 0.01), 2)
  loglik <- 0

  for (t in seq_len(n)) {
    x_t_pred <- state_transition(x_tt, pars)
    Fmat <- matrix(c(exp(-pars["kappa_x"] * dt_monthly), 0, 0, exp(-pars["kappa_y"] * dt_monthly)), 2, 2, byrow = TRUE)
    Qmat <- state_Q(pars)
    P_t_pred <- Fmat %*% P_tt %*% t(Fmat) + Qmat
    P_t_pred <- (P_t_pred + t(P_t_pred)) / 2

    hfun <- function(xx) measurement_function(xx, pars)
    z_t_pred <- hfun(x_t_pred)
    if (any(!is.finite(z_t_pred))) stop("Non-finite measurement prediction")

    Hmat <- num_jacobian(hfun, x_t_pred)
    if (any(!is.finite(Hmat))) stop("Non-finite Jacobian")

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
    P_t_filt <- (I2 - KH) %*% P_t_pred %*% t(I2 - KH) + K_t %*% Rmat %*% t(K_t)
    P_t_filt <- (P_t_filt + t(P_t_filt)) / 2

    innov_mat <- matrix(innov, ncol = 1)
    u <- forwardsolve(t(U), innov_mat, upper.tri = FALSE)
    mahal <- sum(u^2)
    logdet <- 2 * sum(log(diag(U)))
    loglik <- loglik - 0.5 * (m * log(2 * pi) + logdet + mahal)

    x_pred[t, ] <- x_t_pred
    x_filt[t, ] <- x_t_filt
    Z_pred[t, ] <- z_t_pred

    x_tt <- x_t_filt
    P_tt <- P_t_filt
  }

  list(loglik = as.numeric(loglik), x_pred = x_pred, x_filt = x_filt, Z_pred = Z_pred, pars = pars)
}

neg_loglik <- function(theta_vec, Z_obs, model_type) {
  out <- tryCatch(ekf_filter(theta_vec, Z_obs, model_type), error = function(e) NULL)
  if (is.null(out) || !is.finite(out$loglik)) return(1e12)
  -out$loglik
}

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

  denom <- x + x^2
  if (!is.finite(denom) || abs(denom) < 1e-10) return(rep(NA_real_, length(swap_maturities)))

  P_fun <- function(Tau) {
    EX  <- vasicek_mean(x, Tau, kx, thx)
    EX2 <- vasicek_second_moment(x, Tau, kx, thx, sx)
    exp(-alpha * Tau) * (EX + EX2) / denom
  }

  A_period <- function(T_start, delta) {
    EY  <- vasicek_mean(y, T_start, ky, thy)
    EY2 <- vasicek_second_moment(y, T_start, ky, thy, sy)
    delta * exp(-alpha * T_start) * (eta + EY + EY2) / denom
  }

  A_sum <- function(Tn) {
    if (Tn <= 2) {
      starts <- seq(0, Tn - 0.25, by = 0.25)
      deltas <- rep(0.25, length(starts))
    } else {
      starts_3m <- seq(0, 2 - 0.25, by = 0.25)
      starts_6m <- seq(2, Tn - 0.5, by = 0.5)
      starts <- c(starts_3m, starts_6m)
      deltas <- c(rep(0.25, length(starts_3m)), rep(0.5, length(starts_6m)))
    }
    sum(mapply(A_period, starts, deltas))
  }

  sapply(swap_maturities, function(Tn) {
    pay_dates <- seq(0.5, Tn, by = 0.5)
    P_vals <- sapply(pay_dates, P_fun)
    (1 - P_fun(Tn) + A_sum(Tn)) / (fixed_delta * sum(P_vals))
  })
}

forecast_h_steps <- function(last_state, pars, h) {
  state_h <- last_state
  for (i in seq_len(h)) state_h <- state_transition(state_h, pars)
  model_swap_rates(state_h, pars)
}

project_to_bounds <- function(par, bounds, eps = 1e-8) {
  nm <- names(bounds$lower)
  par <- par[nm]
  pmin(pmax(par, bounds$lower + eps), bounds$upper - eps)
}

jitter_start <- function(par, bounds, scale = 0.03) {
  width <- bounds$upper - bounds$lower
  out <- par + rnorm(length(par), mean = 0, sd = scale * width)
  names(out) <- names(bounds$lower)
  project_to_bounds(out, bounds)
}

estimate_model <- function(
    Z_train,
    bounds,
    main_start,
    model_type,
    n_starts = 10,
    maxit = 250
) {
  main_start <- project_to_bounds(main_start, bounds)
  n_jitter <- max(0, n_starts - 1)
  jitter_starts <- if (n_jitter > 0) {
    replicate(n_jitter, jitter_start(main_start, bounds), simplify = FALSE)
  } else {
    list()
  }
  all_starts <- c(list(main_start), jitter_starts)

  fits <- foreach(
    st = all_starts,
    .packages = c("stats"),
    .export = c(
      "neg_loglik", "ekf_filter", "expand_pars", "state_transition", "state_Q",
      "vasicek_mean", "vasicek_var", "vasicek_second_moment",
      "measurement_function", "measurement_R", "num_jacobian",
      "safe_chol", "dt_monthly", "ois_meas_maturities",
      "A_meas_cols", "A_block_grids", "A_block_deltas"
    ),
    .errorhandling = "pass"
  ) %dopar% {
    tryCatch(
      optim(
        par = st,
        fn = neg_loglik,
        Z_obs = Z_train,
        model_type = model_type,
        method = "L-BFGS-B",
        lower = bounds$lower,
        upper = bounds$upper,
        control = list(maxit = maxit)
      ),
      error = function(e) NULL
    )
  }

  fits <- Filter(function(x) !inherits(x, "error") && !is.null(x), fits)
  if (length(fits) == 0) return(NULL)

  vals <- sapply(fits, function(x) x$value)
  best <- fits[[which.min(vals)]]
  ekf <- tryCatch(ekf_filter(best$par, Z_train, model_type), error = function(e) NULL)
  if (is.null(ekf)) return(NULL)

  list(theta_hat = best$par, pars = ekf$pars, ekf = ekf, objective = best$value)
}

# ============================================================
# FORECASTING
# ============================================================

run_rolling_forecast <- function(model_label, bounds, main_start, model_type, window = 60, h = 6) {
  method <- "rolling"
  t_grid <- window:(nrow(Z_obs) - h)
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

    set.seed(100000 + t + ifelse(model_label == "Opt. 1", 1000, 2000))
    train_idx <- (t - window + 1):t
    Z_train <- Z_obs[train_idx, , drop = FALSE]

    fit <- estimate_model(
      Z_train = Z_train,
      bounds = bounds,
      main_start = main_start,
      model_type = model_type,
      n_starts = n_starts_oos,
      maxit = optim_maxit_oos
    )

    if (!is.null(fit)) {
      last_state <- fit$ekf$x_filt[nrow(fit$ekf$x_filt), ]
      pred_swap <- forecast_h_steps(last_state, fit$pars, h)
      actual_swap <- swap_obs_real[t + h, ]

      if (any(!is.finite(pred_swap))) {
        res_list[[i]] <- NULL
        status_msg <- "FAILED - non-finite forecast"
      } else {
        res_list[[i]] <- list(
          t = t,
          date = dates[t + h],
          forecast = pred_swap,
          actual = actual_swap,
          error = actual_swap - pred_swap,
          theta_hat = fit$theta_hat,
          objective = fit$objective
        )
        status_msg <- "SUCCESS"
      }
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
  if (length(res_list) == 0) stop("No successful forecasts for ", model_label, " / ", method)

  forecast_mat <- do.call(rbind, lapply(res_list, `[[`, "forecast"))
  actual_mat   <- do.call(rbind, lapply(res_list, `[[`, "actual"))
  error_mat    <- do.call(rbind, lapply(res_list, `[[`, "error"))

  colnames(forecast_mat) <- swap_cols
  colnames(actual_mat) <- swap_cols
  colnames(error_mat) <- swap_cols

  out <- list(
    model = model_label,
    method = method,
    h = h,
    window = window,
    dates = as.Date(sapply(res_list, `[[`, "date"), origin = "1970-01-01"),
    forecast = forecast_mat,
    actual = actual_mat,
    errors = error_mat,
    rmse = sqrt(colMeans(error_mat^2, na.rm = TRUE)),
    mae = colMeans(abs(error_mat), na.rm = TRUE),
    objectives = sapply(res_list, `[[`, "objective")
  )

  cat("\nFinished", model_label, "with", length(res_list), "successful forecasts out of", n_total, "windows.\n")
  out
}

all_results <- list()
all_results$rolling5y_initial <- run_rolling_forecast("Initial Model", bounds_initial, start_initial, "initial", rolling_window, h)
all_results$rolling5y_opt1    <- run_rolling_forecast("Opt. 1", bounds_opt1, start_opt1, "opt1", rolling_window, h)

# ============================================================
# TABLES AND SAVE RESULTS
# ============================================================

make_metric_table <- function(results, metric = c("rmse", "mae")) {
  metric <- match.arg(metric)
  tab <- data.frame(maturity = swap_cols)
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

write.csv(rmse_table, "forecast_rmse_initial_opt1.csv", row.names = FALSE)
write.csv(mae_table,  "forecast_mae_initial_opt1.csv", row.names = FALSE)
saveRDS(all_results, "forecast_results_initial_opt1.rds")

error_long <- do.call(rbind, lapply(names(all_results), function(nm) {
  res <- all_results[[nm]]
  data.frame(
    result = nm,
    model = res$model,
    method = res$method,
    date = rep(res$dates, times = length(swap_cols)),
    maturity = rep(swap_cols, each = length(res$dates)),
    error = as.vector(res$errors)
  )
}))
write.csv(error_long, "forecast_errors_long_initial_opt1.csv", row.names = FALSE)

cat("\nSaved:\n")
cat("- forecast_rmse_initial_opt1.csv\n")
cat("- forecast_mae_initial_opt1.csv\n")
cat("- forecast_errors_long_initial_opt1.csv\n")
cat("- forecast_results_initial_opt1.rds\n")

# ============================================================
# PLOTS
# ============================================================

make_plot_df <- function(results) {
  do.call(rbind, lapply(names(results), function(nm) {
    res <- results[[nm]]
    do.call(rbind, lapply(seq_along(swap_cols), function(j) {
      data.frame(
        result = nm,
        model = res$model,
        date = res$dates,
        maturity = swap_cols[j],
        actual = res$actual[, j],
        forecast = res$forecast[, j]
      )
    }))
  }))
}

plot_df <- make_plot_df(all_results)
write.csv(plot_df, "forecast_actual_vs_forecast_long_initial_opt1.csv", row.names = FALSE)

# Combined plot: Actual, Initial Model forecast, Opt. 1 forecast in same panel.
actual_df <- plot_df %>%
  select(date, maturity, actual) %>%
  distinct() %>%
  mutate(series = "Actual", value = actual)

forecast_df <- plot_df %>%
  transmute(date, maturity, series = model, value = forecast)

plot_df_three_lines <- bind_rows(actual_df, forecast_df)

maturity_order <- c("Y1", "Y2", "Y3", "Y5", "Y7", "Y10", "Y15", "Y20", "Y30")
plot_df_three_lines <- plot_df_three_lines %>%
  mutate(
    maturity = factor(maturity, levels = maturity_order),
    series = factor(series, levels = c("Actual", "Initial Model", "Opt. 1"))
  )

p_three <- ggplot(
  plot_df_three_lines,
  aes(x = date, y = value, color = series, linetype = series)
) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ maturity, scales = "free_y", ncol = 3) +
  scale_color_manual(
    values = c(
      "Actual" = "black",
      "Initial Model" = "blue",
      "Opt. 1" = "red"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Actual" = "solid",
      "Initial Model" = "solid",
      "Opt. 1" = "solid"
    )
  ) +
  labs(
    title = "Actual vs 6-month-ahead forecasted EURIBOR swap rates",
    subtitle = "5-year rolling estimation window",
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
ggsave("actual_vs_forecast_initial_opt1_colored_ordered.png", p_three, width = 13, height = 9, dpi = 300)

# Forecast error plot.
error_plot_df <- plot_df %>%
  mutate(error = actual - forecast) %>%
  select(date, maturity, model, error) %>%
  mutate(
    maturity = factor(maturity, levels = maturity_order),
    model = factor(model, levels = c("Initial Model", "Opt. 1"))
  )

p_errors <- ggplot(error_plot_df, aes(x = date, y = error, color = model, linetype = model)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ maturity, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Initial Model" = "blue", "Opt. 1" = "red")) +
  labs(
    title = "6-month-ahead forecast errors",
    subtitle = "Error = actual - forecast",
    x = NULL,
    y = "Forecast error",
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

print(p_errors)
ggsave("forecast_errors_initial_opt1_colored_ordered.png", p_errors, width = 13, height = 9, dpi = 300)

# ============================================================
# DONE
# ============================================================
