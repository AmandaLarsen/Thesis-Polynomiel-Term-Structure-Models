###############################################################
# VOLATILITY COMPRESSION:
###############################################################

library(dplyr)

baseline_workspace <- "~/Desktop/combined_vasicek_cir_ekf_simulation_workspace.RData"
initial_workspace  <- "~/Desktop/ekf_workspace_first_model.RData"
ukf_workspace_uden <- "~/Desktop/ukf_workspace_tilfoej_a_uden_udvidet_int.RData"
ukf_workspace_med  <- "~/Desktop/ukf_workspace_tilfoej_a_med_andet_udvidet_int.RData"

OUTPUT_DIR <- "~/Desktop"

base_env <- new.env()
load(baseline_workspace, envir = base_env)

init_env <- new.env()
load(initial_workspace, envir = init_env)

ukf_env_uden <- new.env()
load(ukf_workspace_uden, envir = ukf_env_uden)

ukf_env_med <- new.env()
load(ukf_workspace_med, envir = ukf_env_med)

# ============================================================
#  Helpers
# ============================================================

dt_monthly <- 1 / 12
fixed_delta <- 0.5
swap_maturities <- c(1, 2, 3, 5, 7, 10, 15, 20, 30)

vasicek_mean <- function(x, tau, kappa, theta) {
  theta + exp(-kappa * tau) * (x - theta)
}

vasicek_var <- function(tau, kappa, sigma) {
  if (kappa < 1e-8) {
    return(sigma^2 * tau)
  }
  sigma^2 / (2 * kappa) * (1 - exp(-2 * kappa * tau))
}

vasicek_second_moment <- function(x, tau, kappa, theta, sigma) {
  m <- vasicek_mean(x, tau, kappa, theta)
  v <- vasicek_var(tau, kappa, sigma)
  v + m^2
}

safe_chol <- function(S, jitter = 1e-10, max_tries = 8) {
  S <- (S + t(S)) / 2
  for (i in 0:max_tries) {
    out <- tryCatch(
      chol(S + diag(jitter * 10^i, nrow(S))),
      error = function(e) NULL
    )
    if (!is.null(out)) return(out)
  }
  stop("Failed")
}

num_jacobian <- function(f, x, eps = 1e-6) {
  n <- length(x)
  fx <- f(x)
  m <- length(fx)
  J <- matrix(0, m, n)
  for (i in seq_len(n)) {
    xp <- x
    xm <- x
    xp[i] <- xp[i] + eps
    xm[i] <- xm[i] - eps
    J[, i] <- (f(xp) - f(xm)) / (2 * eps)
  }
  J
}

inject_workspace_helpers <- function(env) {
  helper_names <- c(
    "vasicek_mean",
    "vasicek_var",
    "vasicek_second_moment",
    "safe_chol",
    "num_jacobian"
  )
  
  for (nm in helper_names) {
    assign(nm, get(nm, envir = .GlobalEnv), envir = env)
  }
  
  object_names <- ls(env)
  for (obj_nm in object_names) {
    obj <- get(obj_nm, envir = env)
    if (is.function(obj)) {
      fun_env <- environment(obj)
      for (nm in helper_names) {
        assign(nm, get(nm, envir = .GlobalEnv), envir = fun_env)
      }
    }
  }
  invisible(TRUE)
}

inject_workspace_helpers(base_env)
inject_workspace_helpers(init_env)
inject_workspace_helpers(ukf_env_uden)
inject_workspace_helpers(ukf_env_med)

check_missing_functions <- function(env) {
  if (!requireNamespace("codetools", quietly = TRUE)) {
    stop("Package 'codetools' is needed for this diagnostic.")
  }
  
  fun_names <- ls(env)[vapply(ls(env), function(nm) {
    is.function(get(nm, envir = env))
  }, logical(1))]
  
  missing <- unique(unlist(lapply(fun_names, function(fn) {
    f <- get(fn, envir = env)
    globals <- codetools::findGlobals(f, merge = FALSE)$functions
    globals[!vapply(globals, exists, logical(1), envir = environment(f), inherits = TRUE)]
  })))
  
  sort(missing)
}

# ============================================================
# BASELINE
# ============================================================

Y_real <- base_env$Y_real
dates  <- as.Date(base_env$dates)
taus   <- base_env$taus
dt     <- base_env$dt
pay_freq <- base_env$pay_freq

res_vas_real <- base_env$res_vas_real
res_cir_real <- base_env$res_cir_real

Y_pred_vas <- res_vas_real$filtered$Y_pred
Y_pred_cir <- res_cir_real$filtered$Y_pred

swap_cols <- paste0("Y", taus)


x_max <- 7.5
realized_window <- 12
realized_annualization <- sqrt(12)
model_annualization <- sqrt(12)

rolling_sd_matrix <- function(X, window) {
  out <- matrix(NA_real_, nrow(X), ncol(X))
  for (j in seq_len(ncol(X))) {
    for (i in window:nrow(X)) {
      out[i, j] <- sd(X[(i - window + 1):i, j], na.rm = TRUE)
    }
  }
  out
}

annual_payment_times <- function(T, pay_freq = 1) {
  delta <- 1 / pay_freq
  seq(delta, T, by = delta)
}

# ============================================================
#  VASICEK/CIR MODEL VOLATILITY
# ============================================================

vasicek_B <- function(tau, kappa) {
  (1 - exp(-kappa * tau)) / kappa
}

vasicek_A <- function(tau, kappa, theta, sigma) {
  B <- vasicek_B(tau, kappa)
  (theta - sigma^2 / (2 * kappa^2)) * (B - tau) -
    (sigma^2 / (4 * kappa)) * B^2
}

vasicek_P <- function(tau, r, kappa, theta, sigma) {
  exp(vasicek_A(tau, kappa, theta, sigma) -
        vasicek_B(tau, kappa) * r)
}

swap_rate_dr_vas <- function(r, T, kappa, theta, sigma, pay_freq = 1) {
  delta <- 1 / pay_freq
  pay_times <- annual_payment_times(T, pay_freq)
  
  Pn <- vasicek_P(T, r, kappa, theta, sigma)
  dPn_dr <- -vasicek_B(T, kappa) * Pn
  
  Pk <- vasicek_P(pay_times, r, kappa, theta, sigma)
  dPk_dr <- -vasicek_B(pay_times, kappa) * Pk
  
  N <- 1 - Pn
  D <- sum(delta * Pk)
  
  dN_dr <- -dPn_dr
  dD_dr <- sum(delta * dPk_dr)
  
  (dN_dr * D - N * dD_dr) / D^2
}

vasicek_cond_var <- function(kappa, sigma, dt) {
  sigma^2 / (2 * kappa) * (1 - exp(-2 * kappa * dt))
}

model_vol_vasicek <- function(r_states, taus, params, dt, pay_freq = 1) {
  q <- vasicek_cond_var(params$kappa, params$sigma, dt)
  vol <- matrix(NA_real_, length(r_states), length(taus))
  
  for (i in seq_along(r_states)) {
    H <- sapply(taus, function(T) {
      swap_rate_dr_vas(r_states[i], T,
                       params$kappa, params$theta, params$sigma,
                       pay_freq)
    })
    vol[i, ] <- abs(H) * sqrt(q) * model_annualization * 10000
  }
  vol
}

cir_gamma <- function(kappa, sigma) {
  sqrt(kappa^2 + 2 * sigma^2)
}

cir_B <- function(tau, kappa, sigma) {
  gamma <- cir_gamma(kappa, sigma)
  2 * (exp(gamma * tau) - 1) /
    ((gamma + kappa) * (exp(gamma * tau) - 1) + 2 * gamma)
}

cir_A <- function(tau, kappa, theta, sigma) {
  gamma <- cir_gamma(kappa, sigma)
  denominator <- (gamma + kappa) * (exp(gamma * tau) - 1) + 2 * gamma
  base <- 2 * gamma * exp((kappa + gamma) * tau / 2) / denominator
  power <- 2 * kappa * theta / sigma^2
  base^power
}

cir_P <- function(tau, r, kappa, theta, sigma) {
  r <- pmax(r, 1e-12)
  cir_A(tau, kappa, theta, sigma) *
    exp(-cir_B(tau, kappa, sigma) * r)
}

swap_rate_dr_cir <- function(r, T, kappa, theta, sigma, pay_freq = 1) {
  r <- max(r, 1e-12)
  delta <- 1 / pay_freq
  pay_times <- annual_payment_times(T, pay_freq)
  
  Pn <- cir_P(T, r, kappa, theta, sigma)
  dPn_dr <- -cir_B(T, kappa, sigma) * Pn
  
  Pk <- cir_P(pay_times, r, kappa, theta, sigma)
  dPk_dr <- -cir_B(pay_times, kappa, sigma) * Pk
  
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
  
  r * sigma^2 / kappa * (exp1 - exp2) +
    theta * sigma^2 / (2 * kappa) * (1 - exp1)^2
}

model_vol_cir <- function(r_states, taus, params, dt, pay_freq = 1) {
  vol <- matrix(NA_real_, length(r_states), length(taus))
  
  for (i in seq_along(r_states)) {
    q <- cir_cond_var(r_states[i],
                      params$kappa, params$theta, params$sigma, dt)
    
    H <- sapply(taus, function(T) {
      swap_rate_dr_cir(r_states[i], T,
                       params$kappa, params$theta, params$sigma,
                       pay_freq)
    })
    
    vol[i, ] <- abs(H) * sqrt(q) * model_annualization * 10000
  }
  vol
}

# ============================================================
# BASELINE VOL
# ============================================================

dY_real <- rbind(rep(NA_real_, ncol(Y_real)), diff(Y_real))

realized_vol_real <- rolling_sd_matrix(dY_real, realized_window) *
  realized_annualization * 10000

vol_vas_pred <- model_vol_vasicek(
  res_vas_real$filtered$r_pred,
  taus,
  res_vas_real$params,
  dt,
  pay_freq
)

vol_cir_pred <- model_vol_cir(
  res_cir_real$filtered$r_pred,
  taus,
  res_cir_real$params,
  dt,
  pay_freq
)

# ============================================================
# INITIAL EKF MODEL VOLATILITY
# ============================================================

stopifnot(
  exists("ekf_out_real", envir = init_env),
  exists("swap_pred_real", envir = init_env),
  exists("dates", envir = init_env),
  exists("model_swap_rates", envir = init_env),
  exists("num_jacobian", envir = init_env),
  exists("state_Q", envir = init_env)
)

init_dates <- as.Date(init_env$dates)
Y_pred_init <- init_env$swap_pred_real

vol_init_pred <- matrix(
  NA_real_,
  nrow = nrow(init_env$ekf_out_real$x_pred),
  ncol = length(taus)
)

for (i in seq_len(nrow(init_env$ekf_out_real$x_pred))) {
  
  state_i <- init_env$ekf_out_real$x_pred[i, ]
  
  J_i <- init_env$num_jacobian(
    function(xx) {
      init_env$model_swap_rates(
        state = xx,
        pars  = init_env$ekf_out_real$pars
      )
    },
    state_i
  )
  
  Q_i <- init_env$state_Q(init_env$ekf_out_real$pars)
  
  Sigma_swap_i <- J_i %*% Q_i %*% t(J_i)
  
  vol_init_pred[i, ] <- sqrt(pmax(diag(Sigma_swap_i), 0)) *
    model_annualization * 10000
}

colnames(Y_pred_init) <- swap_cols
colnames(vol_init_pred) <- swap_cols

init_pred_df <- data.frame(Date = init_dates, Y_pred_init)
init_vol_df  <- data.frame(Date = init_dates, vol_init_pred)

main_dates_df <- data.frame(Date = dates)

Y_pred_init_aligned <- main_dates_df %>%
  left_join(init_pred_df, by = "Date") %>%
  select(all_of(swap_cols)) %>%
  as.matrix()

vol_init_pred_aligned <- main_dates_df %>%
  left_join(init_vol_df, by = "Date") %>%
  select(all_of(swap_cols)) %>%
  as.matrix()

# ============================================================
# UKF 
# ============================================================

compute_ukf_vol <- function(env, taus, swap_cols, main_dates_df) {
  
  stopifnot(
    exists("ukf_out_real", envir = env),
    exists("swap_pred_real", envir = env),
    exists("dates", envir = env),
    exists("model_swap_rates", envir = env),
    exists("ukf_sigma_points", envir = env),
    exists("state_Q", envir = env)
  )
  
  ukf_dates <- as.Date(env$dates)
  Y_pred <- env$swap_pred_real
  
  ukf_sigma_swap_cov_local <- function(x_pred, P_pred, pars) {
    
    sig <- env$ukf_sigma_points(
      x = x_pred,
      P = P_pred,
      alpha_ukf = env$ukf_out_real$alpha_ukf,
      beta_ukf  = env$ukf_out_real$beta_ukf,
      kappa_ukf = env$ukf_out_real$kappa_ukf
    )
    
    Xsig <- sig$Xsig
    Wm <- sig$Wm
    Wc <- sig$Wc
    
    Ssig <- t(apply(Xsig, 1, function(xx) {
      env$model_swap_rates(
        state = xx,
        pars  = pars
      )
    }))
    
    Sbar <- as.numeric(colSums(Ssig * Wm))
    
    Pss <- matrix(0, ncol(Ssig), ncol(Ssig))
    
    for (i in seq_len(nrow(Ssig))) {
      ds <- matrix(Ssig[i, ] - Sbar, ncol = 1)
      Pss <- Pss + Wc[i] * (ds %*% t(ds))
    }
    
    (Pss + t(Pss)) / 2
  }
  
  vol_pred <- matrix(
    NA_real_,
    nrow = nrow(env$ukf_out_real$x_pred),
    ncol = length(taus)
  )
  
  for (i in seq_len(nrow(env$ukf_out_real$x_pred))) {
    
    Pss_i <- ukf_sigma_swap_cov_local(
      x_pred = env$ukf_out_real$x_pred[i, ],
      P_pred = env$state_Q(env$ukf_out_real$pars),
      pars   = env$ukf_out_real$pars
    )
    
    vol_pred[i, ] <- sqrt(pmax(diag(Pss_i), 0)) *
      model_annualization * 10000
  }
  
  colnames(Y_pred) <- swap_cols
  colnames(vol_pred) <- swap_cols
  
  pred_df <- data.frame(Date = ukf_dates, Y_pred)
  vol_df  <- data.frame(Date = ukf_dates, vol_pred)
  
  Y_aligned <- main_dates_df %>%
    left_join(pred_df, by = "Date") %>%
    select(all_of(swap_cols)) %>%
    as.matrix()
  
  vol_aligned <- main_dates_df %>%
    left_join(vol_df, by = "Date") %>%
    select(all_of(swap_cols)) %>%
    as.matrix()
  
  list(
    Y_aligned = Y_aligned,
    vol_aligned = vol_aligned
  )
}

ukf_res_uden <- compute_ukf_vol(
  env = ukf_env_uden,
  taus = taus,
  swap_cols = swap_cols,
  main_dates_df = main_dates_df
)

ukf_res_med <- compute_ukf_vol(
  env = ukf_env_med,
  taus = taus,
  swap_cols = swap_cols,
  main_dates_df = main_dates_df
)

Y_pred_ukf_uden_aligned <- ukf_res_uden$Y_aligned
vol_ukf_uden_pred_aligned <- ukf_res_uden$vol_aligned

Y_pred_ukf_med_aligned <- ukf_res_med$Y_aligned
vol_ukf_med_pred_aligned <- ukf_res_med$vol_aligned

# ============================================================
# PLOTS
# ============================================================
op <- par(no.readonly = TRUE)

layout(
  matrix(c(
    1,2,3,
    4,5,6,
    7,8,9,
    10,10,10
  ), nrow = 4, byrow = TRUE),
  heights = c(1,1,1,0.18)
)

cols <- c(
  "Data" = "black",
  "Vasicek" = "#ff69b4",
  "CIR" = "#1f77b4",
  "Initial model" = "orange",
  "Optimization 4 UKF" = "#19B52D",
  "Optimization 5 UKF" = "#8A2BE2"
)

par(mar = c(4,4,3,1))

for (j in seq_along(taus)) {
  
  x_data <- Y_real[, j] * 100
  x_vas  <- Y_pred_vas[, j] * 100
  x_cir  <- Y_pred_cir[, j] * 100
  x_init <- Y_pred_init_aligned[, j] * 100
  x_ukf_uden <- Y_pred_ukf_uden_aligned[, j] * 100
  x_ukf_med  <- Y_pred_ukf_med_aligned[, j] * 100
  
  idx_data <- is.finite(x_data) &
    is.finite(realized_vol_real[, j]) &
    x_data <= x_max
  
  idx_vas <- is.finite(x_vas) &
    is.finite(vol_vas_pred[, j]) &
    x_vas <= x_max
  
  idx_cir <- is.finite(x_cir) &
    is.finite(vol_cir_pred[, j]) &
    x_cir <= x_max
  
  idx_init <- is.finite(x_init) &
    is.finite(vol_init_pred_aligned[, j]) &
    x_init <= x_max
  
  idx_ukf_uden <- is.finite(x_ukf_uden) &
    is.finite(vol_ukf_uden_pred_aligned[, j]) &
    x_ukf_uden <= x_max
  
  idx_ukf_med <- is.finite(x_ukf_med) &
    is.finite(vol_ukf_med_pred_aligned[, j]) &
    x_ukf_med <= x_max
  
  y_all <- c(
    realized_vol_real[idx_data, j],
    vol_vas_pred[idx_vas, j],
    vol_cir_pred[idx_cir, j],
    vol_init_pred_aligned[idx_init, j],
    vol_ukf_uden_pred_aligned[idx_ukf_uden, j],
    vol_ukf_med_pred_aligned[idx_ukf_med, j]
  )
  
  ylim_j <- range(y_all, na.rm = TRUE)
  
  plot(
    x_data[idx_data],
    realized_vol_real[idx_data, j],
    pch = 16,
    cex = 1.1,
    col = adjustcolor(cols["Data"], alpha.f = 0.35),
    xlim = c(-0.5, 4.5),
    ylim = ylim_j,
    main = paste0(taus[j], "Y"),
    xlab = "Swap rate (%)",
    ylab = "Volatility (bp)"
  )
  
  grid(col = "grey85")
  
  points(
    x_vas[idx_vas],
    vol_vas_pred[idx_vas, j],
    pch = 16,
    cex = 1.1,
    col = adjustcolor(cols["Vasicek"], alpha.f = 0.45)
  )
  
  points(
    x_cir[idx_cir],
    vol_cir_pred[idx_cir, j],
    pch = 16,
    cex = 1.1,
    col = adjustcolor(cols["CIR"], alpha.f = 0.45)
  )
  
  points(
    x_init[idx_init],
    vol_init_pred_aligned[idx_init, j],
    pch = 16,
    cex = 1.1,
    col = adjustcolor(cols["Initial model"], alpha.f = 0.45)
  )
  
  points(
    x_ukf_uden[idx_ukf_uden],
    vol_ukf_uden_pred_aligned[idx_ukf_uden, j],
    pch = 16,
    cex = 1.1,
    col = adjustcolor(cols["Optimization 4 UKF"], alpha.f = 0.45)
  )
  
  points(
    x_ukf_med[idx_ukf_med],
    vol_ukf_med_pred_aligned[idx_ukf_med, j],
    pch = 16,
    cex = 1.1,
    col = adjustcolor(cols["Optimization 5 UKF"], alpha.f = 0.45)
  )
}

par(mar = c(0,0,0,0))

plot.new()

legend(
  "center",
  legend = c(
    "Data",
    "Vasicek",
    "CIR",
    "Initial model",
    "Optimization 4 UKF",
    "Optimization 5 UKF"
  ),
  col = c(
    "black",
    "#ff69b4",
    "#1f77b4",
    "orange",
    "#19B52D",
    "#8A2BE2"
  ),
  pch = 16,
  pt.cex = 1.3,
  horiz = TRUE,
  bty = "n"
)

par(op)