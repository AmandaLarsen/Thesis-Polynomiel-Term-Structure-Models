
library(readxl)
library(openxlsx)
library(pracma)
library(ggplot2)
library(dplyr)

df <- EONIA_cleaned_before_bootstrap
df$Date <- as.Date(df$Date)

rates_are_percent <- FALSE
zero_smooth_spar <- 0.65

annual_grid <- 1:30
semi_grid   <- seq(0.5, 30, by = 0.5)

very_short_name <- "EUREON2W="   # evt. "EUREON3W="

short_tenors <- c(
  "EUREON2W="  = 2/52,
  "EUREON3W="  = 3/52,
  "EUREON1M="  = 1/12,
  "EUREON2M="  = 2/12,
  "EUREON3M="  = 3/12,
  "EUREON4M="  = 4/12,
  "EUREON5M="  = 5/12,
  "EUREON6M="  = 6/12,
  "EUREON7M="  = 7/12,
  "EUREON8M="  = 8/12,
  "EUREON9M="  = 9/12,
  "EUREON10M=" = 10/12,
  "EUREON11M=" = 11/12
)

observed_swap_tenors <- c(
  "EUREON1Y="  = 1,
  "EUREON2Y="  = 2,
  "EUREON3Y="  = 3,
  "EUREON4Y="  = 4,
  "EUREON5Y="  = 5,
  "EUREON6Y="  = 6,
  "EUREON7Y="  = 7,
  "EUREON8Y="  = 8,
  "EUREON9Y="  = 9,
  "EUREON10Y=" = 10,
  "EUREON15Y=" = 15,
  "EUREON20Y=" = 20,
  "EUREON30Y=" = 30
)

needed_cols <- c(names(short_tenors), names(observed_swap_tenors))
missing_cols <- setdiff(needed_cols, names(df))
if (length(missing_cols) > 0) {
  stop(paste("Mangler kolonner:", paste(missing_cols, collapse = ", ")))
}

if (!(very_short_name %in% names(short_tenors))) {
  stop("very_short_name findes ikke i short_tenors.")
}


# Discount factor -> zero yield
df_to_zero <- function(df_val, T) {
  ifelse(is.na(df_val) | df_val <= 0, NA_real_, -log(df_val) / T)
}

# Zero yield -> discount factor
zero_to_df <- function(z, T) {
  ifelse(is.na(z), NA_real_, exp(-z * T))
}

# Interpolér observerede swaprenter
interpolate_swap_rates_pchip <- function(rates_row, observed_swap_tenors, full_grid) {
  obs_names <- names(observed_swap_tenors)
  obs_years <- as.numeric(observed_swap_tenors)
  obs_rates <- as.numeric(rates_row[obs_names])
  
  keep <- !(is.na(obs_years) | is.na(obs_rates))
  obs_years <- obs_years[keep]
  obs_rates <- obs_rates[keep]
  
  out <- rep(NA_real_, length(full_grid))
  
  if (length(obs_rates) < 2) {
    names(out) <- paste0("T_", sprintf("%.1f", full_grid))
    return(out)
  }
  
  ord <- order(obs_years)
  obs_years <- obs_years[ord]
  obs_rates <- obs_rates[ord]
  
  if (length(obs_rates) == 2) {
    out <- approx(
      x = obs_years,
      y = obs_rates,
      xout = full_grid,
      method = "linear",
      rule = 2
    )$y
  } else {
    out[full_grid < min(obs_years)] <- obs_rates[1]
    out[full_grid > max(obs_years)] <- obs_rates[length(obs_rates)]
    
    inside <- full_grid >= min(obs_years) & full_grid <= max(obs_years)
    out[inside] <- pracma::pchip(obs_years, obs_rates, full_grid[inside])
  }

  match_idx <- match(obs_years, full_grid)
  out[match_idx[!is.na(match_idx)]] <- obs_rates
  
  names(out) <- paste0("T_", sprintf("%.1f", full_grid))
  out
}

# Bootstrap
bootstrap_dfs_from_par <- function(par_rates, df_first, delta = 1) {
  n <- length(par_rates)
  dfs <- rep(NA_real_, n)
  names(dfs) <- names(par_rates)
  
  dfs[1] <- df_first
  
  if (is.na(df_first) || !is.finite(df_first) || df_first <= 0) {
    return(dfs)
  }
  
  for (k in 2:n) {
    Sk <- par_rates[k]
    
    if (is.na(Sk) || any(is.na(dfs[1:(k - 1)]))) {
      dfs[k] <- NA_real_
    } else {
      dfs[k] <- (1 - Sk * delta * sum(dfs[1:(k - 1)])) / (1 + Sk * delta)
      
      if (!is.finite(dfs[k]) || dfs[k] <= 0) {
        dfs[k] <- NA_real_
      }
    }
  }
  
  dfs
}

# Interpolation
interp_zero_linear <- function(zero_vec, grid, target_grid) {
  keep <- !is.na(zero_vec)
  out <- rep(NA_real_, length(target_grid))
  
  if (sum(keep) < 2) {
    names(out) <- paste0("T_", sprintf("%.2f", target_grid))
    return(out)
  }
  
  out <- approx(
    x = grid[keep],
    y = zero_vec[keep],
    xout = target_grid,
    method = "linear",
    rule = 2
  )$y
  
  names(out) <- paste0("T_", sprintf("%.2f", target_grid))
  out
}

# Smooth zero curve
smooth_zero_curve <- function(zero_vec, grid, spar = 0.65) {
  keep <- !is.na(zero_vec)
  out <- rep(NA_real_, length(grid))
  
  if (sum(keep) < 4) {
    names(out) <- names(zero_vec)
    return(out)
  }
  
  fit <- smooth.spline(x = grid[keep], y = zero_vec[keep], spar = spar)
  sm <- predict(fit, x = grid)$y
  names(sm) <- names(zero_vec)
  sm
}

# Udregn halvårlige par swaprenter fra discount factors
par_from_dfs_semiannual <- function(dfs, delta = 0.5) {
  n <- length(dfs)
  out <- rep(NA_real_, n)
  names(out) <- names(dfs)
  
  out[1] <- NA_real_
  
  for (k in 2:n) {
    if (any(is.na(dfs[1:k]))) next
    denom <- delta * sum(dfs[1:k])
    if (denom <= 0) next
    out[k] <- (1 - dfs[k]) / denom
  }
  
  out
}


build_curve_zero_smooth <- function(rates_row,
                                    short_tenors,
                                    observed_swap_tenors,
                                    annual_grid = 1:30,
                                    semi_grid = seq(0.5, 30, by = 0.5),
                                    zero_smooth_spar = 0.65,
                                    rates_are_percent = FALSE) {
  
  # -------- Korte løbetider (<1 år) --------
  short_rates <- as.numeric(rates_row[names(short_tenors)])
  short_maturities <- as.numeric(short_tenors)
  names(short_rates) <- names(short_tenors)
  
  if (rates_are_percent) short_rates <- short_rates / 100
  
  short_df <- 1 / (1 + short_rates * short_maturities)
  names(short_df) <- names(short_tenors)
  
  short_zero <- -log(short_df) / short_maturities
  names(short_zero) <- names(short_tenors)
  
  # -------- Interpolér årlige punkter --------
  annual_par_raw <- interpolate_swap_rates_pchip(
    rates_row = rates_row,
    observed_swap_tenors = observed_swap_tenors,
    full_grid = annual_grid
  )
  
  if (rates_are_percent) annual_par_raw <- annual_par_raw / 100
  names(annual_par_raw) <- paste0("T_", sprintf("%.1f", annual_grid))
  
  # -------- Bootstrap på årlige punkter --------
  df_1y <- 1 / (1 + annual_par_raw["T_1.0"] * 1)
  
  annual_df_raw <- bootstrap_dfs_from_par(
    par_rates = annual_par_raw,
    df_first = df_1y,
    delta = 1
  )
  
  annual_zero_raw <- mapply(df_to_zero, annual_df_raw, annual_grid)
  names(annual_zero_raw) <- paste0("T_", sprintf("%.1f", annual_grid))
  
  # Interpolér halvårlige punkter
  semi_zero_raw <- interp_zero_linear(
    zero_vec = annual_zero_raw,
    grid = annual_grid,
    target_grid = semi_grid
  )
  names(semi_zero_raw) <- paste0("T_", sprintf("%.1f", semi_grid))
  
  semi_zero_raw["T_0.5"] <- short_zero["EUREON6M="]
  
  # -------- Smooth halvårlig zero curve --------
  semi_zero_smooth <- smooth_zero_curve(
    zero_vec = semi_zero_raw,
    grid = semi_grid,
    spar = zero_smooth_spar
  )
  names(semi_zero_smooth) <- paste0("T_", sprintf("%.1f", semi_grid))
  
  # -------- Discount factors fra smooth zero curve --------
  semi_df_smooth <- mapply(zero_to_df, semi_zero_smooth, semi_grid)
  names(semi_df_smooth) <- paste0("T_", sprintf("%.1f", semi_grid))
  
  # -------- Halvårlige par swaprenter fra smooth DF --------
  semi_par_smooth <- par_from_dfs_semiannual(semi_df_smooth, delta = 0.5)
  
  list(
    short_df = short_df,
    short_zero = short_zero,
    annual_par_raw = annual_par_raw,
    annual_df_raw = annual_df_raw,
    annual_zero_raw = annual_zero_raw,
    semi_zero_raw = semi_zero_raw,
    semi_zero = semi_zero_smooth,
    semi_df = semi_df_smooth,
    semi_par = semi_par_smooth
  )
}

# ============================================
# HELE DATASETTET
# ============================================
curve_list <- apply(df[, needed_cols], 1, function(row) {
  build_curve_zero_smooth(
    rates_row = row,
    short_tenors = short_tenors,
    observed_swap_tenors = observed_swap_tenors,
    annual_grid = annual_grid,
    semi_grid = semi_grid,
    zero_smooth_spar = zero_smooth_spar,
    rates_are_percent = rates_are_percent
  )
})

# ============================================
# GEM OUTPUT
# ============================================
bond_prices_df <- data.frame(Date = df$Date)
zero_rates_df  <- data.frame(Date = df$Date)

for (nm in names(short_tenors)) {
  bond_prices_df[[paste0("P_", nm)]] <- sapply(curve_list, function(x) x$short_df[nm])
  zero_rates_df[[paste0("Z_", nm)]]  <- sapply(curve_list, function(x) x$short_zero[nm])
}

for (tm in annual_grid) {
  nm <- paste0("T_", sprintf("%.1f", tm))
  
  bond_prices_df[[paste0("Praw_", nm)]] <- sapply(curve_list, function(x) x$annual_df_raw[nm])
  zero_rates_df[[paste0("Zraw_", nm)]]  <- sapply(curve_list, function(x) x$annual_zero_raw[nm])
  zero_rates_df[[paste0("Sraw_", nm)]]  <- sapply(curve_list, function(x) x$annual_par_raw[nm])
}

for (tm in semi_grid) {
  nm <- paste0("T_", sprintf("%.1f", tm))
  
  bond_prices_df[[paste0("P_", nm)]] <- sapply(curve_list, function(x) x$semi_df[nm])
  zero_rates_df[[paste0("Z_", nm)]]  <- sapply(curve_list, function(x) x$semi_zero[nm])
  zero_rates_df[[paste0("S_", nm)]]  <- sapply(curve_list, function(x) x$semi_par[nm])
}

full_bond_prices_df <- data.frame(Date = df$Date)
full_zero_rates_df  <- data.frame(Date = df$Date)

full_bond_prices_df[[paste0("P_", very_short_name)]] <- sapply(curve_list, function(x) x$short_df[very_short_name])
full_zero_rates_df[[paste0("Z_", very_short_name)]]  <- sapply(curve_list, function(x) x$short_zero[very_short_name])

full_bond_prices_df[["P_T_0.5"]] <- sapply(curve_list, function(x) x$semi_df["T_0.5"])
full_zero_rates_df[["Z_T_0.5"]]  <- sapply(curve_list, function(x) x$semi_zero["T_0.5"])

for (tm in semi_grid[semi_grid >= 1]) {
  nm <- paste0("T_", sprintf("%.1f", tm))
  full_bond_prices_df[[paste0("P_", nm)]] <- sapply(curve_list, function(x) x$semi_df[nm])
  full_zero_rates_df[[paste0("Z_", nm)]]  <- sapply(curve_list, function(x) x$semi_zero[nm])
}

write.xlsx(bond_prices_df,      "~/Desktop/zero_smooth_bond_prices.xlsx", overwrite = TRUE)
write.xlsx(zero_rates_df,       "~/Desktop/zero_smooth_zero_rates.xlsx", overwrite = TRUE)
write.xlsx(full_bond_prices_df, "~/Desktop/full_curve_bond_prices.xlsx", overwrite = TRUE)
write.xlsx(full_zero_rates_df,  "~/Desktop/full_curve_zero_rates.xlsx", overwrite = TRUE)

# ============================================
# PLOT
# ============================================
i <- 10

very_short_mat <- as.numeric(short_tenors[very_short_name])
plot_grid <- c(very_short_mat, 0.5, semi_grid[semi_grid >= 1])

z_very_short <- curve_list[[i]]$short_zero[very_short_name]
p_very_short <- curve_list[[i]]$short_df[very_short_name]

z_6m_raw <- curve_list[[i]]$semi_zero_raw["T_0.5"]
z_6m_sm  <- curve_list[[i]]$semi_zero["T_0.5"]
p_6m_sm  <- curve_list[[i]]$semi_df["T_0.5"]

z_long_raw <- sapply(semi_grid[semi_grid >= 1], function(tm) {
  curve_list[[i]]$semi_zero_raw[paste0("T_", sprintf("%.1f", tm))]
})
z_long_sm <- sapply(semi_grid[semi_grid >= 1], function(tm) {
  curve_list[[i]]$semi_zero[paste0("T_", sprintf("%.1f", tm))]
})
p_long_sm <- sapply(semi_grid[semi_grid >= 1], function(tm) {
  curve_list[[i]]$semi_df[paste0("T_", sprintf("%.1f", tm))]
})

z_raw_full <- c(z_very_short, z_6m_raw, z_long_raw)
z_sm_full  <- c(z_very_short, z_6m_sm,  z_long_sm)
p_sm_full  <- c(p_very_short, p_6m_sm,  p_long_sm)

par(mfrow = c(1, 2))
plot(plot_grid, z_raw_full, type = "b",
     main = "Zero yields: short point + 6M",
     xlab = "Maturity (years)", ylab = "Zero yield")
lines(plot_grid, z_sm_full, type = "b", lty = 2)
legend("topright", legend = c("Raw/observed", "Smoothed"), lty = c(1, 2), bty = "n")

plot(plot_grid, p_sm_full, type = "b",
     main = "Discount curve",
     xlab = "Maturity (years)", ylab = "Discount factor")
par(mfrow = c(1, 1))

# ============================================
# DISCOUNT CURVE
# ============================================
date1 <- as.Date("2008-10-17")
date2 <- as.Date("2016-09-16")

i1 <- which(df$Date == date1)
i2 <- which(df$Date == date2)

if (length(i1) == 0 || length(i2) == 0) {
  stop("En eller begge datoer findes ikke i datasættet.")
}

get_full_price_curve <- function(curve_obj, very_short_name, short_tenors, semi_grid) {
  c(
    curve_obj$short_df[very_short_name],
    curve_obj$semi_df["T_0.5"],
    sapply(semi_grid[semi_grid >= 1], function(tm) {
      curve_obj$semi_df[paste0("T_", sprintf("%.1f", tm))]
    })
  )
}

get_full_zero_curve_raw <- function(curve_obj, very_short_name, short_tenors, semi_grid) {
  c(
    curve_obj$short_zero[very_short_name],
    curve_obj$semi_zero_raw["T_0.5"],
    sapply(semi_grid[semi_grid >= 1], function(tm) {
      curve_obj$semi_zero_raw[paste0("T_", sprintf("%.1f", tm))]
    })
  )
}

get_full_zero_curve_smoothed <- function(curve_obj, very_short_name, short_tenors, semi_grid) {
  c(
    curve_obj$short_zero[very_short_name],
    curve_obj$semi_zero["T_0.5"],
    sapply(semi_grid[semi_grid >= 1], function(tm) {
      curve_obj$semi_zero[paste0("T_", sprintf("%.1f", tm))]
    })
  )
}

full_maturity_grid <- c(as.numeric(short_tenors[very_short_name]), 0.5, semi_grid[semi_grid >= 1])

plot_df <- bind_rows(
  data.frame(
    Date = df$Date[i1],
    maturity = full_maturity_grid,
    price = as.numeric(get_full_price_curve(curve_list[[i1]], very_short_name, short_tenors, semi_grid))
  ),
  data.frame(
    Date = df$Date[i2],
    maturity = full_maturity_grid,
    price = as.numeric(get_full_price_curve(curve_list[[i2]], very_short_name, short_tenors, semi_grid))
  )
)

plot_df$Date_lab <- paste("Discount curve / bond prices on", plot_df$Date)

ggplot(plot_df, aes(x = maturity, y = price)) +
  geom_line(color = "black", linewidth = 0.5) +
  geom_point(color = "black", size = 1.4) +
  facet_wrap(~ Date_lab, ncol = 1, scales = "free_y") +
  labs(
    x = "Maturity (years)",
    y = "Bond price P(0,T)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# ============================================
# SMOOTHED ZERO CURVES
# ============================================
plot_df_zero_both <- bind_rows(
  data.frame(
    Date = df$Date[i1],
    maturity = full_maturity_grid,
    value = as.numeric(get_full_zero_curve_raw(curve_list[[i1]], very_short_name, short_tenors, semi_grid)),
    Curve = "Raw"
  ),
  data.frame(
    Date = df$Date[i1],
    maturity = full_maturity_grid,
    value = as.numeric(get_full_zero_curve_smoothed(curve_list[[i1]], very_short_name, short_tenors, semi_grid)),
    Curve = "Smoothed"
  ),
  data.frame(
    Date = df$Date[i2],
    maturity = full_maturity_grid,
    value = as.numeric(get_full_zero_curve_raw(curve_list[[i2]], very_short_name, short_tenors, semi_grid)),
    Curve = "Raw"
  ),
  data.frame(
    Date = df$Date[i2],
    maturity = full_maturity_grid,
    value = as.numeric(get_full_zero_curve_smoothed(curve_list[[i2]], very_short_name, short_tenors, semi_grid)),
    Curve = "Smoothed"
  )
)

plot_df_zero_both$Panel <- paste("Zero yield curve on", plot_df_zero_both$Date)

p_zero_both <- ggplot(plot_df_zero_both, aes(x = maturity, y = value, linetype = Curve, shape = Curve)) +
  geom_line(color = "black", linewidth = 0.6) +
  geom_point(color = "black", size = 1.6) +
  facet_wrap(~ Panel, ncol = 1, scales = "free_y") +
  labs(
    title = "Zero yield curves",
    x = "Maturity (years)",
    y = "Zero yield"
  ) +
  scale_linetype_manual(values = c("Raw" = "dashed", "Smoothed" = "solid")) +
  scale_shape_manual(values = c("Raw" = 1, "Smoothed" = 16)) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

print(p_zero_both)
































