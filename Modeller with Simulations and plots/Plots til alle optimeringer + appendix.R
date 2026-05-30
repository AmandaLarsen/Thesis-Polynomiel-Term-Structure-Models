# ============================================================
# PLOTS TIL ALLE OPT
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

workspace_files <- c(
  "Initial model" =
    "~/Desktop/ekf_workspace_first_model.RData",
  
  "Optimization 1" =
    "~/Desktop/ekf_workspace_correlated_model.RData",
  
  "Optimization 2" =
    "~/Desktop/ekf_workspace_udvidet_int_model.RData",
  
  "Optimization 3" =
    "~/Desktop/ekf_workspace_tilfoej_a.RData",
  
  "Optimization 4" =
    "~/Desktop/ekf_workspace_tilfoej_a_uden_udevidet.RData",
  
  "Optimization 5" =
    "~/Desktop/ekf_workspace_tilfoej_a_med_andet_udvidet_int.RData"
)


load_workspace <- function(file, model_name) {
  
  e <- new.env()
  load(file, envir = e)
  
  needed <- c(
    "dates",
    "swap_cols",
    "swap_maturities",
    "swap_obs_real",
    "swap_pred_real",
    "resid_swap_pred_real",
    "ekf_out_real",
    "Z_obs"
  )
  
  missing <- setdiff(needed, ls(e))
  
  if (length(missing) > 0) {
    stop(
      "Mangler i ", model_name, ": ",
      paste(missing, collapse = ", ")
    )
  }
  
  swap_df  <- as.data.frame(e$swap_obs_real)
  pred_df  <- as.data.frame(e$swap_pred_real)
  resid_df <- as.data.frame(e$resid_swap_pred_real)
  
  colnames(swap_df)  <- e$swap_cols
  colnames(pred_df)  <- e$swap_cols
  colnames(resid_df) <- e$swap_cols
  
  plot_df <- tibble(Date = e$dates) %>%
    bind_cols(
      swap_df  %>% rename_with(~ paste0("obs_", .x)),
      pred_df  %>% rename_with(~ paste0("pred_", .x)),
      resid_df %>% rename_with(~ paste0("resid_", .x))
    ) %>%
    slice(3:n()) %>%
    pivot_longer(
      cols = -Date,
      names_to = c("type", "maturity"),
      names_pattern = "(obs|pred|resid)_(.*)",
      values_to = "value"
    ) %>%
    mutate(model = model_name)
  
  param_df <- tibble(
    model = model_name,
    parameter = names(e$ekf_out_real$pars),
    estimate = as.numeric(e$ekf_out_real$pars),
    objective_value = if ("best_fit" %in% ls(e)) e$best_fit$value else NA_real_,
    best_start = if ("best_idx" %in% ls(e)) e$best_idx else NA_integer_
  )
  
  list(
    plot_df = plot_df,
    param_df = param_df,
    env = e
  )
}

all_models <- imap(workspace_files, load_workspace)

plot_all  <- bind_rows(map(all_models, "plot_df"))
param_all <- bind_rows(map(all_models, "param_df"))

ref_env <- all_models[[1]]$env

dates <- ref_env$dates
swap_cols <- ref_env$swap_cols
swap_maturities <- ref_env$swap_maturities
swap_obs_real <- ref_env$swap_obs_real

idx <- 3:length(dates)

# ============================================================
# PARAMETER TABLE
# ============================================================

param_table <- param_all %>%
  select(model, parameter, estimate) %>%
  pivot_wider(
    names_from = model,
    values_from = estimate
  )

print(param_table)

write.csv(
  param_table,
  "~/Desktop/parameter_table_all_models.csv",
  row.names = FALSE
)

# ============================================================
# MODEL FIT TABLE
# ============================================================

model_fit_table <- param_all %>%
  distinct(model, objective_value, best_start) %>%
  arrange(objective_value)

print(model_fit_table)

write.csv(
  model_fit_table,
  "~/Desktop/model_fit_table_all_models.csv",
  row.names = FALSE
)

# ============================================================
# AIC AND BIC TABLE
# ============================================================

aic_bic_table <- map_dfr(all_models, function(mod) {
  
  e <- mod$env
  
  loglik <- e$ekf_out_real$loglik
  k <- length(e$ekf_out_real$pars)
  n <- nrow(e$Z_obs)
  
  tibble(
    model = unique(mod$param_df$model),
    loglik = loglik,
    n_obs = n,
    n_parameters = k,
    AIC = -2 * loglik + 2 * k,
    BIC = -2 * loglik + log(n) * k
  )
})

aic_bic_table <- aic_bic_table %>%
  arrange(AIC)

print(aic_bic_table)

write.csv(
  aic_bic_table,
  "~/Desktop/aic_bic_table_all_models.csv",
  row.names = FALSE
)
# ============================================================
# OVERALL MSE / RMSE TABLE
# ============================================================

mse_model_table <- plot_all %>%
  filter(type == "resid") %>%
  group_by(model) %>%
  summarise(
    MSE = mean(value^2, na.rm = TRUE),
    RMSE = sqrt(MSE),
    .groups = "drop"
  )

print(mse_model_table)


# ============================================================
# RMSE TABLE
# ============================================================

rmse_table <- plot_all %>%
  filter(type == "resid") %>%
  group_by(model, maturity) %>%
  summarise(
    RMSE = sqrt(mean(value^2, na.rm = TRUE)),
    Mean_residual = mean(value, na.rm = TRUE),
    SD_residual = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

print(rmse_table, n = Inf)

write.csv(
  rmse_table,
  "~/Desktop/rmse_table_all_models.csv",
  row.names = FALSE
)

# ============================================================
# RMSE PLOT ALL MODELS
# ============================================================

rmse_plot_all <- rmse_table %>%
  mutate(
    maturity_num = as.numeric(gsub("Y", "", maturity))
  ) %>%
  arrange(model, maturity_num)

model_cols <- c(
  "Optimization 1" = "#F05A4A",
  "Initial model"  = "orange",
  "Optimization 2" = "#F04BE3",
  "Optimization 3" = "#19B52D",
  "Optimization 4" = "#4A7BFF",
  "Optimization 5" = "#27C4C9"
)

p_rmse <- ggplot(
  rmse_plot_all,
  aes(
    x = maturity_num,
    y = RMSE,
    color = model,
    group = model
  )
) +
  geom_line(linewidth = 1, linetype = "dashed") +
  geom_point(size = 2.5) +
  
  scale_color_manual(values = model_cols) +
  
  scale_x_continuous(
    breaks = c(1, 2, 3, 5, 7, 10, 15, 20, 30)
  ) +
  
  labs(
    x = "Maturity (Years)",
    y = "RMSE",
    color = "Model"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.position = "bottom"
  )

print(p_rmse)

print(p_rmse)

########
p_rmse <- ggplot(
  rmse_plot_all,
  aes(
    x = maturity_num,
    y = RMSE,
    color = model,
    group = model
  )
) +
  geom_line(linewidth = 1, linetype = "dashed") +
  geom_point(size = 2.5) +
  scale_x_continuous(
    breaks = c(1, 2, 3, 5, 7, 10, 15, 20, 30)
  ) +
  labs(
    x = "Maturity (Years)",
    y = "RMSE",
    color = "Model"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.position = "bottom"
  )

print(p_rmse)

# ============================================================
# OBSERVED VS PREDICTED SWAPS
# ============================================================
model_names <- names(workspace_files)

cols <- c(
  "black",     # Observed
  "#F05A4A",   # Optimization 1
  "orange",    # Initial model
  "#F04BE3",   # Optimization 2
  "#19B52D",   # Optimization 3
  "#4A7BFF",   # Optimization 4
  "#27C4C9"    # Optimization 5
)

op <- par(no.readonly = TRUE)

par(
  mfrow = c(3, 3),
  mar = c(1.8, 1.8, 1.2, 0.2),
  oma = c(7, 3.5, 1, 0),
  mgp = c(1.4, 0.4, 0),
  tcl = -0.25
)

year_ticks <- seq(
  from = as.Date(format(min(dates, na.rm = TRUE), "%Y-01-01")),
  to   = as.Date(format(max(dates, na.rm = TRUE), "%Y-01-01")),
  by   = "2 years"
)

for (j in seq_along(swap_cols)) {
  
  y_min <- min(
    swap_obs_real[idx, j],
    sapply(all_models, function(x) {
      min(x$env$swap_pred_real[idx, j], na.rm = TRUE)
    }),
    na.rm = TRUE
  )
  
  y_max <- max(
    swap_obs_real[idx, j],
    sapply(all_models, function(x) {
      max(x$env$swap_pred_real[idx, j], na.rm = TRUE)
    }),
    na.rm = TRUE
  )
  
  plot(
    dates[idx],
    swap_obs_real[idx, j],
    type = "l",
    main = paste0(swap_maturities[j], "Y"),
    xlab = "",
    ylab = "",
    lwd = 1.8,
    col = "black",
    ylim = c(y_min, y_max),
    xaxt = "n",
    cex.main = 1.45,
    cex.axis = 1.4
  )
  
  axis(
    1,
    at = year_ticks,
    labels = format(year_ticks, "%Y"),
    cex.axis = 1.4
  )
  
  abline(v = year_ticks, col = "grey85", lwd = 0.8)
  abline(h = axTicks(2), col = "grey85", lwd = 0.8)
  
  lines(
    dates[idx],
    swap_obs_real[idx, j],
    lwd = 1.8,
    col = "black"
  )
  
  for (m in seq_along(all_models)) {
    lines(
      dates[idx],
      all_models[[m]]$env$swap_pred_real[idx, j],
      lwd = 1.5,
      lty = 1,
      col = cols[m + 1]
    )
  }
}

mtext("Date", side = 1, outer = TRUE, line = 1.1, cex = 1.2)
mtext("Swap rate", side = 2, outer = TRUE, line = 0.8, cex = 1.2)

par(
  fig = c(0, 1, 0, 1),
  mar = c(0, 0, 0, 0),
  oma = c(0, 0, 0, 0),
  new = TRUE
)

plot.new()

legend(
  x = "bottom",
  inset = 0.02,
  horiz = TRUE,
  bty = "n",
  cex = 1.25,
  xpd = NA,
  legend = c("Observed", model_names),
  col = cols[seq_len(length(model_names) + 1)],
  lty = 1,
  lwd = c(1.8, rep(1.5, length(model_names)))
)

par(op)
#############
library(tidyverse)
library(patchwork)

plot_fit_terms_gg <- function(selected_terms, plot_title = "") {
  
  model_cols <- c(
    "Observed"       = "black",
    "Initial model"  = "orange",
    "Optimization 1" = "#F05A4A",
    "Optimization 2" = "#F04BE3",
    "Optimization 3" = "#19B52D",
    "Optimization 4" = "#4A7BFF",
    "Optimization 5" = "#27C4C9"
  )
  
  plot_df <- plot_all %>%
    filter(type %in% c("obs", "pred")) %>%
    mutate(
      Term_num = as.numeric(gsub("[^0-9.]", "", maturity)),
      Term = paste0(Term_num, "Y"),
      Model = if_else(type == "obs", "Observed", model)
    ) %>%
    filter(
      Term_num %in% selected_terms,
      Model %in% names(model_cols)
    ) %>%
    distinct(Date, Term_num, Term, Model, value) %>%
    mutate(
      Term = factor(
        Term,
        levels = paste0(selected_terms, "Y")
      ),
      Model = factor(
        Model,
        levels = names(model_cols)
      )
    )
  
  ggplot(
    plot_df,
    aes(
      x = Date,
      y = value,
      color = Model,
      linetype = Model
    )
  ) +
    
    geom_line(linewidth = 0.8) +
    
    facet_wrap(
      ~ Term,
      ncol = 2,
      scales = "free_y"
    ) +
    
    scale_color_manual(
      values = model_cols
    ) +
    
    scale_linetype_manual(
      values = rep("solid", length(model_cols))
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
    
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        keywidth = unit(1.3, "cm")
      ),
      
      linetype = "none"
    ) +
    
    theme_bw(base_size = 14) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 18,
        hjust = 0.5
      ),
      
      strip.background = element_blank(),
      
      strip.text = element_text(
        face = "bold",
        size = 10
      ),
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 12
      ),
      
      axis.text.y = element_text(size = 12),
      
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      
      legend.position = c(0.464, -0.14),
      legend.direction = "horizontal",
      legend.box = "horizontal",
      
      legend.text = element_text(size = 9.5),
      
      legend.key.width = unit(1.3, "cm"),
      legend.spacing.x = unit(0.25, "cm"),
      
      panel.grid.major = element_line(
        color = "grey85",
        linewidth = 0.4
      ),
      
      panel.grid.minor = element_blank(),
      
      plot.margin = margin(
        10,
        10,
        50,
        10
      )
    )
}

plot_fit_terms_gg(
  selected_terms = c(1, 3, 5, 10, 20, 30)
)






# ============================================================
# RESIDUALS
# ============================================================

op <- par(no.readonly = TRUE)

par(
  mfrow = c(3, 3),
  mar = c(1.8, 1.8, 1.2, 0.2),
  oma = c(7, 3.5, 1, 0),
  mgp = c(1.4, 0.4, 0),
  tcl = -0.25
)

for (j in seq_along(swap_cols)) {
  
  y_min <- min(
    sapply(all_models, function(x) {
      min(x$env$resid_swap_pred_real[idx, j], na.rm = TRUE)
    }),
    na.rm = TRUE
  )
  
  y_max <- max(
    sapply(all_models, function(x) {
      max(x$env$resid_swap_pred_real[idx, j], na.rm = TRUE)
    }),
    na.rm = TRUE
  )
  
  plot(
    dates[idx],
    all_models[[1]]$env$resid_swap_pred_real[idx, j],
    type = "l",
    main = paste0(swap_maturities[j], "Y"),
    xlab = "",
    ylab = "",
    lwd = 1.5,
    col = cols[2],
    ylim = c(y_min, y_max),
    xaxt = "n",
    cex.main = 1.45,
    cex.axis = 1.4
  )
  
  axis(
    1,
    at = year_ticks,
    labels = format(year_ticks, "%Y"),
    cex.axis = 1.4
  )
  
  abline(v = year_ticks, col = "grey85", lwd = 0.8)
  abline(h = axTicks(2), col = "grey85", lwd = 0.8)
  abline(h = 0, col = "black", lty = 2, lwd = 1)
  
  for (m in seq_along(all_models)) {
    lines(
      dates[idx],
      all_models[[m]]$env$resid_swap_pred_real[idx, j],
      lwd = 1.5,
      lty = 1,
      col = cols[m + 1]
    )
  }
}

mtext("Date", side = 1, outer = TRUE, line = 1.1, cex = 1.2)
mtext("Residual", side = 2, outer = TRUE, line = 0.8, cex = 1.2)

par(
  fig = c(0, 1, 0, 1),
  mar = c(0, 0, 0, 0),
  oma = c(0, 0, 0, 0),
  new = TRUE
)

plot.new()

legend(
  x = "bottom",
  inset = 0.02,
  horiz = TRUE,
  bty = "n",
  cex = 1.25,
  seg.len = 3,
  xpd = NA,
  legend = model_names,
  col = cols[2:(length(model_names) + 1)],
  lty = 1,
  lwd = 2.2
)

par(op)



############
library(tidyverse)
library(patchwork)

plot_residual_terms_gg <- function(selected_terms, plot_title = "") {
  
  model_cols <- c(
    "Initial model"  = "orange",
    "Optimization 1" = "#F05A4A",
    "Optimization 2" = "#F04BE3",
    "Optimization 3" = "#19B52D",
    "Optimization 4" = "#4A7BFF",
    "Optimization 5" = "#27C4C9"
  )
  
  plot_df <- plot_all %>%
    filter(type == "resid") %>%
    mutate(
      Term_num = as.numeric(gsub("[^0-9.]", "", maturity)),
      Term = paste0(Term_num, "Y"),
      Model = model
    ) %>%
    filter(
      Term_num %in% selected_terms,
      Model %in% names(model_cols)
    ) %>%
    distinct(Date, Term_num, Term, Model, value) %>%
    mutate(
      Term = factor(
        Term,
        levels = paste0(selected_terms, "Y")
      ),
      Model = factor(
        Model,
        levels = names(model_cols)
      )
    )
  
  ggplot(
    plot_df,
    aes(
      x = Date,
      y = value,
      color = Model,
      linetype = Model
    )
  ) +
    
    geom_hline(
      yintercept = 0,
      color = "black",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    geom_line(linewidth = 0.8) +
    
    facet_wrap(
      ~ Term,
      ncol = 2,
      scales = "free_y"
    ) +
    
    scale_color_manual(
      values = model_cols
    ) +
    
    scale_linetype_manual(
      values = rep("solid", length(model_cols))
    ) +
    
    scale_x_date(
      date_breaks = "2 years",
      date_labels = "%Y"
    ) +
    
    labs(
      title = plot_title,
      x = "Date",
      y = "Residual",
      color = NULL,
      linetype = NULL
    ) +
    
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        keywidth = unit(1.3, "cm")
      ),
      
      linetype = "none"
    ) +
    
    theme_bw(base_size = 14) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 18,
        hjust = 0.5
      ),
      
      strip.background = element_blank(),
      
      strip.text = element_text(
        face = "bold",
        size = 10
      ),
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 12
      ),
      
      axis.text.y = element_text(size = 12),
      
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      
      legend.text = element_text(size = 10),
      
      legend.key.width = unit(1.3, "cm"),
      legend.spacing.x = unit(0.25, "cm"),
      
      panel.grid.major = element_line(
        color = "grey85",
        linewidth = 0.4
      ),
      
      panel.grid.minor = element_blank(),
      
      plot.margin = margin(
        10,
        10,
        10,
        10
      )
    )
}

plot_residual_terms_gg(
  selected_terms = c(1, 3, 5, 10, 20, 30)
)


# ============================================================
# STATE PLOTS
# ============================================================

state_df <- map_dfr(all_models, function(mod) {
  
  e <- mod$env
  model_name <- unique(mod$param_df$model)
  
  tibble(
    Date = e$dates,
    model = model_name,
    X_predicted = e$ekf_out_real$x_pred[, 1],
    Y_predicted = e$ekf_out_real$x_pred[, 2]
  )
}) %>%
  slice(3:n())

model_cols <- c(
  "Optimization 1" = "orange",
  "Initial model"  = "#F05A4A",
  "Optimization 2" = "#F04BE3",
  "Optimization 3" = "#19B52D",
  "Optimization 4" = "#4A7BFF",
  "Optimization 5" = "#27C4C9"
)

# ============================================================
# X_t STATE
# ============================================================

p_state_X <- ggplot(
  state_df,
  aes(
    x = Date,
    y = X_predicted,
    color = model,
    group = model
  )
) +
  geom_line(linewidth = 1) +
  
  scale_color_manual(values = model_cols) +
  
  labs(
    x = "Date",
    y = expression(X[t]),
    color = "Model"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.position = "bottom"
  )

print(p_state_X)

# ============================================================
# Y_t STATE
# ============================================================

p_state_Y <- ggplot(
  state_df,
  aes(
    x = Date,
    y = Y_predicted,
    color = model,
    group = model
  )
) +
  geom_line(linewidth = 1) +
  
  scale_color_manual(values = model_cols) +
  
  labs(
    x = "Date",
    y = expression(Y[t]),
    color = "Model"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.position = "bottom"
  )

print(p_state_Y)





combined_states_plot <- p_state_X + p_state_Y +
  plot_layout(ncol = 2, guides = "collect") &
  theme(
    legend.position = "bottom"
  )

print(combined_states_plot)




# ============================================================
# DONE
# ============================================================


# ============================================================
# COMPARE 2 EKF + 2 UKF MODELS
# ============================================================

load_workspace_flexible <- function(file, model_name) {
  
  e <- new.env()
  load(file, envir = e)
  
  # UKF compatibility:
  # If the workspace has ukf_out_real instead of ekf_out_real,
  # copy it to ekf_out_real so the old comparison code still works.
  if (!"ekf_out_real" %in% ls(e) && "ukf_out_real" %in% ls(e)) {
    e$ekf_out_real <- e$ukf_out_real
  }
  
  needed <- c(
    "dates",
    "swap_cols",
    "swap_maturities",
    "swap_obs_real",
    "swap_pred_real",
    "resid_swap_pred_real",
    "ekf_out_real",
    "Z_obs"
  )
  
  missing <- setdiff(needed, ls(e))
  
  if (length(missing) > 0) {
    stop(
      "Mangler i ", model_name, ": ",
      paste(missing, collapse = ", ")
    )
  }
  
  swap_df  <- as.data.frame(e$swap_obs_real)
  pred_df  <- as.data.frame(e$swap_pred_real)
  resid_df <- as.data.frame(e$resid_swap_pred_real)
  
  colnames(swap_df)  <- e$swap_cols
  colnames(pred_df)  <- e$swap_cols
  colnames(resid_df) <- e$swap_cols
  
  plot_df <- tibble(Date = e$dates) %>%
    bind_cols(
      swap_df  %>% rename_with(~ paste0("obs_", .x)),
      pred_df  %>% rename_with(~ paste0("pred_", .x)),
      resid_df %>% rename_with(~ paste0("resid_", .x))
    ) %>%
    slice(3:n()) %>%
    pivot_longer(
      cols = -Date,
      names_to = c("type", "maturity"),
      names_pattern = "(obs|pred|resid)_(.*)",
      values_to = "value"
    ) %>%
    mutate(model = model_name)
  
  param_df <- tibble(
    model = model_name,
    parameter = names(e$ekf_out_real$pars),
    estimate = as.numeric(e$ekf_out_real$pars),
    objective_value = if ("best_fit" %in% ls(e)) e$best_fit$value else NA_real_,
    best_start = if ("best_idx" %in% ls(e)) e$best_idx else NA_integer_
  )
  
  list(
    plot_df = plot_df,
    param_df = param_df,
    env = e
  )
}

workspace_files_4 <- c(
  "Optimization 4 EKF" =
    "~/Desktop/ekf_workspace_tilfoej_a_uden_udevidet.RData",
  
  "Optimization 5 EKF" =
    "~/Desktop/ekf_workspace_tilfoej_a_med_andet_udvidet_int.RData",
  
  "Optimization 4 UKF" =
    "~/Desktop/ukf_workspace_tilfoej_a_uden_udvidet_int.RData",
  
  "Optimization 5 UKF" =
    "~/Desktop/ukf_workspace_tilfoej_a_med_andet_udvidet_int.RData"
)

all_models_4 <- imap(workspace_files_4, load_workspace_flexible)

plot_all_4  <- bind_rows(map(all_models_4, "plot_df"))
param_all_4 <- bind_rows(map(all_models_4, "param_df"))

# ============================================================
# TABLES
# ============================================================

param_table_4 <- param_all_4 %>%
  select(model, parameter, estimate) %>%
  pivot_wider(
    names_from = model,
    values_from = estimate
  )

print(param_table_4)

write.csv(
  param_table_4,
  "~/Desktop/parameter_table_EKF_UKF_4_models.csv",
  row.names = FALSE
)

model_fit_table_4 <- param_all_4 %>%
  distinct(model, objective_value, best_start) %>%
  arrange(objective_value)

print(model_fit_table_4)

write.csv(
  model_fit_table_4,
  "~/Desktop/model_fit_table_EKF_UKF_4_models.csv",
  row.names = FALSE
)

aic_bic_table_4 <- map_dfr(all_models_4, function(mod) {
  
  e <- mod$env
  
  loglik <- e$ekf_out_real$loglik
  k <- length(e$ekf_out_real$pars)
  n <- nrow(e$Z_obs)
  
  tibble(
    model = unique(mod$param_df$model),
    loglik = loglik,
    n_obs = n,
    n_parameters = k,
    AIC = -2 * loglik + 2 * k,
    BIC = -2 * loglik + log(n) * k
  )
}) %>%
  arrange(AIC)

print(aic_bic_table_4)

write.csv(
  aic_bic_table_4,
  "~/Desktop/aic_bic_table_EKF_UKF_4_models.csv",
  row.names = FALSE
)

# ============================================================
# RMSE
# ============================================================

mse_model_table_4 <- plot_all_4 %>%
  filter(type == "resid") %>%
  group_by(model) %>%
  summarise(
    MSE = mean(value^2, na.rm = TRUE),
    RMSE = sqrt(MSE),
    .groups = "drop"
  ) %>%
  arrange(RMSE)

print(mse_model_table_4)

write.csv(
  mse_model_table_4,
  "~/Desktop/overall_rmse_EKF_UKF_4_models.csv",
  row.names = FALSE
)

rmse_table_4 <- plot_all_4 %>%
  filter(type == "resid") %>%
  group_by(model, maturity) %>%
  summarise(
    RMSE = sqrt(mean(value^2, na.rm = TRUE)),
    Mean_residual = mean(value, na.rm = TRUE),
    SD_residual = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

print(rmse_table_4, n = Inf)

write.csv(
  rmse_table_4,
  "~/Desktop/rmse_by_maturity_EKF_UKF_4_models.csv",
  row.names = FALSE
)

# ============================================================
# RMSE PLOT
# ============================================================

model_cols_4 <- c(
  "Optimization 4 EKF" = "#4A7BFF",
  "Optimization 5 EKF" = "#27C4C9",
  "Optimization 4 UKF" = "#F05A4A",
  "Optimization 5 UKF" = "#19B52D"
)

rmse_plot_4 <- rmse_table_4 %>%
  mutate(
    maturity_num = as.numeric(gsub("Y", "", maturity))
  ) %>%
  arrange(model, maturity_num)

p_rmse_4 <- ggplot(
  rmse_plot_4,
  aes(
    x = maturity_num,
    y = RMSE,
    color = model,
    group = model
  )
) +
  geom_line(linewidth = 1, linetype = "dashed") +
  geom_point(size = 2.5) +
  scale_color_manual(values = model_cols_4) +
  scale_x_continuous(
    breaks = c(1, 2, 3, 5, 7, 10, 15, 20, 30)
  ) +
  labs(
    x = "Maturity (Years)",
    y = "RMSE",
    color = "Model"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.position = "bottom"
  )

print(p_rmse_4)

# ============================================================
# STATE PLOTS: X_t AND Y_t
# ============================================================

state_df_4 <- map_dfr(all_models_4, function(mod) {
  
  e <- mod$env
  model_name <- unique(mod$param_df$model)
  
  tibble(
    Date = e$dates,
    model = model_name,
    X_predicted = e$ekf_out_real$x_pred[, 1],
    Y_predicted = e$ekf_out_real$x_pred[, 2]
  )
}) %>%
  slice(3:n())
ylim_shared <- c(-1.3, 1.5)

p_state_X_4 <- ggplot(
  state_df_4,
  aes(
    x = Date,
    y = X_predicted,
    color = model,
    group = model
  )
) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = model_cols_4) +
  coord_cartesian(ylim = ylim_shared) +
  labs(
    x = "Date",
    y = expression(X[t]),
    color = "Model"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.position = "bottom"
  )

print(p_state_X_4)

p_state_Y_4 <- ggplot(
  state_df_4,
  aes(
    x = Date,
    y = Y_predicted,
    color = model,
    group = model
  )
) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = model_cols_4) +
  coord_cartesian(ylim = ylim_shared) +
  labs(
    x = "Date",
    y = expression(Y[t]),
    color = "Model"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.position = "bottom"
  )

print(p_state_Y_4)

combined_states_plot_4 <- p_state_X_4 + p_state_Y_4 +
  plot_layout(ncol = 2, guides = "collect") &
  theme(
    legend.position = "bottom"
  )

print(combined_states_plot_4)

# ============================================================
# SWAP PLOTS
# ============================================================

ref_env_4 <- all_models_4[[1]]$env

dates_4 <- ref_env_4$dates
swap_cols_4 <- ref_env_4$swap_cols
swap_maturities_4 <- ref_env_4$swap_maturities
swap_obs_real_4 <- ref_env_4$swap_obs_real

idx_4 <- 3:length(dates_4)

model_names_4 <- names(workspace_files_4)

cols_4 <- c(
  "black",
  "#4A7BFF",
  "#27C4C9",
  "#F05A4A",
  "#19B52D"
)

op <- par(no.readonly = TRUE)

par(
  mfrow = c(3, 3),
  mar = c(1.8, 1.8, 1.2, 0.2),
  oma = c(7, 3.5, 1, 0),
  mgp = c(1.4, 0.4, 0),
  tcl = -0.25
)

year_ticks_4 <- seq(
  from = as.Date(format(min(dates_4, na.rm = TRUE), "%Y-01-01")),
  to   = as.Date(format(max(dates_4, na.rm = TRUE), "%Y-01-01")),
  by   = "2 years"
)

for (j in seq_along(swap_cols_4)) {
  
  y_min <- min(
    swap_obs_real_4[idx_4, j],
    sapply(all_models_4, function(x) {
      min(x$env$swap_pred_real[idx_4, j], na.rm = TRUE)
    }),
    na.rm = TRUE
  )
  
  y_max <- max(
    swap_obs_real_4[idx_4, j],
    sapply(all_models_4, function(x) {
      max(x$env$swap_pred_real[idx_4, j], na.rm = TRUE)
    }),
    na.rm = TRUE
  )
  
  plot(
    dates_4[idx_4],
    swap_obs_real_4[idx_4, j],
    type = "l",
    main = paste0(swap_maturities_4[j], "Y"),
    xlab = "",
    ylab = "",
    lwd = 1.8,
    col = "black",
    ylim = c(y_min, y_max),
    xaxt = "n",
    cex.main = 1.45,
    cex.axis = 1.4
  )
  
  axis(
    1,
    at = year_ticks_4,
    labels = format(year_ticks_4, "%Y"),
    cex.axis = 1.4
  )
  
  abline(v = year_ticks_4, col = "grey85", lwd = 0.8)
  abline(h = axTicks(2), col = "grey85", lwd = 0.8)
  
  lines(
    dates_4[idx_4],
    swap_obs_real_4[idx_4, j],
    lwd = 1.8,
    col = "black"
  )
  
  for (m in seq_along(all_models_4)) {
    lines(
      dates_4[idx_4],
      all_models_4[[m]]$env$swap_pred_real[idx_4, j],
      lwd = 1.5,
      col = cols_4[m + 1]
    )
  }
}

mtext("Date", side = 1, outer = TRUE, line = 1.1, cex = 1.2)
mtext("Swap rate", side = 2, outer = TRUE, line = 0.8, cex = 1.2)

par(
  fig = c(0, 1, 0, 1),
  mar = c(0, 0, 0, 0),
  oma = c(0, 0, 0, 0),
  new = TRUE
)

plot.new()

legend(
  x = "bottom",
  inset = 0.02,
  horiz = TRUE,
  bty = "n",
  cex = 1.15,
  xpd = NA,
  legend = c("Observed", model_names_4),
  col = cols_4,
  lty = 1,
  lwd = c(1.8, rep(1.5, length(model_names_4)))
)

par(op)


# ============================================================
# RESIDUAL PLOTS
# ============================================================

op <- par(no.readonly = TRUE)

par(
  mfrow = c(3, 3),
  mar = c(1.8, 1.8, 1.2, 0.2),
  oma = c(7, 3.5, 1, 0),
  mgp = c(1.4, 0.4, 0),
  tcl = -0.25
)

for (j in seq_along(swap_cols_4)) {
  
  y_min <- min(
    sapply(all_models_4, function(x) {
      min(x$env$resid_swap_pred_real[idx_4, j], na.rm = TRUE)
    }),
    na.rm = TRUE
  )
  
  y_max <- max(
    sapply(all_models_4, function(x) {
      max(x$env$resid_swap_pred_real[idx_4, j], na.rm = TRUE)
    }),
    na.rm = TRUE
  )
  
  plot(
    dates_4[idx_4],
    all_models_4[[1]]$env$resid_swap_pred_real[idx_4, j],
    type = "l",
    main = paste0(swap_maturities_4[j], "Y"),
    xlab = "",
    ylab = "",
    lwd = 1.5,
    col = cols_4[2],
    ylim = c(y_min, y_max),
    xaxt = "n",
    cex.main = 1.45,
    cex.axis = 1.4
  )
  
  axis(
    1,
    at = year_ticks_4,
    labels = format(year_ticks_4, "%Y"),
    cex.axis = 1.4
  )
  
  abline(v = year_ticks_4, col = "grey85", lwd = 0.8)
  abline(h = axTicks(2), col = "grey85", lwd = 0.8)
  abline(h = 0, col = "black", lty = 2, lwd = 1)
  
  for (m in seq_along(all_models_4)) {
    lines(
      dates_4[idx_4],
      all_models_4[[m]]$env$resid_swap_pred_real[idx_4, j],
      lwd = 1.5,
      col = cols_4[m + 1]
    )
  }
}

mtext("Date", side = 1, outer = TRUE, line = 1.1, cex = 1.2)
mtext("Residual", side = 2, outer = TRUE, line = 0.8, cex = 1.2)

par(
  fig = c(0, 1, 0, 1),
  mar = c(0, 0, 0, 0),
  oma = c(0, 0, 0, 0),
  new = TRUE
)

plot.new()

legend(
  x = "bottom",
  inset = 0.02,
  horiz = TRUE,
  bty = "n",
  cex = 1.15,
  seg.len = 3,
  xpd = NA,
  legend = model_names_4,
  col = cols_4[2:5],
  lty = 1,
  lwd = 2.2
)

par(op)




# ============================================================
# RESIDUAL PLOTS
# ============================================================

op <- par(no.readonly = TRUE)

# Global y-axis limits across all maturities and all models
global_y_min <- min(
  sapply(all_models_4, function(x) {
    min(x$env$resid_swap_pred_real[idx_4, ], na.rm = TRUE)
  }),
  na.rm = TRUE
)

global_y_max <- max(
  sapply(all_models_4, function(x) {
    max(x$env$resid_swap_pred_real[idx_4, ], na.rm = TRUE)
  }),
  na.rm = TRUE
)

par(
  mfrow = c(3, 3),
  mar = c(1.8, 1.8, 1.2, 0.2),
  oma = c(7, 3.5, 1, 0),
  mgp = c(1.4, 0.4, 0),
  tcl = -0.25
)

for (j in seq_along(swap_cols_4)) {
  
  plot(
    dates_4[idx_4],
    all_models_4[[1]]$env$resid_swap_pred_real[idx_4, j],
    type = "l",
    main = paste0(swap_maturities_4[j], "Y"),
    xlab = "",
    ylab = "",
    lwd = 1.5,
    col = cols_4[2],
    ylim = c(global_y_min, global_y_max),
    xaxt = "n",
    cex.main = 1.45,
    cex.axis = 1.4
  )
  
  axis(
    1,
    at = year_ticks_4,
    labels = format(year_ticks_4, "%Y"),
    cex.axis = 1.4
  )
  
  abline(v = year_ticks_4, col = "grey85", lwd = 0.8)
  abline(h = axTicks(2), col = "grey85", lwd = 0.8)
  abline(h = 0, col = "black", lty = 2, lwd = 1)
  
  for (m in seq_along(all_models_4)) {
    lines(
      dates_4[idx_4],
      all_models_4[[m]]$env$resid_swap_pred_real[idx_4, j],
      lwd = 1.5,
      col = cols_4[m + 1]
    )
  }
}

mtext("Date", side = 1, outer = TRUE, line = 1.1, cex = 1.2)
mtext("Residual", side = 2, outer = TRUE, line = 0.8, cex = 1.2)

par(
  fig = c(0, 1, 0, 1),
  mar = c(0, 0, 0, 0),
  oma = c(0, 0, 0, 0),
  new = TRUE
)

plot.new()

legend(
  x = "bottom",
  inset = 0.02,
  horiz = TRUE,
  bty = "n",
  cex = 1.15,
  seg.len = 3,
  xpd = NA,
  legend = model_names_4,
  col = cols_4[2:5],
  lty = 1,
  lwd = 2.2
)

par(op)




############## appandix

library(purrr)

workspace_files_final <- c(
  "Initial model" =
    "~/Desktop/ekf_workspace_first_model.RData",
  
  "Optimization 1" =
    "~/Desktop/ekf_workspace_correlated_model.RData",
  
  "Optimization 2" =
    "~/Desktop/ekf_workspace_udvidet_int_model.RData",
  
  "Optimization 3" =
    "~/Desktop/ekf_workspace_tilfoej_a.RData",
  
  "Optimization 4" =
    "~/Desktop/ekf_workspace_tilfoej_a_uden_udevidet.RData",
  
  "Optimization 5" =
    "~/Desktop/ekf_workspace_tilfoej_a_med_andet_udvidet_int.RData",
  
  "Optimization 4 UKF" =
    "~/Desktop/ukf_workspace_tilfoej_a_uden_udvidet_int.RData",
  
  "Optimization 5 UKF" =
    "~/Desktop/ukf_workspace_tilfoej_a_med_andet_udvidet_int.RData"
)

load_workspace_final <- function(file, model_name) {
  
  e <- new.env()
  load(file, envir = e)
  
  if (!"ekf_out_real" %in% ls(e) && "ukf_out_real" %in% ls(e)) {
    e$ekf_out_real <- e$ukf_out_real
  }
  
  needed <- c(
    "dates",
    "swap_maturities",
    "swap_obs_real",
    "swap_pred_real",
    "resid_swap_pred_real"
  )
  
  missing <- setdiff(needed, ls(e))
  
  if (length(missing) > 0) {
    stop(
      "Mangler i ", model_name, ": ",
      paste(missing, collapse = ", ")
    )
  }
  
  list(
    model = model_name,
    dates = e$dates,
    taus = e$swap_maturities,
    obs = e$swap_obs_real,
    pred = e$swap_pred_real,
    resid = e$resid_swap_pred_real
  )
}

workspace_models <- imap(workspace_files_final, load_workspace_final)

# ------------------------------------------------------------
# Add baseline models
# ------------------------------------------------------------

baseline_models <- list(
  "Vasicek" = list(
    model = "Vasicek",
    dates = dates,
    taus = taus,
    obs = Y_real,
    pred = Y_pred_vas_real,
    resid = resid_pred_vas_real
  ),
  
  "CIR" = list(
    model = "CIR",
    dates = dates,
    taus = taus,
    obs = Y_real,
    pred = Y_pred_cir_real,
    resid = resid_pred_cir_real
  )
)

all_models_final <- c(workspace_models, baseline_models)

ref <- all_models_final[[1]]

dates_final <- ref$dates
taus_final <- as.numeric(ref$taus)
obs_final <- ref$obs
idx_final <- 3:length(dates_final)

model_names_final <- names(all_models_final)

cols_final <- c(
  "Observed" = "black",
  "Initial model" = "orange",
  "Optimization 1" = "#F05A4A",
  "Optimization 2" = "#F04BE3",
  "Optimization 3" = "#19B52D",
  "Optimization 4" = "#4A7BFF",
  "Optimization 5" = "#27C4C9",
  "Optimization 4 UKF" = "#B22222",
  "Optimization 5 UKF" = "#228B22",
  "Vasicek" = "#ff69b4",
  "CIR" = "#1f77b4"
)

year_ticks_final <- seq(
  from = as.Date(format(min(dates_final, na.rm = TRUE), "%Y-01-01")),
  to   = as.Date(format(max(dates_final, na.rm = TRUE), "%Y-01-01")),
  by   = "2 years"
)

# ============================================================
# SWAP PLOT
# ============================================================

op <- par(no.readonly = TRUE)

par(
  mfrow = c(3, 3),
  mar = c(1.8, 1.8, 1.2, 0.2),
  oma = c(7, 3.5, 1, 0),
  mgp = c(1.4, 0.4, 0),
  tcl = -0.25
)

for (j in seq_along(taus_final)) {
  
  y_min <- min(
    obs_final[idx_final, j],
    sapply(all_models_final, function(x) {
      min(x$pred[idx_final, j], na.rm = TRUE)
    }),
    na.rm = TRUE
  )
  
  y_max <- max(
    obs_final[idx_final, j],
    sapply(all_models_final, function(x) {
      max(x$pred[idx_final, j], na.rm = TRUE)
    }),
    na.rm = TRUE
  )
  
  plot(
    dates_final[idx_final],
    obs_final[idx_final, j],
    type = "l",
    main = paste0(taus_final[j], "Y"),
    xlab = "",
    ylab = "",
    lwd = 1.8,
    col = "black",
    ylim = c(y_min, y_max),
    xaxt = "n",
    cex.main = 1.45,
    cex.axis = 1.4
  )
  
  axis(
    1,
    at = year_ticks_final,
    labels = format(year_ticks_final, "%Y"),
    cex.axis = 1.4
  )
  
  abline(v = year_ticks_final, col = "grey85", lwd = 0.8)
  abline(h = axTicks(2), col = "grey85", lwd = 0.8)
  
  lines(dates_final[idx_final], obs_final[idx_final, j], lwd = 1.8, col = "black")
  
  for (model_name in model_names_final) {
    lines(
      dates_final[idx_final],
      all_models_final[[model_name]]$pred[idx_final, j],
      lwd = 1.5,
      col = cols_final[model_name]
    )
  }
}

mtext("Date", side = 1, outer = TRUE, line = 1.1, cex = 1.2)
mtext("Swap rate", side = 2, outer = TRUE, line = 0.8, cex = 1.2)

par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE)
plot.new()

legend(
  x = "bottom",
  inset = 0.02,
  horiz = FALSE,
  ncol = ceiling(length(c("Observed", model_names_final)) / 2),
  bty = "n",
  cex = 0.85,
  seg.len = 2.3,
  xpd = NA,
  legend = c("Observed", model_names_final),
  col = cols_final[c("Observed", model_names_final)],
  lty = 1,
  lwd = c(1.8, rep(1.5, length(model_names_final)))
)

par(op)

# ============================================================
# RESIDUAL PLOT
# ============================================================

op <- par(no.readonly = TRUE)

par(
  mfrow = c(3, 3),
  mar = c(1.8, 1.8, 1.2, 0.2),
  oma = c(7, 3.5, 1, 0),
  mgp = c(1.4, 0.4, 0),
  tcl = -0.25
)

for (j in seq_along(taus_final)) {
  
  y_min <- min(
    sapply(all_models_final, function(x) {
      min(x$resid[idx_final, j], na.rm = TRUE)
    }),
    na.rm = TRUE
  )
  
  y_max <- max(
    sapply(all_models_final, function(x) {
      max(x$resid[idx_final, j], na.rm = TRUE)
    }),
    na.rm = TRUE
  )
  
  first_model <- model_names_final[1]
  
  plot(
    dates_final[idx_final],
    all_models_final[[first_model]]$resid[idx_final, j],
    type = "l",
    main = paste0(taus_final[j], "Y"),
    xlab = "",
    ylab = "",
    lwd = 1.5,
    col = cols_final[first_model],
    ylim = c(y_min, y_max),
    xaxt = "n",
    cex.main = 1.45,
    cex.axis = 1.4
  )
  
  axis(
    1,
    at = year_ticks_final,
    labels = format(year_ticks_final, "%Y"),
    cex.axis = 1.4
  )
  
  abline(v = year_ticks_final, col = "grey85", lwd = 0.8)
  abline(h = axTicks(2), col = "grey85", lwd = 0.8)
  abline(h = 0, col = "black", lty = 2, lwd = 1)
  
  for (model_name in model_names_final) {
    lines(
      dates_final[idx_final],
      all_models_final[[model_name]]$resid[idx_final, j],
      lwd = 1.5,
      col = cols_final[model_name]
    )
  }
}

mtext("Date", side = 1, outer = TRUE, line = 1.1, cex = 1.2)
mtext("Residual", side = 2, outer = TRUE, line = 0.8, cex = 1.2)

par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE)
plot.new()

legend(
  x = "bottom",
  inset = 0.02,
  horiz = FALSE,
  ncol = ceiling(length(model_names_final) / 2),
  bty = "n",
  cex = 0.85,
  seg.len = 2.3,
  xpd = NA,
  legend = model_names_final,
  col = cols_final[model_names_final],
  lty = 1,
  lwd = 2
)

par(op)


dev.off()
# ============================================================
# RMSE OVER TERM
# ============================================================

rmse_final <- do.call(
  rbind,
  lapply(model_names_final, function(model_name) {
    
    resid_mat <- all_models_final[[model_name]]$resid
    
    data.frame(
      model = model_name,
      maturity = taus_final,
      RMSE = sqrt(colMeans(resid_mat[idx_final, , drop = FALSE]^2, na.rm = TRUE))
    )
  })
)

op <- par(no.readonly = TRUE)

plot(
  NA,
  xlim = range(taus_final),
  ylim = range(rmse_final$RMSE, na.rm = TRUE),
  xlab = "Maturity (Years)",
  ylab = "RMSE",
  xaxt = "n",
  cex.axis = 1.3,
  cex.lab = 1.4
)

axis(
  1,
  at = c(1, 2, 3, 5, 7, 10, 15, 20, 30),
  labels = c(1, 2, 3, 5, 7, 10, 15, 20, 30),
  cex.axis = 1.3
)

abline(h = axTicks(2), col = "grey85", lwd = 0.8)
abline(v = c(1, 2, 3, 5, 7, 10, 15, 20, 30), col = "grey85", lwd = 0.8)

for (model_name in model_names_final) {
  
  tmp <- rmse_final[rmse_final$model == model_name, ]
  
  lines(
    tmp$maturity,
    tmp$RMSE,
    col = cols_final[model_name],
    lwd = 2,
    lty = 2
  )
  
  points(
    tmp$maturity,
    tmp$RMSE,
    col = cols_final[model_name],
    pch = 16,
    cex = 1.2
  )
}

legend(
  "topright",
  legend = model_names_final,
  col = cols_final[model_names_final],
  lty = 2,
  lwd = 2,
  pch = 16,
  bty = "n",
  cex = 0.85
)

par(op)

print(rmse_final)