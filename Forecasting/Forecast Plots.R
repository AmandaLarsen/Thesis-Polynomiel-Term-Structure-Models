

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
  library(grid)
})

input_files <- c(
  initial_opt1 = "Forecasting til int og 1 opt.RData",
  ekf_opt4_opt5 = "forecasting.RData",
  ukf_opt1_opt2 = "Forcasting_UKF.RData",
  baseline = "baselinemodeler_forecasting.RData"
)

save_tables <- FALSE
output_dir <- "samlede_forecast_tables"

exclude_object_names <- c("rolling5y_opt4")
exclude_model_regex <- c("^EKF\\s*Opt\\.?\\s*4$")

maturity_order <- c("Y1", "Y2", "Y3", "Y5", "Y7", "Y10", "Y15", "Y20", "Y30")
maturity_breaks <- c(1, 2, 3, 5, 7, 10, 15, 20, 30)

model_cols <- c(
  "Observed"      = "black",
  "Initial Model" = "orange",
  "Opt. 1 EKF"    = "#F05A4A",
  "Opt. 2 EKF"    = "#F04BE3",
  "Opt. 3 EKF"    = "#434323",
  "Opt. 5 EKF"    = "#27C4C9",
  "Opt. 4 UKF"    = "#19B52D",
  "Opt. 5 UKF"    = "#8A2BE2",
  "Vasicek"       = "#ff69b4",
  "CIR"           = "#1f77b4"
)


load_env <- function(file) {
  env <- new.env(parent = emptyenv())
  load(file, envir = env)
  env
}

first_existing <- function(x, names_vec) {
  for (nm in names_vec) {
    if (!is.null(x[[nm]])) return(x[[nm]])
  }
  NULL
}

clean_model_label <- function(object_name, model_label = NULL, source_name = NULL) {
  raw <- if (!is.null(model_label) && length(model_label) == 1 && !is.na(model_label)) {
    as.character(model_label)
  } else {
    object_name
  }
  
  label <- raw %>%
    str_replace_all("rolling5y_", "") %>%
    str_replace_all("_", " ") %>%
    str_squish()
  
  source_low <- tolower(source_name)
  object_low <- tolower(object_name)
  label_low  <- tolower(label)
  
  if (str_detect(source_low, "baseline")) {
    label <- str_replace_all(label, regex("rolling5y\\s+", ignore_case = TRUE), "")
    label <- str_replace_all(label, regex("vasicek", ignore_case = TRUE), "Vasicek")
    label <- str_replace_all(label, regex("cir", ignore_case = TRUE), "CIR")
    return(str_squish(label))
  }
  
  if (str_detect(source_low, "ukf") || str_detect(object_low, "ukf")) {
    prefix <- "UKF"
  } else if (str_detect(label_low, "initial")) {
    prefix <- ""
  } else {
    prefix <- "EKF"
  }
  
  label <- label %>%
    str_replace_all(regex("initial", ignore_case = TRUE), "Initial Model") %>%
    str_replace_all(regex("model", ignore_case = TRUE), "Model") %>%
    str_replace_all(regex("opt\\.?\\s*1", ignore_case = TRUE), "Opt. 1") %>%
    str_replace_all(regex("opt\\.?\\s*2", ignore_case = TRUE), "Opt. 2") %>%
    str_replace_all(regex("opt\\.?\\s*3", ignore_case = TRUE), "Opt. 3") %>%
    str_replace_all(regex("opt\\.?\\s*4", ignore_case = TRUE), "Opt. 4") %>%
    str_replace_all(regex("opt\\.?\\s*5", ignore_case = TRUE), "Opt. 5") %>%
    str_squish()
  
  if (label == "Initial Model" || str_detect(label, regex("Initial Model", ignore_case = TRUE))) {
    return("Initial Model")
  }
  
  if (!str_detect(label, regex("^(EKF|UKF)\\b", ignore_case = TRUE))) {
    label <- paste(prefix, label)
  }
  
  str_squish(label)
}

is_excluded_model <- function(object_name, model_label) {
  if (object_name %in% exclude_object_names) return(TRUE)
  any(str_detect(model_label, regex(paste(exclude_model_regex, collapse = "|"), ignore_case = TRUE)))
}

result_to_long <- function(res, object_name, source_name) {
  forecast <- first_existing(res, c("forecast", "forecast_mat", "pred", "pred_swap", "predicted"))
  actual   <- first_existing(res, c("actual", "actual_mat", "actual_swap", "observed"))
  dates    <- first_existing(res, c("dates", "Date", "forecast_dates"))
  errors   <- first_existing(res, c("errors", "error", "error_mat", "residuals", "residual"))
  rmse     <- first_existing(res, c("rmse", "RMSE"))
  
  if (is.null(forecast) || is.null(actual) || is.null(dates)) {
    stop("Mangler forecast/actual/dates i objektet: ", object_name)
  }
  
  forecast <- as.data.frame(forecast)
  actual   <- as.data.frame(actual)
  
  if (ncol(forecast) != ncol(actual)) {
    stop("forecast og actual har forskelligt antal kolonner i: ", object_name)
  }
  
  if (is.null(colnames(forecast)) || any(colnames(forecast) == "")) {
    colnames(forecast) <- maturity_order[seq_len(ncol(forecast))]
  }
  if (is.null(colnames(actual)) || any(colnames(actual) == "")) {
    colnames(actual) <- colnames(forecast)
  }
  
  model_label_raw <- first_existing(res, c("model_label", "model", "label", "method"))
  model_label <- clean_model_label(object_name, model_label_raw, source_name)
  
  if (is_excluded_model(object_name, model_label)) return(NULL)
  
  # Residualer som i dine forecast-scripts: actual - forecast.
  if (is.null(errors)) {
    errors <- as.matrix(actual) - as.matrix(forecast)
  }
  errors <- as.data.frame(errors)
  colnames(errors) <- colnames(forecast)
  
  out <- tibble(Date = as.Date(dates)) %>%
    bind_cols(
      forecast %>% setNames(paste0("forecast__", colnames(forecast))),
      actual   %>% setNames(paste0("actual__", colnames(actual))),
      errors   %>% setNames(paste0("residual__", colnames(errors)))
    ) %>%
    pivot_longer(
      cols = -Date,
      names_to = c("type", "maturity"),
      names_sep = "__",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = type, values_from = value) %>%
    mutate(
      source = source_name,
      object_name = object_name,
      model = model_label,
      residual = actual - forecast
    )
  
  if (!is.null(rmse)) {
    rmse_vec <- as.numeric(rmse)
    if (length(rmse_vec) == ncol(forecast)) {
      rmse_df <- tibble(
        model = model_label,
        source = source_name,
        object_name = object_name,
        maturity = colnames(forecast),
        RMSE = rmse_vec
      )
    } else {
      rmse_df <- out %>%
        group_by(model, source, object_name, maturity) %>%
        summarise(RMSE = sqrt(mean(residual^2, na.rm = TRUE)), .groups = "drop")
    }
  } else {
    rmse_df <- out %>%
      group_by(model, source, object_name, maturity) %>%
      summarise(RMSE = sqrt(mean(residual^2, na.rm = TRUE)), .groups = "drop")
  }
  
  list(long = out, rmse = rmse_df)
}

extract_results_from_file <- function(file, source_name) {
  if (!file.exists(file)) stop("Filen findes ikke: ", file)
  env <- load_env(file)
  
  if (exists("all_results", envir = env, inherits = FALSE)) {
    results <- get("all_results", envir = env)
  } else {
    candidates <- ls(env)
    results <- list()
    for (nm in candidates) {
      obj <- get(nm, envir = env)
      if (is.list(obj) && all(c("forecast", "actual", "dates") %in% names(obj))) {
        results[[nm]] <- obj
      }
    }
    if (length(results) == 0) {
      stop("Kunne ikke finde all_results eller forecast-resultatobjekter i: ", file)
    }
  }
  
  parts <- lapply(names(results), function(nm) result_to_long(results[[nm]], nm, source_name))
  parts <- Filter(Negate(is.null), parts)
  
  if (length(parts) == 0) {
    warning("Ingen modeller tilbage efter eksklusion i: ", file)
    return(NULL)
  }
  
  list(
    long = bind_rows(lapply(parts, `[[`, "long")),
    rmse = bind_rows(lapply(parts, `[[`, "rmse"))
  )
}

# -----------------------------
# SAMMENLÆG 
# -----------------------------
all_parts <- lapply(names(input_files), function(source_name) {
  extract_results_from_file(input_files[[source_name]], source_name)
})
all_parts <- Filter(Negate(is.null), all_parts)

forecast_long <- bind_rows(lapply(all_parts, `[[`, "long")) %>%
  mutate(
    maturity = factor(maturity, levels = maturity_order),
    maturity_num = as.numeric(gsub("Y", "", as.character(maturity)))
  )

rmse_table <- bind_rows(lapply(all_parts, `[[`, "rmse")) %>%
  mutate(
    maturity = factor(maturity, levels = maturity_order),
    maturity_num = as.numeric(gsub("Y", "", as.character(maturity)))
  )

model_name_map <- c(
  "EKF Opt. 1" = "Opt. 1 EKF",
  "EKF Opt. 2" = "Opt. 2 EKF",
  "EKF Opt. 3" = "Opt. 3 EKF",
  "EKF Opt. 5" = "Opt. 5 EKF",
  "UKF Opt. 1" = "Opt. 4 UKF",
  "UKF Opt. 2" = "Opt. 5 UKF",
  "Initial Model" = "Initial Model",
  "Vasicek" = "Vasicek",
  "CIR" = "CIR"
)

rename_model <- function(x) {
  x <- as.character(x)
  dplyr::recode(x, !!!model_name_map, .default = x)
}

forecast_long <- forecast_long %>%
  mutate(model = rename_model(model))

rmse_table <- rmse_table %>%
  mutate(model = rename_model(model))

model_order <- c(
  "Initial Model",
  "Opt. 1 EKF",
  "Opt. 5 EKF",
  "Opt. 4 UKF",
  "Opt. 5 UKF",
  "Vasicek",
  "CIR"
)

used_models <- model_order[model_order %in% unique(as.character(forecast_long$model))]
used_cols <- model_cols[c("Observed", used_models)]
used_model_cols <- model_cols[used_models]
used_model_cols <- used_model_cols[!is.na(used_model_cols)]

# Sæt faktor-levels efter used_models er fastlagt.
forecast_long <- forecast_long %>%
  mutate(model = factor(as.character(model), levels = used_models))

rmse_table <- rmse_table %>%
  mutate(model = factor(as.character(model), levels = used_models)) %>%
  arrange(model, maturity_num)

rmse_by_model <- forecast_long %>%
  group_by(model) %>%
  summarise(
    MSE = mean(residual^2, na.rm = TRUE),
    RMSE = sqrt(MSE),
    .groups = "drop"
  ) %>%
  arrange(RMSE)


if (isTRUE(save_tables)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  write_csv(forecast_long, file.path(output_dir, "combined_forecast_predicted_actual_residuals_long.csv"))
  write_csv(rmse_table, file.path(output_dir, "combined_rmse_by_model_and_maturity.csv"))
  write_csv(rmse_by_model, file.path(output_dir, "combined_rmse_by_model.csv"))
}


plot_predicted_base <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  dates_all <- forecast_long$Date
  year_ticks <- seq(
    from = as.Date(format(min(dates_all, na.rm = TRUE), "%Y-01-01")),
    to   = as.Date(format(max(dates_all, na.rm = TRUE), "%Y-01-01")),
    by   = "2 years"
  )
  
  par(
    mfrow = c(3, 3),
    mar = c(1.8, 1.8, 1.2, 0.2),
    oma = c(7, 3.5, 1, 0),
    mgp = c(1.4, 0.4, 0),
    tcl = -0.25
  )
  
  for (mat in maturity_order) {
    df_m <- forecast_long %>% filter(as.character(maturity) == mat)
    actual_m <- df_m %>% distinct(Date, actual) %>% arrange(Date)
    
    y_min <- min(c(actual_m$actual, df_m$forecast), na.rm = TRUE)
    y_max <- max(c(actual_m$actual, df_m$forecast), na.rm = TRUE)
    
    plot(
      actual_m$Date,
      actual_m$actual,
      type = "l",
      main = mat,
      xlab = "",
      ylab = "",
      lwd = 1.8,
      col = "black",
      ylim = c(y_min, y_max),
      xaxt = "n",
      cex.main = 1.45,
      cex.axis = 1.4
    )
    
    axis(1, at = year_ticks, labels = format(year_ticks, "%Y"), cex.axis = 1.4)
    abline(v = year_ticks, col = "grey85", lwd = 0.8)
    abline(h = axTicks(2), col = "grey85", lwd = 0.8)
    lines(actual_m$Date, actual_m$actual, lwd = 1.8, col = "black")
    
    for (mod in used_models) {
      df_line <- df_m %>% filter(as.character(model) == mod) %>% arrange(Date)
      lines(
        df_line$Date,
        df_line$forecast,
        lwd = 1.5,
        lty = 1,
        col = used_model_cols[mod]
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
    horiz = TRUE,
    bty = "n",
    cex = 1.05,
    xpd = NA,
    legend = c("Observed", used_models),
    col = used_cols,
    lty = 1,
    lwd = c(1.8, rep(1.5, length(used_models)))
  )
}





plot_residuals_base <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  dates_all <- forecast_long$Date
  year_ticks <- seq(
    from = as.Date(format(min(dates_all, na.rm = TRUE), "%Y-01-01")),
    to   = as.Date(format(max(dates_all, na.rm = TRUE), "%Y-01-01")),
    by   = "2 years"
  )
  
  par(
    mfrow = c(3, 3),
    mar = c(1.8, 1.8, 1.2, 0.2),
    oma = c(7, 3.5, 1, 0),
    mgp = c(1.4, 0.4, 0),
    tcl = -0.25
  )
  
  for (mat in maturity_order) {
    df_m <- forecast_long %>% filter(as.character(maturity) == mat)
    
    y_min <- min(df_m$residual, na.rm = TRUE)
    y_max <- max(df_m$residual, na.rm = TRUE)
    
    first_mod <- used_models[1]
    first_df <- df_m %>% filter(as.character(model) == first_mod) %>% arrange(Date)
    
    plot(
      first_df$Date,
      first_df$residual,
      type = "l",
      main = mat,
      xlab = "",
      ylab = "",
      lwd = 1.5,
      col = used_model_cols[first_mod],
      ylim = c(y_min, y_max),
      xaxt = "n",
      cex.main = 1.45,
      cex.axis = 1.4
    )
    
    axis(1, at = year_ticks, labels = format(year_ticks, "%Y"), cex.axis = 1.4)
    abline(v = year_ticks, col = "grey85", lwd = 0.8)
    abline(h = axTicks(2), col = "grey85", lwd = 0.8)
    abline(h = 0, col = "black", lty = 2, lwd = 1)
    
    for (mod in used_models) {
      df_line <- df_m %>% filter(as.character(model) == mod) %>% arrange(Date)
      lines(
        df_line$Date,
        df_line$residual,
        lwd = 1.5,
        lty = 1,
        col = used_model_cols[mod]
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
    horiz = TRUE,
    bty = "n",
    cex = 1.05,
    seg.len = 3,
    xpd = NA,
    legend = used_models,
    col = used_model_cols,
    lty = 1,
    lwd = 2.2
  )
}

p_rmse <- ggplot(
  rmse_table,
  aes(
    x = maturity_num,
    y = RMSE,
    color = model,
    group = model
  )
) +
  geom_line(linewidth = 1, linetype = "dashed", na.rm = TRUE) +
  geom_point(size = 2.5, na.rm = TRUE) +
  scale_color_manual(values = used_model_cols) +
  scale_x_continuous(breaks = maturity_breaks) +
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


plot_fit_terms_gg <- function(selected_terms = c(1, 3, 5, 10, 20, 30), plot_title = "") {
  actual_df <- forecast_long %>%
    distinct(Date, maturity_num, actual) %>%
    transmute(Date, Term_num = maturity_num, Term = paste0(maturity_num, "Y"), Model = "Observed", value = actual)
  
  forecast_df <- forecast_long %>%
    transmute(Date, Term_num = maturity_num, Term = paste0(maturity_num, "Y"), Model = as.character(model), value = forecast)
  
  plot_df <- bind_rows(actual_df, forecast_df) %>%
    filter(Term_num %in% selected_terms) %>%
    mutate(
      Term = factor(Term, levels = paste0(selected_terms, "Y")),
      Model = factor(Model, levels = c("Observed", used_models))
    )
  
  ggplot(plot_df, aes(x = Date, y = value, color = Model, linetype = Model)) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    facet_wrap(~ Term, ncol = 2, scales = "free_y") +
    scale_color_manual(values = used_cols) +
    scale_linetype_manual(values = rep("solid", length(used_cols))) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(title = plot_title, x = "Date", y = "Swap rate", color = NULL, linetype = NULL) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, keywidth = unit(1.3, "cm")), linetype = "none") +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = element_text(size = 9.5),
      legend.key.width = unit(1.3, "cm"),
      legend.spacing.x = unit(0.25, "cm"),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}

plot_residual_terms_gg <- function(selected_terms = c(1, 3, 5, 10, 20, 30), plot_title = "") {
  plot_df <- forecast_long %>%
    transmute(Date, Term_num = maturity_num, Term = paste0(maturity_num, "Y"), Model = as.character(model), value = residual) %>%
    filter(Term_num %in% selected_terms) %>%
    mutate(
      Term = factor(Term, levels = paste0(selected_terms, "Y")),
      Model = factor(Model, levels = used_models)
    )
  
  ggplot(plot_df, aes(x = Date, y = value, color = Model, linetype = Model)) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.5) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    facet_wrap(~ Term, ncol = 2, scales = "free_y") +
    scale_color_manual(values = used_model_cols) +
    scale_linetype_manual(values = rep("solid", length(used_model_cols))) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(title = plot_title, x = "Date", y = "Residual", color = NULL, linetype = NULL) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, keywidth = unit(1.3, "cm")), linetype = "none") +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = element_text(size = 10),
      legend.key.width = unit(1.3, "cm"),
      legend.spacing.x = unit(0.25, "cm"),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}


plot_predicted_base()
plot_residuals_base()
print(p_rmse)

cat("\nFærdig. Plots er vist i plot-vinduet. Ingen billeder er gemt automatisk.\n")
cat("\nModeller inkluderet:\n")
print(used_models)
cat("\nSamlet RMSE pr. model:\n")
print(rmse_by_model)
cat("\nRMSE pr. model og maturity:\n")
print(rmse_table, n = Inf)