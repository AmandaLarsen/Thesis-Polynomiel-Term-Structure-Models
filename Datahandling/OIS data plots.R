


library(readxl)


yield <- OIS_Zero_Yield
bond  <- OIS_Zero_Bond_

names(yield)[1] <- "Date"
names(bond)[1]  <- "Date"

yield <- yield %>%
  mutate(Date = as.Date(Date)) %>%
  filter(!is.na(Date)) %>%
  arrange(Date)


df <- df %>%
  mutate(across(matches("^Z_T|^Z_EUROON2W"), ~ . * 100))


bond <- bond %>%
  mutate(Date = as.Date(Date)) %>%
  filter(!is.na(Date)) %>%
  arrange(Date)






library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# -------------------------
# Zero yields
# -------------------------

plot_window_yield <- yield %>%
  filter(Date >= as.Date("2019-06-01"),
         Date <= as.Date("2020-03-01")) %>%
  select(Date, Z_T_1.0, Z_T_5.0, Z_T_10.0, Z_T_15.0, Z_T_20.0, Z_T_25.0, Z_T_30.0) %>%
  pivot_longer(
    cols = -Date,
    names_to = "Tenor",
    values_to = "ZeroRate"
  ) %>%
  mutate(
    Tenor = factor(
      Tenor,
      levels = c("Z_T_1.0", "Z_T_5.0","Z_T_10.0","Z_T_15.0","Z_T_20.0", "Z_T_25.0", "Z_T_30.0"),
      labels = c("1Y", "5Y", "10Y", "15Y", "20Y", "25Y", "30Y")
    )
  )

ggplot(plot_window_yield, aes(x = Date, y = ZeroRate, color = Tenor)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = as.Date("2019-10-02"),
             linetype = "dashed") +
  scale_color_viridis_d(option = "I", direction = -1) +
  theme_light() +
  labs(
    x = NULL,
    y = "Zero rate (%)",
    color = "Term"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )




#### Andet plot
ggplot(plot_window_yield, aes(x = Date, y = ZeroRate, color = Tenor)) +
  geom_line(linewidth = 1.2) +
  geom_vline(
    xintercept = as.Date("2019-10-02"),
    linetype = "dashed",
    linewidth = 0.9
  ) +
  scale_color_viridis_d(option = "I", direction = -1) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b %Y"
  ) +
  labs(
    x = NULL,
    y = "Zero rate (%)",
    color = "Term"
  ) +
  theme_light(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    
    axis.title = element_text(size = 20),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 14
    ),
    axis.text.y = element_text(size = 16),
    
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    
    plot.title = element_text(
      face = "bold",
      size = 22
    )
  )



# -------------------------
# discount factors
# -------------------------

plot_window_bond <- bond %>%
  filter(Date >= as.Date("2019-06-01"),
         Date <= as.Date("2020-03-01")) %>%
  select(Date, P_T_1.0, P_T_5.0, P_T_10.0, P_T_15.0, P_T_20.0, P_T_25.0, P_T_30.0) %>%
  pivot_longer(
    cols = -Date,
    names_to = "Tenor",
    values_to = "DiscountFactor"
  ) %>%
  mutate(
    Tenor = factor(
      Tenor,
      levels = c("P_T_1.0", "P_T_5.0", "P_T_10.0", "P_T_15.0", "P_T_20.0", "P_T_25.0", "P_T_30.0"),
      labels = c("1Y", "5Y", "10Y", "15Y", "20Y", "25Y", "30Y")
    )
  )

ggplot(plot_window_bond, aes(x = Date, y = DiscountFactor, color = Tenor)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = as.Date("2019-10-02"),
             linetype = "dashed") +
  scale_color_viridis_d(option = "I", direction = -1) +
  theme_light() +
  labs(
    x = NULL,
    y = "Discount factor",
    color = "Term"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )


##### Andet plot


ggplot(
  plot_window_bond,
  aes(x = Date, y = DiscountFactor, color = Tenor)
) +
  
  geom_line(linewidth = 1.2) +
  
  geom_vline(
    xintercept = as.Date("2019-10-02"),
    linetype = "dashed",
    linewidth = 0.9
  ) +
  
  scale_color_viridis_d(
    option = "I",
    direction = -1
  ) +
  
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b %Y"
  ) +
  
  labs(
    x = NULL,
    y = "Discount factor",
    color = "Term"
  ) +
  
  theme_light(base_size = 18) +
  
  theme(
    panel.grid.minor = element_blank(),
    
    axis.title = element_text(size = 20),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 14
    ),
    
    axis.text.y = element_text(size = 16),
    
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    
    plot.title = element_text(
      face = "bold",
      size = 22
    )
  )


# -------------------------
# alle år
# -------------------------



plot_all_bond <- bond %>%
  select(Date, P_T_1.0, P_T_5.0, P_T_10.0, P_T_15.0, P_T_20.0, P_T_25.0, P_T_30.0) %>%
  pivot_longer(
    cols = -Date,
    names_to = "Tenor",
    values_to = "DiscountFactor"
  ) %>%
  mutate(
    Tenor = factor(
      Tenor,
      levels = c("P_T_1.0", "P_T_5.0", "P_T_10.0", "P_T_15.0", "P_T_20.0", "P_T_25.0", "P_T_30.0"),
      labels = c("1Y", "5Y", "10Y", "15Y", "20Y", "25Y", "30Y")
    )
  )

ggplot(plot_all_bond, aes(x = Date, y = DiscountFactor, color = Tenor)) +
  geom_line(linewidth = 0.8) +
  scale_color_viridis_d(option = "I", direction = -1) +
  theme_light() +
  labs(
    x = NULL,
    y = "Discount factor",
    color = "Term"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
#### andet plot

ggplot(plot_all_bond,
       aes(x = Date,
           y = DiscountFactor,
           color = Tenor)) +
  
  geom_line(linewidth = 1.2) +
  
  scale_color_viridis_d(
    option = "I",
    direction = -1
  ) +
  
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  
  labs(
    x = "Date",
    y = "Discount factor",
    color = "Term"
  ) +
  
  theme_light(base_size = 18) +
  
  theme(
    panel.grid.minor = element_blank(),
    
    axis.title = element_text(size = 20),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 14
    ),
    
    axis.text.y = element_text(size = 16),
    
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    
    plot.title = element_text(
      face = "bold",
      size = 22
    ),
    
    legend.position = "right"
  )



dates_to_plot <- as.Date(c(
  "2008-10-17",
  "2011-03-10",
  "2016-09-16",
  "2018-12-12"
))

idx <- match(dates_to_plot, df$Date)

if (any(is.na(idx))) {
  stop(
    "Følgende datoer findes ikke i datasættet: ",
    paste(dates_to_plot[is.na(idx)], collapse = ", ")
  )
}

plot_df <- bind_rows(
  lapply(idx, function(i) {
    data.frame(
      Date = df$Date[i],
      maturity = full_maturity_grid,
      price = as.numeric(
        get_full_price_curve(
          curve_list[[i]],
          very_short_name,
          short_tenors,
          semi_grid
        )
      )
    )
  })
)

plot_df$Date_lab <- paste("", plot_df$Date)

ggplot(plot_df, aes(x = maturity, y = price)) +
  geom_line(color = "black", linewidth = 0.5) +
  geom_point(color = "black", size = 1.4) +
  facet_wrap(~ Date_lab, ncol = 2, scales = "free_y") +
  labs(
    x = "Maturity",
    y = "Discount Curve"
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

plot_df_zero_both$Panel <- paste("", plot_df_zero_both$Date)

p_zero_both <- ggplot(plot_df_zero_both, aes(x = maturity, y = value, linetype = Curve, shape = Curve)) +
  geom_line(color = "black", linewidth = 0.6) +
  geom_point(color = "black", size = 1.6) +
  facet_wrap(~ Panel, ncol = 1, scales = "free_y") +
  labs(
    title = "",
    x = "Maturity",
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






library(dplyr)
library(ggplot2)

dates_to_plot <- as.Date(c(
  "2008-10-17",
  "2011-03-10",
  "2016-09-16",
  "2018-12-12"
))

if (!all(dates_to_plot %in% bond_prices_df$Date)) {
  stop("En eller flere datoer findes ikke i datasættet.")
}

plot_df_multi <- bind_rows(
  lapply(dates_to_plot, function(d) {
    data.frame(
      Date = d,
      maturity = 1:30,
      price = as.numeric(
        bond_prices_df[bond_prices_df$Date == d, paste0("P_Y", 1:30)]
      )
    )
  })
)

plot_df_multi$Panel <- paste("", plot_df_multi$Date)

ggplot(plot_df_multi, aes(x = maturity, y = price)) +
  geom_line(color = "black", linewidth = 0.6) +
  geom_point(color = "black", size = 1.5) +
  facet_wrap(~ Panel, ncol = 2, scales = "free_y") +
  labs(
    x = "Maturity",
    y = "Discound curve"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
















