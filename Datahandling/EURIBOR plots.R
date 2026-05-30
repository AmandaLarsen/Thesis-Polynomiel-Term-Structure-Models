

library(tidyverse)

df_long <- EURIBOR_swap_data%>%
  mutate(Date = as.Date(Date)) %>%
  pivot_longer(
    cols = matches("^Y\\d+"),
    names_to = "Term",
    values_to = "Yield"
  ) %>%
  mutate(
    Term = factor(Term, levels = c("Y1","Y2","Y3","Y5","Y7","Y10","Y15","Y20","Y30"))
  )



df_long <- df_long %>%
  mutate(
    Term = factor(Term, levels = c("Y1","Y2","Y3","Y5","Y7","Y10","Y15","Y20","Y30"))
  )


ggplot(df_long, aes(x = Date, y = Yield, color = Term)) +
  geom_line(linewidth = 0.9) +
  theme_light() +
  scale_color_viridis_d(option = "I", direction = -1) +
  labs(
    x = NULL,
    y = "Swap Rate",
    color = "Term"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 15),
    panel.grid.minor = element_blank()
  )

## Ændrede plots
ggplot(df_long, aes(x = Date, y = Yield, color = Term)) +
  
  geom_line(linewidth = 1.2) +
  
  scale_color_viridis_d(option = "I", direction = -1) +
  
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  
  labs(
    x = "Date",
    y = "Swap Rate",
    color = "Term"
  ) +
  
  theme_light(base_size = 18) +
  
  theme(
    legend.position = "right",
    
    axis.title = element_text(size = 20),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 14
    ),
    
    axis.text.y = element_text(size = 16),
    
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    
    panel.grid.minor = element_blank()
  )



