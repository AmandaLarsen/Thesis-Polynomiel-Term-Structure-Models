load("~/Desktop/ukf_workspace_tilfoej_a_uden_udvidet_int.RData")








#######################
#### testttt        ###
#######################

payer_swap_value <- function(state, pars, Tn, K) {
  
  x <- state[1]
  y <- state[2]
  
  kx <- pars["kappa_x"]
  thx <- pars["theta_x"]
  sx <- pars["sigma_x"]
  
  ky <- pars["kappa_y"]
  thy <- pars["theta_y"]
  sy <- pars["sigma_y"]
  
  alpha <- pars["alpha"]
  eta <- pars["eta"]
  a <- pars["a"]
  
  denom <- a + x + x^2
  
  P_fun <- function(Tau) {
    EX  <- vasicek_mean(x, Tau, kx, thx)
    EX2 <- vasicek_second_moment(x, Tau, kx, thx, sx)
    exp(-alpha * Tau) * (a + EX + EX2) / denom
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
      deltas <- c(rep(0.25, length(starts_3m)),
                  rep(0.5, length(starts_6m)))
    }
    sum(mapply(A_period, starts, deltas))
  }
  
  pay_dates <- seq(0.5, Tn, by = 0.5)
  annuity <- fixed_delta * sum(sapply(pay_dates, P_fun))
  
  value <- 1 - P_fun(Tn) - K * annuity + A_sum(Tn)
  return(value)
}



state_t <- tail(ukf_out_real$x_filt, 1)
pars <- ukf_out_real$pars



S_model <- model_swap_rates(state_t, pars)

fair_values <- sapply(seq_along(swap_maturities), function(i) {
  payer_swap_value(
    state = state_t,
    pars = pars,
    Tn = swap_maturities[i],
    K = S_model[i]
  )
})

names(fair_values) <- swap_cols
fair_values




state_t <- tail(ukf_out_real$x_pred, 1)
pars <- ukf_out_real$pars

S_pred <- as.numeric(tail(swap_pred_real, 1))

fair_values_pred <- sapply(seq_along(swap_maturities), function(i) {
  payer_swap_value(
    state = state_t,
    pars  = pars,
    Tn    = swap_maturities[i],
    K     = S_pred[i]
  )
})

names(fair_values_pred) <- swap_cols
fair_values_pred


##################
## noget andet  ##
##################




historical_payer_swap_mtm <- function(t0, Tn = 5) {
  
  j <- which(swap_maturities == Tn)
  if (length(j) != 1) stop("Tn must be one of swap_maturities")
  
  # Model-implied par swap rate at inception
  K <- model_swap_rates(
    state = ukf_out_real$x_filt[t0, ],
    pars  = ukf_out_real$pars
  )[j]
  
  mtm_values <- sapply(t0:nrow(ukf_out_real$x_filt), function(i) {
    payer_swap_value(
      state = ukf_out_real$x_filt[i, ],
      pars  = ukf_out_real$pars,
      Tn    = Tn,
      K     = K
    )
  })
  
  data.frame(
    Date = dates[t0:nrow(ukf_out_real$x_filt)],
    Maturity = paste0(Tn, "Y"),
    K = K,
    MTM = mtm_values
  )
}


start_years <- c(2016, 2018, 2020, 2022)

t0_list <- sapply(start_years, function(yy) {
  which.min(abs(dates - as.Date(paste0(yy, "-01-01"))))
})

mtm_multi_5y <- do.call(rbind, lapply(t0_list, function(t0) {
  
  out <- historical_payer_swap_mtm(
    t0 = t0,
    Tn = 5
  )
  
  out$StartDate <- dates[t0]
  out$Contract <- paste0(" ", format(dates[t0], "%B %Y"))
  
  out
}))

ggplot(mtm_multi_5y, aes(x = Date, y = MTM, color = Contract)) +
  geom_line(linewidth = 1) +
  geom_point(data = mtm_multi_5y |> dplyr::group_by(Contract) |> dplyr::slice(1),
             size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "",
    x = NULL,
    y = "Swap value",
    color = "Contract Start"
  ) +
  theme_light() +
  theme(legend.position = "bottom")







###

start_years <- c(2016, 2018, 2020, 2022)

t0_list <- sapply(start_years, function(yy) {
  which.min(abs(dates - as.Date(paste0(yy, "-01-01"))))
})



######

all_mtm <- do.call(rbind,
                   
                   lapply(swap_maturities, function(Tn) {
                     
                     mtm_Tn <- do.call(rbind,
                                       
                                       lapply(t0_list, function(t0) {
                                         
                                         out <- historical_payer_swap_mtm(
                                           t0 = t0,
                                           Tn = Tn
                                         )
                                         
                                         out$Contract <- format(
                                           dates[t0],
                                           "%B %Y"
                                         )
                                         
                                         out
                                       })
                                       
                     )
                     
                     mtm_Tn$Maturity <- paste0(Tn, "Y")
                     
                     mtm_Tn
                   })
                   
)


library(ggplot2)


all_mtm$Maturity <- factor(
  all_mtm$Maturity,
  levels = c(
    "1Y",
    "2Y",
    "3Y",
    "5Y",
    "7Y",
    "10Y",
    "15Y",
    "20Y",
    "30Y"
  )
)


ggplot(
  all_mtm,
  aes(
    x = Date,
    y = MTM,
    colour = Contract
  )
) +
  geom_line(linewidth = 1) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~ Maturity,
    ncol = 3,
    scales = "free_y"
  ) +
  labs(
    title = "",
    x = NULL,
    y = "Swap value",
    colour = "Contract start"
  ) +
  theme_light() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )



