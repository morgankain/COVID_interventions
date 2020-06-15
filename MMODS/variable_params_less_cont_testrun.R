fixed_params        <- fixed_params[-c(1, 5, 11, 12)]

test_simstart <- seq(
    as.numeric((as.Date(county.data[min(which(county.data$cases > 0)), ]$date) - 17) - as.Date("2020-01-01"))
  , as.numeric((as.Date(county.data[min(which(county.data$cases > 0)), ]$date) - 3) - as.Date("2020-01-01"))
  , by = 1
  )

## parameters that will vary
variable_params <- data.frame(
  sim_start = test_simstart
, Ca        = as.numeric(rep(params[params$Parameter == "Ca", ]$est, length(test_simstart)))
, alpha     = as.numeric(rep(params[params$Parameter == "alpha", ]$est, length(test_simstart)))
, delta     = as.numeric(rep(params[params$Parameter == "delta", ]$est, length(test_simstart)))
, mu        = as.numeric(rep(params[params$Parameter == "mu", ]$est, length(test_simstart)))
) %>% mutate(sim_start   = round(sim_start))

variable_params <- variable_params %>% mutate(
    sim_start  = as.Date(sim_start, origin = as.Date("2020-01-01"))
 ) %>% mutate(paramset = seq(1, length(test_simstart))) 

variable_params <- variable_params %>% 
  mutate(
    log_lik        = 0
  , log_lik_se     = 0
  , beta0est       = 0
  , E_init         = 0
  , detect_k       = 0
  , detect_mid     = 0
  , theta          = 0
  , theta2         = 0
  , beta_min       = 0
  ) 

nparams           <- nrow(variable_params)
