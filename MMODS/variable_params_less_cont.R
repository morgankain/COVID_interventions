fixed_params        <- fixed_params[-c(1, 5, 11, 12)]

set.seed(1022)

## parameters that will vary
variable_params <- sobolDesign(
  lower = c(
    E0          = 3
  , sim_start   = as.numeric(as.Date("2020-03-03") - as.Date("2020-01-01"))
  , Ca          = 0.4
  , alpha       = 0.3
  , delta       = 0.1
  , mu          = 0.925
  )
, upper = c(
    E0          = 6
  , sim_start   = as.numeric(as.Date("2020-03-15") - as.Date("2020-01-01"))
  , Ca          = 0.8
  , alpha       = 0.5
  , delta       = 0.3
  , mu          = 0.975
)
, nseq  = nparams
) %>% mutate(
  sim_start   = round(sim_start)
, E0 = round(E0))

if (fit.E0) {
variable_params <- variable_params %>% dplyr::select(-E0)
}

variable_params <- variable_params %>% mutate(
    sim_start  = as.Date(sim_start, origin = as.Date("2020-01-01"))
 ) %>% mutate(paramset = seq(1, nparams)) 

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
