## Required packages to run this code
needed_packages <- c(
    "pomp"
  , "plyr"
  , "dplyr"
  , "ggplot2"
  , "magrittr"
  , "scales"
  , "lubridate"
  , "tidyr"
  , "foreach"
  , "doParallel"
  , "data.table"
  , "doRNG")

## load packages. Install all packages that return "FALSE"
lapply(needed_packages, require, character.only = TRUE)

if (focal.county == "Fulton") {

deaths <- read.csv("us-counties.txt") %>%
  mutate(date = as.Date(date)) %>%
  dplyr::filter((county == focal.county | county == "DeKalb") & state == focal.state) %>%
  group_by(date) %>% summarize(cases = sum(cases), deaths = sum(deaths)) %>%
  mutate(county = focal.county, state = focal.state)

} else {
  
deaths  <- read.csv("us-counties.txt") %>% 
  mutate(date = as.Date(date)) %>% 
  dplyr::filter(county == focal.county & state == focal.state) %>% 
  dplyr::filter(date < max(date) - fit.minus) %>% 
  filter(state == focal.state)
  
}
  
## Load the previously saved fits
 ## If COVID_fit_cont.R was just run, use parameers already stored in the global env 
if (use.rds) {
prev.fit         <- readRDS(rds.name)
variable_params  <- prev.fit$mifs_local_v2
covid_mobility   <- prev.fit[["covid_mobility"]] # get pomp object
# variable_params <- prev.fit[["variable_params"]]
# fixed_params    <- prev.fit[["fixed_params"]]
}
  
if (loglik.max) {
  
variable_params <- variable_params %>% 
  filter(log_lik != 0) %>% 
  filter(log_lik == max(log_lik))

print(variable_params$paramset)
print(variable_params$log_lik)

} else {
  
  if (is.na(loglik.num)) {
    
variable_params <- variable_params %>% 
filter(log_lik != 0) %>% 
filter(log_lik > (max(log_lik) - loglik.thresh)) 

print(variable_params$paramset)
print(variable_params$log_lik) 

  } else {
    
variable_params <- variable_params %>% 
filter(log_lik != 0) %>% 
arrange(desc(log_lik)) %>%
slice(1:loglik.num)

print(variable_params$paramset)
print(variable_params$log_lik)

}
}

pomp_data        <- data.frame(covid_mobility) %>% 
  full_join(as.data.frame(cbind(day = covid_mobility@covar@times
    , t(covid_mobility@covar@table)))) %>%
    arrange(day) %>% 
    distinct_at(vars(day), .keep_all = T) # sometimes the pomp data and covariate data disagree? not sure why..., for now default to takign the pomp data
  dt             <- covid_mobility@rprocess@delta.t
  date_origin    <- prev.fit[["date_origin"]][1]

if (int.catch_eff_post == 0) {

SIP_post <- check_R0(
    beta0est     = variable_params$beta0
  , beta_min     = variable_params$beta_min
  , fixed_params = variable_params
  , sd_strength  = 1
  , prop_S       = 1
  , desired_R    = desired.R
  )

} else {
  
SIP_post <- sip_trunc_combns(
                              beta_catch      = int.beta_catch_post
                            , beta_catch_type = "pct"
                            , catch_eff       = int.catch_eff_post
                            , k               = variable_params$beta0_k
                            , beta0           = variable_params$beta0
                            , beta_min        = variable_params$beta_min
                            , d               = variable_params$d
                            , dt              = dt
                            , desired_R       = desired.R)
  
}

if (focal.county == "Fulton") {
  mobility <- readRDS(mobility.file) %>% 
    dplyr::filter((county_name == focal.county | county_name == "DeKalb") & (state_abbr == focal.state_abbr)) %>%
    dplyr::group_by(datestr) %>%
    dplyr::summarize(sip_prop = mean(sip_prop), .groups = "drop") 
} else {
  mobility <- readRDS(mobility.file) %>% 
    dplyr::filter(county_name == focal.county & state_abbr == focal.state_abbr)
}

# clean up mobility data by taking 3 day moving median and adding days since date_origin
mobility %<>% 
  dplyr::select(datestr, sip_prop) %>% 
  filter(sip_prop != 0) %>% 
  {full_join(., # join with full list of dates so NAs will occur for missing days
             data.frame(datestr = seq(pull(., datestr) %>% min, 
                                      pull(., datestr) %>% max, 
                                      by = "day")))} %>% 
  mutate(day = as.numeric(datestr - date_origin)) %>% 
  arrange(day) %>% 
  # three day moving median
  mutate(sip_prop = mutate(., 
                           sip_prop_lag = lag(sip_prop),
                           sip_prop_lead = lead(sip_prop)) %>% 
           select(contains("sip_prop")) %>%
           purrr::pmap_dbl(~median(c(...), na.rm = T))) %>% 
  dplyr::filter(!is.na(sip_prop)) 

# backfill mobility to date_origin
if (min(mobility$day) > 1) {
  mobility <- rbind(
    data.frame(
      datestr  = seq.Date(date_origin + 1, 
                          min(mobility$datestr) - 1,
                          "day")
      , sip_prop = mean(mobility$sip_prop[1:10])
    ) %>% mutate(day = as.numeric(datestr - date_origin)) 
    , mobility
  )  
}

## Adjust the cases/deaths data to anchor dates from date_origin
county.data <- deaths %>% 
  mutate(day = as.numeric(date - date_origin)) %>% 
  dplyr::select(day, date, deaths, cases) %>% 
  arrange(date) %>%
  # convert to daily deaths and cases
  mutate(deaths = deaths - lag(deaths),
         cases = cases - lag(cases)) %>%
  # convert negatve cases and deaths to zeros
  dplyr::rowwise() %>%
  mutate(cases = max(0, cases),
         deaths = max(0, deaths)) %>%
  ungroup() %>%
  dplyr::filter(!is.na(deaths), !is.na(cases)) %>%
  # fill with NAs from the latest considered sim start date to the data, 
  # NAs are to prevent any oddities with the accumulator variables
  # if the last sim_start date considered occurs when or after the data start, problems will still occur
  {if(max(variable_params[, "sim_start"]) < min(pull(., "date"))) {rbind(
    data.frame(
      date   = seq.Date(variable_params[, "sim_start"] %>% max, 
                        min(pull(., "date")) - 1,
                        "day")
      , deaths = NA
      , cases  = NA
    ) %>% 
      mutate(day = as.numeric(date - date_origin))
    , .
  )} else{ # if the data start before the last sim start date, trim the data to the last sim start and add one row of NAs to deal with accumulator variable
    filter(., date >=  variable_params[, "sim_start"] %>% max) %>% 
      {rbind(data.frame(date = min(pull(., "date")) - 1,
                       deaths = NA, 
                       cases = NA) %>% 
               mutate(day = as.numeric(date - date_origin)),
             .)}
  ## Remove dates after one week after movement data ends
  }} %>% dplyr::filter(date < (max(mobility$datestr) + 7))

Reff   <- data.frame(paramset = 0, date = 0, Reff = 0)  ; Reff   <- Reff[-1, ]
betat  <- data.frame(paramset = 0, date = 0, betat = 0) ; betat  <- betat[-1, ]
detect <- data.frame(paramset = 0, date = 0, detect = 0); detect <- detect[-1, ]

## container for the extra parameters needed
post_params <- c()
