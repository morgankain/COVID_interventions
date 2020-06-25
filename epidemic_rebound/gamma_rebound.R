### Simulate trajectories under a gamma distributed beta distribution and compare it to the mean

source("epidemic_rebound/gamma_rebound_setup.R", local = T)

####
## Simulate: loop over rows of variable_params
####

for (i in 1:nrow(variable_params)) {

print(rds.name)
print(variable_params$paramset[i])
print(variable_params$log_lik[i])

if (counter.factual) {
  
if (cf.type == "no_int") {
mob.covtab <- covariate_table(
   sip_prop         = rep(mean(head(mobility$sip_prop, 7)), sim_length)
 , order            = "constant"
 , times            = seq(1, sim_length, by = 1)
 , detect_t0        = min(county.data[!is.na(county.data$cases), ]$day) - 1
 , iso_mild_level   = iso_mild_level
 , iso_severe_level = iso_severe_level
 , intervention     = rep(0, sim_length))  
} else if (cf.type == "delay") {
mob.covtab <- covariate_table(
   sip_prop         = c(rep(mean(head(mobility$sip_prop, 7)), 7), mobility$sip_prop
     , rep(mean(tail(mobility$sip_prop, 3)), sim_length - (length(mobility$sip_prop) + 7))) 
 , order            = "constant"
 , times            = seq(1, sim_length, by = 1)
 , detect_t0        = min(county.data[!is.na(county.data$cases), ]$day) - 1
 , iso_mild_level   = iso_mild_level
 , iso_severe_level = iso_severe_level
 , intervention     = rep(0, sim_length))  
} else if (cf.type == "may1") { 
  
int.begin    <- as.numeric((as.Date("2020-05-01") - variable_params[i, ]$sim_start))
int.phase1   <- nrow(mobility) + (int.begin - nrow(mobility))
  
mob.covtab <- covariate_table(
   sip_prop         = c(
      mobility$sip_prop[1:int.begin]
    , rep(mean(head(mobility$sip_prop, 7)), sim_length - int.phase1)
      )
 , order            = "constant"
 , times            = seq(1, sim_length, by = 1)
 , detect_t0        = min(county.data[!is.na(county.data$cases), ]$day) - 1
 , iso_mild_level   = iso_mild_level
 , iso_severe_level = iso_severe_level
 , intervention     = rep(0, sim_length))  
}
  
} else {
  
####
## Set up intervention scenario
####

## Covariate table with all of the intervention scenarios
source("epidemic_rebound/gamma_rebound_covar.R", local = T)

}

covid_mobility <- pomp(
  data         = county.data %>% select(day, cases, deaths)
  , times      = "day"
  , t0         = 1
  , covar      = mob.covtab
  , rprocess   = euler(sir_step_mobility, delta.t = dt)
  , rmeasure   = {if(con_theta){rmeas_multi_logis_con}else{rmeas_multi_logis_ind}}
  , dmeasure   = {if(con_theta){dmeas_multi_logis_con}else{dmeas_multi_logis_ind}}
  , rinit      = sir_init
  , partrans   = {if(con_theta){par_trans_con}else{par_trans_ind}}
  , accumvars  = accum_names
  , paramnames = param_names
  , statenames = state_names
  , globals    = globs
)  

## Second bit of the intervention will affect these parameters:
variable_params[i, "beta0_k"]      <- int.beta0_k  ## heterogeneity value
variable_params[i, "beta_catch"]   <- {                ## beta0 values caught by intervention
  if(int.beta_catch_type == "pct") {
    max(qgamma(p = (1 - int.beta_catch), 
                                shape = int.beta0_k*dt/variable_params[i, "d"], 
                                scale = variable_params[i, "beta0"]/int.beta0_k),
                         0.001)
  } else{
    max(int.beta_catch, 0.001) # crude way to ensure int.beta_catch > 0, note that even this nonzero low beta_catch will be problematic but will eventually finish
  }
}
variable_params[i, "catch_eff"] <- int.catch_eff    ## beta0 values caught by intervention

post_params["sip_prop_post"]    <- SIP_post

post_params["beta_catch_post"]  <- {                ## beta0 values caught by intervention
  if(int.beta_catch_type == "pct") {
    max(qgamma(p = (1 - int.beta_catch_post), 
                                shape = int.beta0_k_post*dt/variable_params[i, "d"], 
                                scale = variable_params[i, "beta0"]/int.beta0_k_post),
                         0.001)
  } else{
   max(int.beta_catch, 0.001) # crude way to ensure int.beta_catch > 0, note that even this nonzero low beta_catch will be problematic but will eventually finish
  }
}
post_params["catch_eff_post"]   <- int.catch_eff_post
post_params["beta0_k_post"]     <- int.beta0_k_post  ## heterogeneity value

print(int.type)
print(post_params["catch_eff_post"])
print(post_params["beta_catch_post"])

checktime <- system.time({

sim_times <- seq(as.numeric(variable_params[i,]$sim_start - date_origin), max(mob.covtab@times), 1) 

## Simulate from fits
    SEIR.sim <- do.call(
      pomp::simulate
      , list(
        object   = covid_mobility
        , t0     = as.numeric(variable_params[i,]$sim_start - date_origin)
        , times  = sim_times
        , params = c(variable_params[i,] %>% unlist %>% 
                       inset("beta_catch", value = beta_catch_val) %>% 
                       inset("beta0_k", value = int.beta0_k) %>% 
                       inset("catch_eff", value = int.catch_eff), unlist(post_params))
        , nsim         = nsim
        , format       = "d"
        , include.data = F
        , seed         = 1002))
  
})[3]

print(checktime)
      
SEIR.sim <- SEIR.sim %>%
  mutate(
      date     = day + date_origin
    , paramset = i)

if (ci.epidemic) {
  epi_ids <- SEIR.sim %>% 
    group_by(.id) %>% 
    summarise(total_infect = max(D + R), .groups = "drop") %>% 
    filter(total_infect > ci.epidemic_cut*ceiling(variable_params[i, "E_init"])) %>% 
    filter(total_infect > ci.epidemic_cut) %>% #*ceiling(variable_params[i, "E_init"])) %>% 
    pull(.id)
  print(paste0("limiting to epidemics, including ", length(epi_ids), " simulations"))
  SEIR.sim %<>% filter(.id %in% epi_ids)
}

SEIR.sim %<>% {
  rbind(.,
        group_by(., day) %>%
          dplyr::select(-.id) %>%
          summarise_all(median, .groups = "drop") %>%
          mutate(.id = "median"))
}

SEIR.sim.s <- SEIR.sim %>% 
  dplyr::filter(.id != "median") %>%
  group_by(day) %>%
  summarise(S = mean(S), D = mean(D), .groups = "drop") %>%
  mutate(date = day + date_origin)

betat.t <- data.frame(
  paramset = i
  , date     = mob.covtab@times + date_origin
  , betat    = variable_params[i, "beta0"] * exp(log(variable_params[i, "beta_min"])*mob.covtab@table[which(dimnames(mob.covtab@table)[[1]] == "sip_prop"), ])
)

# not the cleanest, but join the sim summary and betat to get the dates aligned first, then calculate Reff
Reff.t <- full_join(SEIR.sim.s, betat.t) %>% arrange(day) %>% 
  filter(!is.na(S), !is.na(D), !is.na(betat)) %>% 
  {data.frame(
      paramset = i
      , date     = pull(., date)
      , Reff     = covid_R0(
        beta0est      = pull(., betat)
        , fixed_params  = unlist(variable_params[i, ])
        , sd_strength   = 1
        , prop_S        = pull(.,S) / (variable_params[i, "N"] - pull(.,D))
      )
    )
    }

detect.t <- data.frame(
  paramset = i
  , date     = mob.covtab@times + date_origin
  , detect   = c(rep(0, min(pomp_data[!is.na(pomp_data$cases), "day"]) - 1)
                 , (variable_params[i, "detect_max"] /
                      (1 + exp(-variable_params[i, "detect_k"] * (mob.covtab@times - variable_params[i, "detect_mid"])))
                 )[-seq(1, min(pomp_data[!is.na(pomp_data$cases), ]$day) - 1)]))

betat  <- rbind(betat, betat.t)
Reff   <- rbind(Reff, Reff.t)
detect <- rbind(detect, detect.t)

## Stich together output
if (i == 1) {
  SEIR.sim.f <- SEIR.sim  
} else {
  SEIR.sim.f <- rbind(SEIR.sim.f, SEIR.sim) 
}  

## Keep track of progress
if (((i / 20) %% 1) == 0) {
  print(paste(round(i / nrow(variable_params), 2)*100, "% Complete", sep = ""))
}

}

####
## Summary for plotting
####

SEIR.sim.f.t <- SEIR.sim.f %>% 
  dplyr::select(date, day, .id
  #  , paramset
    , any_of(plot_vars))  %>%
  pivot_longer(any_of(plot_vars))

SEIR.sim.f.ci <- SEIR.sim.f.t %>% 
# fig4_data <- fig4_data %>% 
  group_by(date
  #  , paramset
    , name
#    , intervention
    ) %>% 
  summarize(
      lwr = quantile(value, 0 + ci.stoch, na.rm = T)
    , mid = {
      if (plot.median) {
      quantile(value, c(0.500), na.rm = T)  
      } else {
      mean(value, na.rm = T) 
      }
    }
    , upr = quantile(value, 1 - ci.stoch, na.rm = T)) 

# SEIR.sim.f %>% filter(.id == 10) %>% ggplot() + geom_line(aes(date, betat))
# SEIR.sim.f %>% filter(.id == 10) %>% ggplot() + geom_line(aes(date, I)) + geom_hline(yintercept = 1) + scale_y_log10()
# fig4_data %>% filter(.id != "median", name == "betat", intervention == "Continue Shelter in Place") %>% ggplot() + geom_line(aes(date, value))
# fig4_data %>% filter(.id == 4, name == "betat", intervention == "Continue Shelter in Place") %>% ggplot() + geom_line(aes(date, value))
# fig4_data %>% filter(.id != "median", name == "I", intervention == "Continue Shelter in Place") %>% ggplot() + geom_line(aes(date, value)) + scale_y_log10()
# fig4_data %>% filter(.id == 10, name == "I", intervention == "Continue Shelter in Place") %>% ggplot() + geom_line(aes(date, value)) + scale_y_log10()
# fig4_data %>% filter(.id != "median", name == "I", intervention == "Continue Shelter in Place") %>% ggplot() +
#   geom_line(aes(date, value, group = .id), alpha = 0.1) + scale_y_log10()
