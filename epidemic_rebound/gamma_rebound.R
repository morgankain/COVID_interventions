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

int.begin    <- as.numeric((as.Date(int.init) - date_origin))
int.duration <- min(as.Date(int.end), variable_params[i, ]$sim_start + sim_length) - as.Date(int.init)

int.phase1 <- nrow(mobility) + (int.begin - nrow(mobility))
int.phase2 <- as.numeric(int.duration)
int.phase3 <- sim_length - int.phase1 - int.phase2

## Covariate table with all of the intervention scenarios
{
mob.covtab <- covariate_table(
  ## First, can add whatever sip_prop into the future
  
   sip_prop      = {
     if (int.movement == "pre") {
     c(
      mobility$sip_prop
    , rep(mean(tail(mobility$sip_prop, 3)), int.phase1 - length(mobility$sip_prop))
    , rep(mean(head(mobility$sip_prop, 7)), sim_length - int.phase1)
      )       
     } else if (int.movement == "post") {
     c(
      mobility$sip_prop
    , rep(mean(tail(mobility$sip_prop, 3)), int.phase1 - length(mobility$sip_prop))
    , rep(mean(tail(mobility$sip_prop, 3)), sim_length - int.phase1)
      )
     } else if (int.movement == "mid") {
     c(
      mobility$sip_prop
    , rep(mean(tail(mobility$sip_prop, 3)), int.phase1 - length(mobility$sip_prop))
    , rep((mean(head(mobility$sip_prop, 7)) + mean(tail(mobility$sip_prop, 3))) / 2, sim_length - int.phase1)
      )         
     }
    
   }
  
 , order                 = "constant"
 , times                 = seq(1, sim_length, by = 1)
 , thresh_inf            = thresh_inf.val
 , check_int             = c(rep(1, int.phase1), rep(0, int.phase2 + int.phase3))
 , iso_mild_level        = iso_mild_level
 , iso_severe_level      = iso_severe_level
 , iso_mild_level_post   = iso_mild_level_post
 , iso_severe_level_post = iso_severe_level_post
  
  ## Next can add an intervention that stops large gatherings
  
 , intervention  = {
   if (int.type == "none") {
   c(
    rep(0, int.phase1)
  , rep(0, int.phase2)
  , rep(0, int.phase3)
  )     
   } else if (int.type == "tail") {
   c(
    rep(0, int.phase1)
  , rep(1, int.phase2)
  , rep(1, int.phase3)
   ) 
   } else if (int.type == "inf_iso") {
   c(
    rep(0, int.phase1)
  , rep(2, int.phase2)
  , rep(2, int.phase3)
   )   
  }
   }
    )
}

}

covid_mobility <- pomp(
   data       = county.data %>% select(day, cases, deaths)
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
fixed_params["beta0_k"]      <- int.beta0_k  ## heterogeneity value
fixed_params["beta_catch"]  <- {                ## beta0 values caught by intervention
  if(int.beta_catch_type == "pct") {
    max(qgamma(p = (1 - int.beta_catch), 
               shape = int.beta0_k, 
               scale = variable_params[i, "beta0"]/int.beta0_k),
        0.001)
  } else{
    max(int.beta_catch, 0.001) # crude way to ensure int.beta_catch > 0, note that even this nonzero low beta_catch will be problematic but will eventually finish
  }
}
fixed_params["catch_eff"]        <- int.catch_eff    ## beta0 values caught by intervention
fixed_params["sip_prop_post"]    <- SIP_post

fixed_params["beta_catch_post"]  <- {                ## beta0 values caught by intervention
  if(int.beta_catch_type == "pct") {
    max(qgamma(p = (1 - int.beta_catch_post), 
               shape = int.beta0_k_post, 
               scale = variable_params[i, "beta0"]/int.beta0_k_post),
        0.001)
  } else{
    max(int.beta_catch_post, 0.001) # crude way to ensure int.beta_catch > 0, note that even this nonzero low beta_catch will be problematic but will eventually finish
  }
}
fixed_params["catch_eff_post"]   <- int.catch_eff_post
fixed_params["beta0_k_post"]     <- int.beta0_k_post  ## heterogeneity value

print(int.type)
print(fixed_params["catch_eff_post"])
print(fixed_params["beta_catch_post"])

checktime <- system.time({

## Simulate from fits
SEIR.sim <- do.call(
  pomp::simulate
  , list(
    object   = covid_mobility
    , t0     = min(as.numeric(variable_params[i, ]$sim_start - date_origin), county.data$day)
  # , times  = 1:sim_length
    , times  = c(county.data$day, max(county.data$day):(max(county.data$day) + sim_length - length(county.data$day)))
    , params = c(fixed_params %>% t() %>% 
                   as.data.frame() %>% 
                   mutate(d  = variable_params[i, "alpha"] * lambda_a + 
                       (1 - variable_params[i, "alpha"]) * variable_params[i, "mu"] * (lambda_p + lambda_m) + 
                       (1 - variable_params[i, "alpha"])*(1 - variable_params[i, "mu"])*(lambda_p + lambda_s)
                   ) %>%  
                   # convert periods to rates
                   mutate(gamma    = -1/(nE*dt)*log(1-nE*dt/gamma),
                          lambda_a = -1/(nIa*dt)*log(1-nIa*dt/lambda_a),
                          lambda_p = -1/(nIp*dt)*log(1-nIp*dt/lambda_p), 
                          lambda_m = -1/(nIm*dt)*log(1-nIm*dt/lambda_m),
                          lambda_s = -1/(nIs*dt)*log(1-nIs*dt/lambda_s),
                          rho_r    = -1/dt*log(1-dt/rho_r)) %>%
                   unlist
      , c(
        beta0      = variable_params[i, "beta0"]
      , E_init     = variable_params[i, "E_init"]
      , detect_max = variable_params[i, "detect_max"]
      , detect_mid = variable_params[i, "detect_mid"]
      , detect_k   = variable_params[i, "detect_k"]
      , theta      = variable_params[i, "theta"]
      , theta2     = variable_params[i, "theta2"]
      , beta_min   = variable_params[i, "beta_min"]
      , Ca         = variable_params[i, "Ca"]
      , alpha      = variable_params[i, "alpha"]
      , delta      = variable_params[i, "delta"]
      , mu         = variable_params[i, "mu"]
      , E_init     = variable_params[i, "E_init"]
      , rho_d      = variable_params[i, "rho_d"]
      , detect_t0  = min(county.data[!is.na(county.data$cases), ]$day) - 1
        ))
    , nsim         = nsim
    , format       = "d"
    , include.data = F
    , seed         = 1001))
      
})[3]

print(checktime)
      
SEIR.sim <- SEIR.sim %>%
  mutate(
      date     = round(as.Date(day, origin = variable_params[i, ]$sim_start))
    , paramset = variable_params[i, ]$paramset)

if(ci.epidemic){
  epi_ids <- SEIR.sim %>% 
    group_by(.id) %>% 
    summarise(total_infect = max(D + R)) %>% 
    filter(total_infect > ci.epidemic_cut*ceiling(variable_params[i, "E_init"])) %>% 
    pull(.id)
  print(paste0("limiting to epidemics, including ", length(epi_ids), " simulations"))
  SEIR.sim %<>% filter(.id %in% epi_ids)
}

SEIR.sim %<>% {
  rbind(.,
        group_by(., day) %>%
          dplyr::select(-.id) %>%
          summarise_all(median) %>%
          mutate(.id = "median"))
  
}

SEIR.sim.s <- SEIR.sim %>% 
  dplyr::filter(.id != "median") %>%
  group_by(day) %>%
  summarize(S = mean(S), D = mean(D))

betat.t      <- variable_params[i, "beta0"] * exp(log(variable_params[i, "beta_min"])*mob.covtab@table[which(dimnames(mob.covtab@table)[[1]] == "sip_prop"), ])

#Reff.t       <- with(variable_params[i, ], covid_R0(
#   beta0est      = betat.t
# , fixed_params  = c(fixed_params, unlist(variable_params[i, ]))
# , sd_strength   = 1
# , prop_S        = SEIR.sim.s$S / (variable_params[i, "N"] - SEIR.sim.s$D)
#  )
#  )

#detect.t <- c(rep(0, mob.covtab@table[which(dimnames(mob.covtab@table)[[1]] == "detect_t0"), 1])
#  , (variable_params[i, "detect_max"] / 
#    (1 + exp(-variable_params[i, "detect_k"] * (mob.covtab@times - variable_params[i, "detect_mid"])))
#    )[-seq(1, mob.covtab@table[which(dimnames(mob.covtab@table)[[1]] == "detect_t0"), 1])]
#)

#betat.t <- data.frame(
#  paramset = variable_params[i, ]$paramset
#, date     = seq(variable_params[i, ]$sim_start
#, variable_params[i, ]$sim_start + (sim_length - 1), by = 1)
#, betat    = betat.t
#  )

#Reff.t <- data.frame(
#  paramset = variable_params[i, ]$paramset
#, date     = seq(variable_params[i, ]$sim_start
#, variable_params[i, ]$sim_start + (sim_length - 1), by = 1)
#, Reff     = Reff.t
#  )

#detect.t <- data.frame(
#  paramset = variable_params[i, ]$paramset
#, date     = seq(variable_params[i, ]$sim_start
#, variable_params[i, ]$sim_start + (sim_length - 1), by = 1)
#, detect   = detect.t
#  )

#if (sir_init.mid) {
#SEIR.sim.c <- left_join(SEIR.sim.c, detect.t[, c(2, 3)], by = "date")
#SEIR.sim.c <- SEIR.sim.c %>% mutate(cases = I_new_sympt * detect)
#}

# betat  <- rbind(betat, betat.t)
# Reff   <- rbind(Reff, Reff.t)
# detect <- rbind(detect, detect.t)

## Stich together output
#if (sir_init.mid) {
#if (i == 1) {
# SEIR.sim.f <- SEIR.sim.c 
#} else {
# SEIR.sim.f <- rbind(SEIR.sim.f, SEIR.sim.c) 
#}
#} else {
if (i == 1) {
 SEIR.sim.f <- SEIR.sim  
} else {
 SEIR.sim.f <- rbind(SEIR.sim.f, SEIR.sim) 
}  
#}

## Keep track of progress
if (((i / 20) %% 1) == 0) {
  print(paste(round(i / nrow(variable_params), 2)*100, "% Complete", sep = ""))
}

}

SEIR.sim.f <- SEIR.sim.f %>% filter(date < min(variable_params$sim_start + sim_length))

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
