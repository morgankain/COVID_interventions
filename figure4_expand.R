source("needed_packages.R")
source("ggplot_theme.R")
source("epidemic_rebound/gamma_rebound_params.R")
source("epidemic_rebound/gamma_rebound_pomp2.R")

## truncated gamma distribution simulator
rtgamma = function(n, shape, scale, lower, upper, limits = "pct"){
  if(limits == "pct"){
    upper <- qgamma(upper, shape, scale = scale)
    lower <- qgamma(lower, shape, scale = scale)
  }
  out <- rgamma(n, shape, scale = scale)
  trunc <- which(out < lower | out > upper)
  while(length(trunc) > 0){
    out[trunc] <- rgamma(length(trunc), shape, scale = scale)
    trunc <- which(out < lower | out > upper)
  }
  return(out)
}


# here's function to calculate SIP prop needed to maintain specified R with a given tail truncation
# works with vector catch_eff and beta_catch. output is beta_catch X catch_eff in dimensions
sip_trunc_combns = function(beta_catch, 
                            beta_catch_type = "pct", 
                            catch_eff, 
                            k, beta0, beta_min, d, dt,
                            desired_R){
  shp = k*dt/d
  scl = beta0/k
  if(beta_catch_type == "pct"){
    upper = qgamma(1 - beta_catch, shape = shp, scale = scl)
  } else{
    upper = beta_catch
  }
  
  # calculate expected value of truncated gamma dist when truncation with 100% efficacy
  E_trunc <- beta0/k*pgamma(upper/scl, shp+1)*gamma(shp+1)/(pgamma(upper/scl, shp)*gamma(shp))
  
  # expected value is weighted sum of truncated and untruncated distributions depdending on efficacy
  out <- (log(outer(E_trunc, catch_eff, "*") + outer(rep(shp*scl, length(E_trunc)), 1-catch_eff, "*")) - 
            log(desired_R) - log(dt) + log(d))/-log(beta_min)
  return(out)
}

counties_list <- {
  list(
  FC = list(
    focal.county     = "Fulton"
  , focal.state      = "Georgia"
  , focal.state_abbr = "GA"
  , rds.name         = "output/Fulton_0_2020-06-25_cont_round2.rds"
  , con_theta        = F
  )
)
}

nsim            <- 200
plot_vars       <- c("cases", "deaths")
ci.stoch        <- 0.1
counter.factual <- FALSE
int.movement    <- c("post", "mid")
int.type        <- "tail"
int.init        <- "2020-07-01"
sim_end         <- "2020-12-01"
sim_title       <- "Minor tail chop"
thresh_inf.val  <- 5
int.beta_catch_type <- "pct"
int.catch_eff   <- 0.75
int.beta_catch  <- 0.005
int.beta0_k     <- 0.16
int.beta0_k_post <- 0.16
loglik.max      <- TRUE
ci.stoch        <- 0.10
ci.epidemic     <- T
plot.median     <- F
dt              <- 1/6

trunc_params <- expand.grid(
    beta_catch_vals = seq(0, 0.006, by = 0.002)
 # , catch_eff_vals = c(0.50, 0.75, 1.0)
   , catch_eff_vals = 0.75
)

fig4_data <- adply(1:length(counties_list), 1, 
  
                   function (i) {
                     
                     variable_params <- readRDS(counties_list[[i]]$rds.name)$mifs_local_v2 %>% 
                       filter(log_lik < 0) %>% 
                       filter(log_lik == max(log_lik))
                    
                     params_row     <- unlist(variable_params[1,])
                     
                     # now loop over all truncation interventions
                     adply(1:nrow(trunc_params), 1, function(j) {
                            with(c(counties_list[[i]]), {
                    
                       
                     SIP_post <- sip_trunc_combns(
                              beta_catch      = trunc_params[j, 1]
                            , beta_catch_type = "pct"
                            , catch_eff       = trunc_params[j, 2]
                            , k               = params_row["beta0_k"]
                            , beta0           = params_row["beta0"]
                            , beta_min        = params_row["beta_min"]
                            , d               = params_row["d"]
                            , dt              = dt
                            , desired_R       = 2.0)
                    
                     int.beta_catch_post  <- trunc_params[j, 1]
                     int.catch_eff_post   <- trunc_params[j, 2]
   
                     source("epidemic_rebound/gamma_rebound.R", local = T)
                     
                     SEIR.sim.f.t %<>% 
                       full_join(county.data %>%
                       arrange(date) %>%
                       mutate(D = cumsum(ifelse(is.na(deaths), 0, deaths))*
                               ifelse(is.na(deaths), NA, 1)) %>%
                       select(date, any_of(plot_vars)) %>%
                       pivot_longer(any_of(plot_vars), values_to = "data")) %>% 
                       mutate(
                           county          = focal.county
                         , state           = focal.state_abbr
                         , beta_catch      = trunc_params[j, 1]
                         , catch_eff       = trunc_params[j, 2]
                       )
                     
         return(SEIR.sim.f.t)
                            })
        
                     }, .id = NULL)
     }, .id = NULL)

###
## CI on the range of results after 3 weeks of hitting 0
###

fig4_data$name <- sapply(fig4_data$name, simpleCap)

fig4.top95 <- fig4_data %>% 
filter(.id != "median", name == "Cases") %>%
  mutate(reach1 = ifelse(value <= 1, 1, 0)) %>%
  filter(date > "2020-05-01") %>%
  group_by(.id, beta_catch, catch_eff) 

new_index   <- fig4.top95 %>% group_by(.id, beta_catch, catch_eff) %>% filter(reach1 == 1) %>% slice(1) %>%
  dplyr::select(day, .id, beta_catch, catch_eff) 

names(new_index)[1] <- "day2"

fig4.top95 <- fig4.top95 %>% left_join(., new_index)

###
## needed for prop extinct
###
fig4.prop_extinct <- fig4.top95 %>% group_by(.id, beta_catch, catch_eff) %>% filter(day >= (max(day) - 10)) %>% summarize(total_new_I = sum(value))
fig4.prop_extinct <- fig4.prop_extinct %>% group_by(beta_catch, catch_eff) %>% summarize(prop_extinct = length(which(total_new_I == 0))/nsim)

fig4.top95.gg <- fig4.top95 %>% filter(day > (day2 - 21)) %>% filter(day <= (day2 + 21))
fig4.top95    <- fig4.top95 %>% filter(day > day2) %>% filter(day <= (day2 + 21))

fig4.top95 <- fig4.top95 %>% group_by(beta_catch, catch_eff) %>%
  summarize(
      lwr = quantile(value, 0.025)
    , est = quantile(value, 0.500) 
    , upr = quantile(value, 0.975)
    )

###
## ggplot: Proportion that go extinct 
###

fig4.prop_extinct %>% 
  mutate(., catch_eff_label = factor(paste0(as.numeric(catch_eff)*100, "% efficiency"),
                                     levels = paste0(as.numeric(unique(pull(., catch_eff)))*100, "% efficiency"))) %>%
  ggplot(aes(x = 1- beta_catch, y = prop_extinct, color = as.factor(catch_eff))) + 
  geom_line(alpha = 1, lwd = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  xlab("Percentile of superspreading averted") +
  ylab("Proportion of epidemic simulations that go extinct") + 
  scale_color_brewer(
    palette="Dark2"
  , name   = "Probability each superspreading\nevent is averted")

###
## ggplot: top 95% for epidemic rebound
###

fig4.top95 %>% 
  mutate(., catch_eff_label = factor(paste0(as.numeric(catch_eff)*100, "% efficiency"),
                                     levels = paste0(as.numeric(unique(pull(., catch_eff)))*100, "% efficiency"))) %>%
  ggplot(aes(x = 1- beta_catch, y = upr, color = as.factor(catch_eff))) + 
  geom_line(alpha = 1, lwd = 1) +
 # scale_y_continuous(trans = "log10") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  xlab("Percentile of superspreading averted") +
  ylab("Upper 95% CI on new cases ") + 
  scale_color_brewer(
    palette="Dark2"
  , name   = "Probability each superspreading\nevent is averted")

###
## Summarize the original data to get envelopes on I
###

fig4_data.gg <- fig4_data %>% 
  filter(name == "I") %>%
  group_by(date, beta_catch, catch_eff) %>% 
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

fig4_data.gg <- fig4_data.gg %>% filter(catch_eff == 1)

check_date <- (fig4_data %>% filter(name == "I", .id == "median", catch_eff == "0.05", date > "2020-05-01") %>% 
    filter(value == min(value)) %>% 
    filter(date == min(date)))$date

fig4_data.gg$beta_catch <- as.factor(fig4_data.gg$beta_catch)

###
## ggplot: Envelopes on I
###

fig4_data.gg %>% filter(date > "2020-07-01") %>%
  ggplot(aes(x = date
    , y = mid
    , ymin = lwr, ymax = upr
    , fill  = beta_catch
    , color = beta_catch
    )) +
  geom_ribbon(alpha = 0.50, colour = NA) +
  geom_line() +   
  scale_fill_brewer( palette = "Dark2", name = "Percentile of superspreading averted") +
  scale_color_brewer(palette = "Dark2", name = "Percentile of superspreading averted") +
  ylab("") + 
  theme(
    strip.background = element_blank()
    , legend.position  = c(0.250, 0.885)
    , legend.background = element_blank()
    , strip.placement = "outside"
    , strip.text = element_text(size = 16)
    , axis.text.x = element_text(size = 14)) +
  xlab("Date") +
  scale_y_log10()

fig4.top95.gg %>% ggplot(aes(date, value)) + geom_line(aes(group = .id), alpha = 0.25) + 
  xlab("Date") + ylab("Total Infected") +
  facet_wrap(~beta_catch) + scale_y_log10()
 # geom_line(data = (fig4.top95.gg %>% filter(.id == 228)), lwd = 1, colour = "firebrick3") +
 # 
