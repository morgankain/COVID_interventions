source("ggplot_theme.R")
source("needed_packages.R")

#####
## Figure 1: model fits ----
#####

counties_list = {list(
  SC = list(
    focal.county = "Santa Clara",
    focal.state = "California",
    focal.state_abbr = "CA",
    rds.name = "./output/Santa Clara_TRUE_FALSE_0_2020-06-15_cont_temp_SC_ind_theta_continue.rds",
    con_theta = F
  ),
  LA = list(
    focal.county = "Los Angeles",
    focal.state = "California",
    focal.state_abbr = "CA",
    rds.name = "output/Los Angeles_independent_theta_200_2020-06-15.rds",
    con_theta = F
  ),
  KC = list(
    focal.county = "King",
    focal.state = "Washington",
    focal.state_abbr = "WA",
    rds.name = "output/King_TRUE_FALSE_0_2020-06-15_cont_temp_KC_ind_theta_continue.rds",
    con_theta = F
  ),
  FC = list(
    focal.county = "Fulton",
    focal.state = "Georgia",
    focal.state_abbr = "GA",
    rds.name = "output/Fulton_independent_theta_200_2020-06-15.rds",
    con_theta = F
  ),
  MD = list(
    focal.county = "Miami-Dade",
    focal.state = "Florida",
    focal.state_abbr = "FL",
    rds.name = "output/Miami-Dade_independent_theta_200_2020-06-15.rds",
    con_theta = F
  ),
  CC = list(
    focal.county = "Contra Costa",
    focal.state = "California",
    focal.state_abbr = "CA",
    rds.name = "output/Contra Costa_independent_theta_200_2020-06-15.rds",
    con_theta = F
  )
)}[-6]
int_vars <- list(
  counter.factual    = FALSE,
  int.movement       = "post",
  int.type           = "none",
  int.init           = "2020-06-08",
  int.end            = "2020-08-01", ## Doesnt matter, ignored for all int.type != "tail"
  sim_title          = "Reality"                                        
)
# set parameters
source("COVID_simulate_cont_params.R")
loglik.max     <- T
loglik.thresh  <- 2
ci.stoc        <- 0.025
ci.epidemic    <- T
nsim           <- 500
plot_vars      <- c("cases", "deaths")

# run the simulations for all locations
fig1_data <- plyr::adply(1:length(counties_list), 1, 
      function(i){
        with(c(counties_list[[i]], int_vars), {
        source("./COVID_simulate_cont.R", local = T)
        SEIR.sim.f.ci %<>% 
          full_join(county.data %>%
                      arrange(date) %>%
                      mutate(D = cumsum(ifelse(is.na(deaths), 0, deaths))*
                               ifelse(is.na(deaths), NA, 1)) %>%
                      select(date, any_of(plot_vars)) %>%
                      pivot_longer(any_of(plot_vars), values_to = "data")) %>% 
#          rbind(Reff %>% mutate(name = "Reff", mid = Reff, lwr = NA, upr = NA,
#                                data = NA) %>% 
#                  select(-Reff)) %>% 
          mutate(intervention = sim_title,
                 county = focal.county,
                 state = focal.state_abbr)
        return(SEIR.sim.f.ci)
})}, .id = NULL)

# plot the fits and data
fig1_colors = c("#FF0000",
                "#00A08A",
                "#F2AD00",
                "#F98400",
                "#5BBCD6")[c(1,4,5,2,3)]

fig1_data$name <- sapply(fig1_data$name, simpleCap)

fig1_data %>% 
  filter(county != "Los Angeles") %>%
  mutate(county = paste0(county, " County,", state)) %>%
  filter(date < as.Date("2020-06-08")) %>%
  filter(date >= as.Date("2020-02-10")) %>%
  filter(name != "Reff") %>%
  group_by(county) %>% 
  mutate(nparams = 0.5/length(unique(paramset))) %>% 
  ggplot(aes(x = date, y = mid, ymin = lwr, ymax = upr, 
             fill  = county, color = county, 
             group = interaction(county,paramset))) +
  geom_ribbon(aes(alpha = I(nparams)), colour = NA) +
  geom_line() +
  geom_point(aes(x = date, y = data), 
             color = "black", size = 0.75) + 
  scale_y_continuous(trans = "sqrt") + 
  scale_fill_manual(guide = F, values = fig1_colors) +
  scale_color_manual(guide = F, values = fig1_colors) +
  facet_grid(name ~ county, scales = "free_y", switch = "y")  +
  ylab("") + 
  theme(
    strip.background = element_blank()
    , strip.placement = "outside"
    , strip.text = element_text(size = 16)
    , axis.text.x = element_text(size = 12)) +
  xlab("Date")

fig1_data %>% 
  filter(date < as.Date("2020-06-12")) %>%
  filter(date >= as.Date("2020-02-10")) %>%
  filter(name == "Reff") %>% 
  ggplot(aes(x = date, y = mid, color = county,
             group = interaction(paramset, county))) + 
  geom_line() + 
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(values = fig1_colors) + 
  ylab("Reff")

#####
## Figure 2: Interventions ----
#####

counties_list <- {
  SC = list(
  list(
    focal.county = "Santa Clara",
    focal.state = "California",
    focal.state_abbr = "CA",
    rds.name = "./output/Santa Clara_TRUE_FALSE_0_2020-06-12_cont_temp_SC_ind_theta.rds",
    con_theta = F
  ),
  LA = list(
    focal.county = "Los Angeles",
    focal.state = "California",
    focal.state_abbr = "CA",
    rds.name = "output/Los Angeles_TRUE_FALSE_0_2020-06-12_cont_temp_LA_ind_theta.rds",
    con_theta = F
  )
)
}

int_vars <- {
  lift_all = list(
  list(
  counter.factual    = FALSE,
  int.movement       = "pre",
  int.type           = "none",
  int.init           = "2020-06-08",
  int.end            = "2020-08-01", ## Doesnt matter, ignored for all int.type != "tail"
  sim_title          = "Lift all interventions"                                    
),
  continue_SIP = list(
  counter.factual    = FALSE,
  int.movement       = "post",
  int.type           = "none",
  int.init           = "2020-06-08",
  int.end            = "2020-08-01", ## Doesnt matter, ignored for all int.type != "tail"
  sim_title          = "Continue Shelter in Place"                                        
),
  inf_iso = list(
  counter.factual    = FALSE,
  int.movement       = "mid",
  int.type           = "inf_iso",
  int.init           = "2020-06-08",
  int.end            = "2020-08-01", ## Doesnt matter, ignored for all int.type != "tail"
  sim_title          = "Infected Isolation",
  iso_mild_level     = 0.1,
iso_severe_level     = 0.1
),
  tail_chop = list(
  counter.factual    = FALSE,
  int.movement       = "mid",
  int.type           = "tail",
  int.init           = "2020-06-08",
  int.end            = "2020-08-01", ## Doesnt matter, ignored for all int.type != "tail"
  sim_title          = "Superspreading Averted",
  int.beta_catch     = 0.10,
  int.catch_eff      = 1
)
)
}

source("COVID_simulate_cont_params.R")

loglik.max     <- TRUE
plot_vars      <- c("cases", "deaths")
ci.stoch       <- 0.05
ci.epidemic    <- T
nsim           <- 20

fig2_data <- adply(1:length(int_vars), 1, 
      function(j) {
        adply(1:length(counties_list), 1,
      function(i) {
        with(c(counties_list[[i]], int_vars[[j]]), {
        source("./COVID_simulate_cont.R", local = T)
        SEIR.sim.f.ci %<>% 
          full_join(county.data %>%
                      arrange(date) %>%
                      mutate(D = cumsum(ifelse(is.na(deaths), 0, deaths))*
                               ifelse(is.na(deaths), NA, 1)) %>%
                      select(date, any_of(plot_vars)) %>%
                      pivot_longer(any_of(plot_vars), values_to = "data")) %>% 
          mutate(intervention = sim_title,
                 county = focal.county,
                 state = focal.state_abbr) 
        return(SEIR.sim.f.ci)
        
        }
        )
      }
        )
      }
        , .id = NULL)

fig2_data$name         <- sapply(fig2_data$name, simpleCap)
fig2_data$intervention <- factor(fig2_data$intervention, levels = unique(fig2_data$intervention))

fig2_colors <- c("#D67236", "dodgerblue4", "#0b775e", "magenta4")

fig2_data %>% 
  filter(date >= as.Date("2020-02-10")) %>%
  filter(name != "Reff") %>%
  ggplot(aes(x = date, y = mid, ymin = lwr, ymax = upr, fill  = intervention, color = intervention)) +
  geom_ribbon(data = (fig2_data %>% filter(date >= as.Date("2020-06-08")))
   , alpha = 0.50, colour = NA) +
  geom_line(data = (fig2_data %>% filter(date >= as.Date("2020-06-08")))
   ) +
  geom_ribbon(data = (fig2_data %>% filter(date <= as.Date("2020-06-08"), intervention == "Lift all interventions"))
  , alpha = 0.50, colour = NA, fill = "black") +
  geom_line(data = (fig2_data %>% filter(date <= as.Date("2020-06-08"), intervention == "Lift all interventions"))
  , colour = "black") +
  geom_vline(xintercept = as.Date("2020-06-08"), linetype = "dashed", lwd = 0.5) + 
  geom_point(aes(x = date, y = data), 
             color = "black", size = 1) + 
  scale_y_continuous(trans = "log10", labels = trans_format('log10', math_format(10^.x))) +
  scale_fill_manual(values = fig2_colors, name = "Intervention") +
  scale_color_manual(values = fig2_colors, name = "Intervention") +
  facet_grid(name ~ county, scales = "free_y", switch = "y")  +
  ylab("") + 
  theme(
    strip.background = element_blank()
    , legend.position  = c(0.650, 0.885)
    , legend.background = element_blank()
    , strip.placement = "outside"
    , strip.text = element_text(size = 16)
    , axis.text.x = element_text(size = 14)) +
  xlab("Date")

#####
## Figure 3: SIP and truncation combinations ----
#####  

# here's function to calculate SIP prop needed to maintain R=1 with a given tail truncation
# works with vector catch_eff and beta_catch. output is beta_catch X catch_eff in dimensions
sip_trunc_combns = function(beta_catch, beta_catch_type = "pct", 
                            catch_eff, k, beta0, beta_min, fixed_params){
  if(beta_catch_type == "pct"){
    upper = qgamma(beta_catch, k, scale = beta0/k)
  } else{
    upper = beta_catch
  }
  
  # from fixed params calculated average duration of infection
  d <- {(                
    ## proportion * time in asymptomatic
    fixed_params["alpha"] * fixed_params["Ca"] * (1/fixed_params["lambda_a"]) +                  
      ## proportion * time in mildly symptomatic
      (1 - fixed_params["alpha"]) * fixed_params["mu"] * ((1/fixed_params["lambda_p"]) + (1/fixed_params["lambda_m"])) +    
      ## proportion * time in severely symptomatic
      (1 - fixed_params["alpha"]) * (1 - fixed_params["mu"]) * ((1/fixed_params["lambda_p"]) + (1/fixed_params["lambda_s"]))      
  )}
  
  # calculate expected value of truncated gamma dist when truncation with 100% efficacy
  E_trunc <- beta0/k*pgamma(upper*k/beta0, k+1)*gamma(k+1)/(pgamma(upper*k/beta0, k)*gamma(k))
  # out <- (log(catch_eff*E_trunc + (1- catch_eff)*beta0) + log(d))/-log(beta_min)
  out <- -(log(outer(E_trunc, catch_eff, "*") + 
                 matrix(rep(outer(beta0, 1- catch_eff,"*"), length(beta_catch)), 
                        byrow = T, nrow = length(beta_catch))) + log(d))/log(beta_min)
  return(out)
}

beta_catch_vals <- seq(0.5, 1, by = 0.01)
catch_eff_vals <- seq(0.5, 1, by = 0.1)
log_lik_thresh <- 0
fig3_data <- adply(1:length(counties_list), 1, 
                   function(i){
                     fixed_params = readRDS(counties_list[[i]]$rds.name)$fixed_params
                     variable_params = readRDS(counties_list[[i]]$rds.name)$variable_params %>% 
                       filter(log_lik < 0) %>% 
                       filter(log_lik >= max(log_lik) - log_lik_thresh) %>%
                       arrange(desc(log_lik))
                     # now loop over all the top fits
                     adply(1:nrow(variable_params), 1, function(j){
                       params_row = unlist(variable_params[j,])
                       data.frame(
                         beta_catch_pct = beta_catch_vals,
                         sip = sip_trunc_combns(beta_catch_vals, 
                                                "pct",
                                                catch_eff_vals, 
                                                k = fixed_params["beta0_sigma"],
                                                beta0 = params_row["beta0est"],
                                                beta_min = params_row["beta_min"],
                                                fixed_params = c(fixed_params, params_row)) %>% 
                           magrittr::set_colnames(paste0("catch_eff_", catch_eff_vals)),
                         county = counties_list[[i]]$focal.county,
                         paramset = params_row["paramset"]) %>% return
                     }, .id = NULL) %>% return
                   }, .id = NULL) %>% 
  # make catch efficiency a column and pivot data longer
  pivot_longer(starts_with("sip"), names_to = "catch_eff", 
               names_prefix = "sip.catch_eff_", values_to = "sip")

fig3_data %>% 
  filter(sip > 0) %>%
  mutate(., catch_eff_label = factor(paste0(as.numeric(catch_eff)*100, "% efficiency"),
                                     levels = paste0(as.numeric(unique(pull(., catch_eff)))*100, "% efficiency"))) %>%
  ggplot(aes(x = 1- beta_catch_pct, y = sip, color = county, 
             group = interaction(paramset,county))) + 
  geom_line(alpha =1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  xlab("Percentile of superspreading averted") +
  ylab("Proportion sheltering-in-place") + 
  scale_color_manual(values = fig1_colors) +
  facet_wrap(~catch_eff_label)
