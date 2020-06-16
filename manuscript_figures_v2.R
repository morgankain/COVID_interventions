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
    rds.name = "./output/Santa Clara_independent_theta_200_2020-06-15.rds",
    con_theta = F
  ),
  LA = list(
    focal.county = "Los Angeles",
    focal.state = "California",
    focal.state_abbr = "CA",
#   rds.name = "output/Los Angeles_independent_theta_200_2020-06-15.rds",
    rds.name = "output/Los Angeles_constrained_theta_200_2020-06-14.rds",
    con_theta = T
  ),
  KC = list(
    focal.county = "King",
    focal.state = "Washington",
    focal.state_abbr = "WA",
    rds.name = "output/King_independent_theta_200_2020-06-16.rds",
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
  # rds.name = "output/Miami-Dade_independent_theta_200_2020-06-15.rds",
    rds.name = "output/Miami-Dade_constrained_theta_200_2020-06-14.rds",
    con_theta = T
  )#,
#  CC = list(
#    focal.county = "Contra Costa",
#    focal.state = "California",
#    focal.state_abbr = "CA",
#    rds.name = "output/Contra Costa_independent_theta_200_2020-06-15.rds",
#    con_theta = F
#  )
)}
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
loglik.max     <- F
loglik.num     <- 10
loglik.thresh  <- 2
ci.stoc        <- 0.025
ci.epidemic    <- T
nsim           <- 500
plot_vars      <- c("cases", "deaths")
plot.median    <- F

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
         rbind(Reff %>% mutate(name = "Reff", mid = Reff, lwr = NA, upr = NA,
                               data = NA) %>%
                 select(-Reff)) %>%
          rbind(detect %>% mutate(name = "detect", mid = detect, lwr = NA, upr = NA,
                                data = NA) %>%
                  select(-detect)) %>%
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

## LA with same axis
{
fig1_data %>% 
 # filter(county != "Los Angeles") %>%
  mutate(county = paste0(county, " County, ", state)) %>%
  filter(date < as.Date("2020-06-08")) %>%
  filter(date >= as.Date("2020-02-10")) %>%
  # filter(!(name %in% c("Reff", "Detect"))) %>%
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
  facet_grid(name ~ county, scales = "free", switch = "y")  +
  ylab("") + 
  theme(
    strip.background = element_blank()
    , strip.placement = "outside"
    , strip.text = element_text(size = 16)
    , axis.text.x = element_text(size = 12)) +
  xlab("Date")
}

## LA with a different axis
{
fig1.1 <- fig1_data %>% 
  filter(county != "Los Angeles") %>%
  mutate(county = paste0(county, " County, ", state)) %>%
  filter(date < as.Date("2020-06-08")) %>%
  filter(date >= as.Date("2020-02-10")) %>%
  filter(!(name %in% c("Reff", "Detect"))) %>%
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
    , strip.text.x = element_text(size = 11)
    , strip.text.y = element_text(size = 16)
    , axis.text.x = element_text(size = 12)
    , axis.title.y = element_text(margin = margin(l = 1))
    , plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")
    ) + 
    xlab("") +
    labs(tag = "C.") +
    theme(plot.tag.position = "topleft"
          , plot.tag = element_text(face = "bold" 
                                   , margin = margin(b = -16, 
                                                     r = -45)))

fig1.2 <- fig1_data %>% 
  filter(county == "Los Angeles") %>%
  mutate(county = paste0(county, " County, ", state)) %>%
  filter(date < as.Date("2020-06-08")) %>%
  filter(date >= as.Date("2020-02-10")) %>%
  filter(!(name %in% c("Reff", "Detect"))) %>%
  group_by(county) %>% 
  mutate(nparams = 0.5/length(unique(paramset))) %>% 
  ggplot(aes(x = date, y = mid, ymin = lwr, ymax = upr, 
             fill  = county, color = county, 
             group = interaction(county,paramset))) +
  geom_ribbon(aes(alpha = I(nparams)), colour = NA) +
  geom_line() +
  geom_point(aes(x = date, y = data), 
             color = "black", size = 0.75) + 
  scale_y_continuous(trans = "sqrt", position = "right") + 
  scale_fill_manual(guide = F, values = fig1_colors[5]) +
  scale_color_manual(guide = F, values = fig1_colors[5]) +
  facet_grid(name ~ county, scales = "free_y", switch = "y")  +
  ylab("") + 
  theme(
    strip.background = element_blank()
    , strip.placement = "outside"
    , strip.text.y = element_blank()
    , strip.text.x = element_text(size = 11)
    , axis.text.x = element_text(size = 12)
    , plot.margin = unit(c(0.25,0.25,0.25,0), "cm")) +
  xlab("") 
  # xlab("Date")

fig1.3 <- fig1_data %>% 
  filter(date < as.Date("2020-06-08")) %>%
  filter(date >= as.Date("2020-02-10")) %>%
  filter(name == "Reff") %>%
  # ugly way to calculate 7 day smoothing 
  arrange(county, paramset, date) %>%
  group_by(county, paramset) %>%
  mutate(
    lag3 = lag(mid, 3),
    lag2 = lag(mid, 2),
    lag1 = lag(mid, 1),
    lead1 = lead(mid, 1),
    lead2 = lead(mid, 2),
    lead3 = lead(mid, 3)
  ) %>% 
  rowwise() %>% 
  mutate(mid = mean(c(lag1, lag2, lag3, mid, lead1, lead2, lead3), na.rm = T)) %>% 
  select(-starts_with("lag"), -starts_with("lead")) %>% 
  mutate(paramset = as.character(paramset)) %>%
  group_by(county, state, date) %>% 
  {rbind(., 
         dplyr::summarise(., 
                          med = mean(mid), 
                          lwr = min(mid),
                          upr = max(mid),
                          .groups = "drop") %>% 
           rename(mid = med) %>% 
           mutate(paramset = "summary",
                  name = "Detect",
                  intervention = "Reality",
                  data = NA))} %>% 
  mutate(width = ifelse(paramset == "summary", 1.5, 0.25),
         alpha = ifelse(paramset == "summary", 1, 0.2)) %>%
  ggplot(aes(x = date, y = mid, ymin = lwr, ymax = upr, 
             fill  = county, color = county, 
             group = interaction(county,paramset))) +
  geom_ribbon(alpha = 0.25, color = NA) +
  geom_line(aes(size = I(width), alpha = I(alpha))) +
  scale_color_manual(guide= F, values = fig1_colors[c(1, 2, 5, 3, 4)]) +
  scale_fill_manual(guide= F, values = fig1_colors[c(1, 2, 5, 3, 4)]) +
  geom_hline(yintercept  = 1, linetype = "dashed") + 
  ylab(expression("R"[e])) + 
  xlab("") +
  scale_y_continuous(position = "left") + 
  ggtitle("A. Reproduction number")+
  theme(plot.title.position = "plot",
        plot.title = element_text(face = "bold"))

fig1.4 <-
  fig1_data %>%  
  filter(date < as.Date("2020-06-08")) %>%
  filter(date >= as.Date("2020-02-10")) %>%
  filter(name == "Detect") %>%
  mutate(paramset = as.character(paramset)) %>%
  group_by(county, state, date) %>% 
  {rbind(., 
         dplyr::summarise(., 
                          med = mean(mid), 
                          lwr = min(mid),
                          upr = max(mid),
                          .groups = "drop") %>% 
           rename(mid = med) %>% 
           mutate(paramset = "summary",
                  name = "Detect",
                  intervention = "Reality",
                  data = NA))} %>% 
  mutate(width = ifelse(paramset == "summary", 1.5, 0.25),
         alpha = ifelse(paramset == "summary", 1, 0.2)) %>%
  ggplot(aes(x = date, y = mid, ymin = lwr, ymax = upr, 
             group = interaction(county, paramset),
             fill  = county, color = county)) +
  geom_ribbon(alpha = 0.25, color = NA) +
  geom_line(aes(size = I(width), alpha = I(alpha))) +
  scale_color_manual(guide= F, values = fig1_colors[c(1, 2, 5, 3, 4)]) +
  scale_fill_manual(guide= F, values = fig1_colors[c(1, 2, 5, 3, 4)]) +
  # geom_hline(yintercept  = 1, linetype = "dashed") + 
  ylab("Detection probability") + 
  xlab("") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     position = "right") +
  ggtitle("B. Syptomatic detection")+
  theme(plot.title.position = "plot",
        plot.title = element_text(face = "bold"))
  
# gridExtra::grid.arrange(fig1.1, fig1.2, ncol = 2, widths = c(4, 1.2))
fig1 <- gridExtra::arrangeGrob(fig1.1, fig1.2, fig1.3,fig1.4,  
                       layout_matrix = matrix(c(3, 4, 4, 1, 1, 2), 
                                              byrow  = T, nrow = 2),
                       widths = c(2.6, 1.3, 1.3), heights = c(1.5, 2.5)) 
ggsave("figures/Manuscript2/Figure1.pdf", 
       fig1, 
       device = "pdf",
       width = 12,
       height = 8)
}

#####
## Figure 2: Interventions ----
#####

counties_list <- {
  list(
  SC =   list(
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
list(
  lift_all = list(
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
nsim           <- 100

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

# check calculation for truncated gamma distribution mean
## truncated gamma distribution 
# rtgamma = function(n, shape, scale, lower, upper, limits = "pct"){
#   if(limits == "pct"){
#     upper <- qgamma(upper, shape, scale = scale)
#     lower <- qgamma(lower, shape, scale = scale)
#   }
#   out <- rgamma(n, shape, scale = scale)
#   trunc <- which(out < lower | out > upper)
#   while(length(trunc) > 0){
#     out[trunc] <- rgamma(length(trunc), shape, scale = scale)
#     trunc <- which(out < lower | out > upper)
#   }
#   return(out)
# }
# 
# test <- adply(seq(1, 5, by = 0.1), 1, 
#               function(beta0){
#                 k = 0.16
#                 upper = 3
#                 sim_est = rtgamma(10000, 
#                                   shape = k,
#                                   scale = beta0/k, 
#                                   0, upper, limits = "abs") %>% mean
#                 alg_est = beta0/k*pgamma(upper*k/beta0, k+1)*gamma(k+1)/(pgamma(upper*k/beta0, k)*gamma(k))
#                 return(data.frame(beta0 = beta0, sim = sim_est, alg = alg_est))
#               })
# plot(test$sim, test$alg)
# abline(a = 0, b = 1, col = "red", lty = "dashed")

# here's function to calculate SIP prop needed to maintain R=1 with a given tail truncation
# works with vector catch_eff and beta_catch. output is beta_catch X catch_eff in dimensions
sip_trunc_combns = function(beta_catch, beta_catch_type = "pct", 
                            catch_eff, k, beta0, beta_min, fixed_params, desired_R){
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
  
  # expected value is weighted sum of truncated and untruncated distributions depdending on efficacy
  out <- -(log(outer(E_trunc, catch_eff, "*") + 
                 matrix(rep(outer(beta0, 1- catch_eff,"*"), length(beta_catch)), 
                        byrow = T, nrow = length(beta_catch))) + log(d) - log(desired_R))/log(beta_min)
  return(out)
}

counties_list = {list(
  SC = list(
    focal.county = "Santa Clara",
    focal.state = "California",
    focal.state_abbr = "CA",
    rds.name = "./output/Santa Clara_independent_theta_200_2020-06-15.rds",
    con_theta = F
  ),
  LA = list(
    focal.county = "Los Angeles",
    focal.state = "California",
    focal.state_abbr = "CA",
    # rds.name = "output/Los Angeles_independent_theta_200_2020-06-15.rds",
    rds.name = "output/Los Angeles_constrained_theta_200_2020-06-14.rds",
    con_theta = T
  ),
  KC = list(
    focal.county = "King",
    focal.state = "Washington",
    focal.state_abbr = "WA",
    rds.name = "output/King_independent_theta_200_2020-06-16.rds",
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
    # rds.name = "output/Miami-Dade_independent_theta_200_2020-06-15.rds",
    rds.name = "output/Miami-Dade_constrained_theta_200_2020-06-14.rds",
    con_theta = T
  )#,
  #  CC = list(
  #    focal.county = "Contra Costa",
  #    focal.state = "California",
  #    focal.state_abbr = "CA",
  #    rds.name = "output/Contra Costa_independent_theta_200_2020-06-15.rds",
  #    con_theta = F
  #  )
)}

beta_catch_vals <- seq(0.5, 1, by = 0.01)
catch_eff_vals <- seq(0.5, 1, by = 0.25)
loglik.max     <- F
loglik.num     <- 10
loglik.thresh  <- 2
fig3_data <- adply(1:length(counties_list), 1, 
                   function(i){
                     fixed_params = readRDS(counties_list[[i]]$rds.name)$fixed_params
                     variable_params = readRDS(counties_list[[i]]$rds.name)$variable_params 
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
                     mobility <- with(counties_list[[i]], {
                       if(focal.county == "Fulton"){
                       mobility <- readRDS(mobility.file) %>% 
                         dplyr::filter((county_name == focal.county | county_name == "DeKalb") & (state_abbr == focal.state_abbr)) %>%
                         dplyr::group_by(datestr) %>%
                         dplyr::summarize(sip_prop = mean(sip_prop)) %>%
                         dplyr::select(datestr, sip_prop) %>% 
                         dplyr::filter(!is.na(sip_prop)) 
                       } else{
                         mobility <- readRDS(mobility.file) %>% 
                           dplyr::filter(county_name == focal.county & state_abbr == focal.state_abbr)  %>%
                           dplyr::select(datestr, sip_prop) %>% 
                           dplyr::filter(!is.na(sip_prop)) %>%
                           mutate(day = as.numeric(datestr - as.Date("2019-12-31")))
                       }
                       mobility
                       })
                       
                     sip_start <- arrange(mobility, datestr) %>% head(7) %>% pull(sip_prop) %>% mean
                     sip_max <- arrange(mobility, datestr) %>% pull(sip_prop) %>% max
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
                                                fixed_params = c(fixed_params, params_row),
                                                desired_R    = 1.5) %>% 
                           magrittr::set_colnames(paste0("catch_eff_", catch_eff_vals)),
                         county = counties_list[[i]]$focal.county,
                         paramset = params_row["paramset"]) %>% return
                     }, .id = NULL) %>%
                       # make catch efficiency a column and pivot data longer
                       pivot_longer(starts_with("sip"), names_to = "catch_eff", 
                                    names_prefix = "sip.catch_eff_", values_to = "sip") %>%
                       # add two rows with data of sip observed
                       rbind(data.frame(
                         beta_catch_pct = 1,
                         county = counties_list[[i]]$focal.county,
                         paramset = "data",
                         sip = c(sip_start, sip_max),
                         catch_eff = 0
                       )) %>% 
                       return
                   }, .id = NULL) 
fig3_data %>% 
  filter(sip > 0.1) %>%
  group_by(county, catch_eff, beta_catch_pct) %>% 
  {rbind(.,
         filter(., paramset != "data") %>% 
           summarise(., 
                     lwr = min(sip),
                     upr = max(sip),
                     sip = mean(sip)) %>% 
           mutate(paramset = "summary"))} %>% 
  ungroup() %>%
  mutate(., catch_eff_label = factor(paste0(as.numeric(catch_eff)*100, "% efficiency"),
                                     levels = paste0(as.numeric(unique(pull(., catch_eff)))*100, "% efficiency"))) %>% 
  {ggplot(data = filter(., paramset == "summary"),
         mapping = aes(x = 1- beta_catch_pct, y = sip, 
             ymin = lwr, ymax = upr,
             color = county, 
             fill = county, 
             group = interaction(paramset,county))) + 
  geom_ribbon(alpha = 0.4, color = NA) +
  geom_line(alpha = 0.75, size = 1.5) + 
  geom_point(aes(shape = type, size = I(size)), 
             position = position_dodge(width = 0.05),
             alpha = 0.75, 
             data = filter(., paramset == "data") %>% 
               select(-catch_eff, -catch_eff_label) %>%
               group_by(county) %>% 
               mutate(type = ifelse(sip == max(sip), "max", "start")) %>%
               {right_join(., 
                           expand.grid(county = pull(., county) %>% unique, 
                                       catch_eff = catch_eff_vals))}  %>%
               mutate(., catch_eff_label = factor(paste0(as.numeric(catch_eff)*100, "% efficiency"),
                                                  levels = paste0(as.numeric(unique(pull(., catch_eff)))*100, "% efficiency"))) %>% 
               ungroup %>% 
               mutate(size = ifelse(type == "start", 2, 2))
             ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  xlab("Percentile of superspreading averted") +
  ylab("Proportion sheltering-in-place") + 
  scale_color_manual(values = fig1_colors) +
  scale_fill_manual(values = fig1_colors) +
      scale_shape_manual(guide = F, values = c(19, 17)) + #95)) + 
  facet_wrap(~catch_eff_label)}

#####
## Figure 4: Epidemic rebound when rare ----
#####  

## Have to run figure 3 in advance, which is annoying and likely should be fixed

counties_list <- {
  list(
  SC =   list(
    focal.county = "Santa Clara",
    focal.state = "California",
    focal.state_abbr = "CA",
    rds.name = "output/Santa Clara_independent_theta_200_2020-06-15.rds",
    con_theta = F
  )#,
#  KC = list(
#    focal.county = "King",
#    focal.state = "Washington",
#    focal.state_abbr = "WA",
#    rds.name = "output/King_independent_theta_200_2020-06-16.rds",
#    con_theta = F
#  )
)
}

int_vars <- {
list(
  lift_all = list(
  counter.factual    = FALSE,
  int.movement       = "mid",
  int.type           = "tail",
  int.init           = "2020-06-08",
  int.end            = "2020-08-01", ## Doesnt matter, ignored for all int.type != "tail"
  sim_title          = "Continue Shelter in Place",
  thresh_inf.val     = 5,
  int.catch_eff      = 1.0,
  int.beta_catch       = 0.10, #0.10,
  int.beta0_sigma      = 0.16,
  int.beta0_sigma_post = 0.16,
  int.beta_catch_post  = 0.05,
  int.catch_eff_post   = 0.0,
    ## Needed SIP scaling with "mid" to get R to the desired R
  int.sip_prop_scaling = check_R0(
    beta0est     = variable_params$beta0est
  , beta_min     = variable_params$beta_min
  , fixed_params = unlist(c(fixed_params, variable_params))
  , sd_strength  = 1
  , prop_S       = 1
  , desired_R    = 1.5) / 0.3463429
),
  continue_SIP = list(
  counter.factual    = FALSE,
  int.movement       = "mid",
  int.type           = "tail",
  int.init           = "2020-06-08",
  int.end            = "2020-08-01", ## Doesnt matter, ignored for all int.type != "tail"
  sim_title          = "Minimial Tail Cut",
  thresh_inf.val     = 5,
  int.catch_eff      = 1.0,
  int.beta_catch     = 0.10,
  int.beta0_sigma    = 0.16,
  int.beta0_sigma_post = 0.16,
  int.beta_catch_post  = 0.05,
  int.catch_eff_post   = 1.0,
  int.sip_prop_scaling = (fig3_data %>% filter(county == "Santa Clara", beta_catch_pct == 0.95, catch_eff == 1))$sip / 0.3463429                                       
),
  inf_iso = list(
  counter.factual    = FALSE,
  int.movement       = "mid",
  int.type           = "tail",
  int.init           = "2020-06-08",
  int.end            = "2020-08-01", ## Doesnt matter, ignored for all int.type != "tail"
  sim_title          = "Large Tail Cut",
  thresh_inf.val     = 5,
  int.catch_eff      = 1.0,
  int.beta_catch     = 0.10,
  int.beta0_sigma    = 0.16,
  int.beta0_sigma_post = 0.16,
  int.beta_catch_post  = 0.10,
  int.catch_eff_post   = 1.0,
  int.sip_prop_scaling = (fig3_data %>% filter(county == "Santa Clara", beta_catch_pct == 0.90, catch_eff == 1))$sip / 0.3463429
)
)}

source("epidemic_rebound/gamma_rebound_pomp.R")
source("ggplot_theme.R")
source("epidemic_rebound/gamma_rebound_params.R")

nsim               <- 500
thresh_inf.val     <- 10
sim_length         <- 275     ## How many days to run the simulation

fig4_data <- adply(1:length(int_vars), 1, 
      function(j) {
        adply(1:length(counties_list), 1,
      function(i) {
        with(c(counties_list[[i]], int_vars[[j]]), {
        source("epidemic_rebound/gamma_rebound.R", local = T)
        SEIR.sim.f.t %<>% 
          full_join(county.data %>%
                      arrange(date) %>%
                      mutate(D = cumsum(ifelse(is.na(deaths), 0, deaths))*
                               ifelse(is.na(deaths), NA, 1)) %>%
                      select(date, any_of(plot_vars)) %>%
                      pivot_longer(any_of(plot_vars), values_to = "data")) %>% 
          mutate(intervention = sim_title,
                 county = focal.county,
                 state = focal.state_abbr) 
    #    return(SEIR.sim.f.ci)
        return(SEIR.sim.f.t)
        
        }
        )
      }
        )
      }
        , .id = NULL)

fig4_data$name         <- sapply(fig4_data$name, simpleCap)
fig4_data$intervention <- factor(fig4_data$intervention, levels = unique(fig4_data$intervention))

fig4_colors <- c("#D67236", "dodgerblue4", "#0b775e", "magenta4")
fig4_data   <- fig4_data %>% filter(state == "CA")

check_date <- (fig4_data %>% filter(name == "Cases", mid == 0, intervention == "Continue Shelter in Place", date > "2020-05-01") %>%
  filter(date == min(date)))$date

fig4_data %>% 
  filter(date >= as.Date("2020-02-10")) %>%
  ggplot(aes(x = date, y = value
 #  , ymin = lwr, ymax = upr, fill  = intervention
    , color = intervention
    , group = .id)) +
#  geom_ribbon(data = (fig4_data %>% filter(date >= check_date))
#   , alpha = 0.50, colour = NA) +
#  geom_line(data = (fig4_data %>% filter(date >= check_date))
#   ) + 
  geom_line(data = (fig4_data %>% filter(date >= check_date, .id != "median"))
    , aes(group = interaction(.id, intervention))
    , lwd = 0.25, alpha = 0.15
   ) +
  geom_line(data = (fig4_data %>% filter(date >= check_date, .id == "median"))
   , lwd = 1.5) +
#  geom_ribbon(data = (fig4_data %>% filter(date <= check_date, intervention == "Continue Shelter in Place"))
#  , alpha = 0.50, colour = NA, fill = "black") +
#  geom_line(data = (fig4_data %>% filter(date <= check_date, intervention == "Continue Shelter in Place"))
#  , colour = "black") +
  geom_line(data = (fig4_data %>% filter(date <= check_date, intervention == "Continue Shelter in Place", .id != "median"))
    , lwd = 0.25, alpha = 0.15, colour = "black"
   ) +
  geom_line(data = (fig4_data %>% filter(date <= check_date, intervention == "Continue Shelter in Place", .id == "median"))
    , colour = "black"
    , lwd = 1.5) +
  geom_vline(xintercept = check_date, linetype = "dashed", lwd = 0.5) + 
  geom_point(aes(x = date, y = data), 
             color = "black", size = 1) + 
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = fig4_colors, name = "Intervention") +
  scale_color_manual(values = fig4_colors, name = "Intervention") +
  facet_grid(name ~ county, scales = "free_y", switch = "y")  +
  ylab("") + 
  theme(
    strip.background = element_blank()
    , legend.position  = c(0.850, 0.885)
    , legend.background = element_blank()
    , strip.placement = "outside"
    , strip.text = element_text(size = 16)
    , axis.text.x = element_text(size = 14)) +
  xlab("Date")

droplevels(fig4_data) %>% filter(state == "CA", .id != "median") %>% filter(name == "Cases") %>% filter(day > 260) %>% 
  group_by(.id, intervention) %>%
  summarize(new_cases = sum(value)) %>% 
  filter(new_cases == 0) %>% 
  group_by(intervention) %>%
  summarize(prop_extinct = n() / nsim)

