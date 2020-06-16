source("ggplot_theme.R")
source("needed_packages.R")

#####
## Figure 1: model fits ---
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
    rds.name = "output/Los Angeles_independent_theta_200_2020-06-15.rds",
#   rds.name = "output/Los Angeles_constrained_theta_200_2020-06-14.rds",
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
#   rds.name = "output/Miami-Dade_constrained_theta_200_2020-06-14.rds",
    con_theta = F
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

## LA with same axis
{
fig1_data %>% 
 # filter(county != "Los Angeles") %>%
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
}

## LA with a different axis
{
fig1.1 <- fig1_data %>% 
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
    , strip.text.x = element_text(size = 11)
    , strip.text.y = element_text(size = 16)
    , axis.text.x = element_text(size = 12)
    , plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) +
  xlab("Date")

fig1.2 <- fig1_data %>% 
  filter(county == "Los Angeles") %>%
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
    , plot.margin = unit(c(0.25,0.25,0.25,-0.25), "cm")) +
  xlab("Date")

gridExtra::grid.arrange(fig1.1, fig1.2, ncol = 2, widths = c(4, 1.2))
}

#####
## Figure 2: Interventions ---
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

  
