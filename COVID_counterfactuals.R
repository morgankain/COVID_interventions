# set the base params that we will then update
source("./COVID_simulate_cont_params.R")

focal.county       <- "Contra Costa"
focal.state_abbr   <- "CA"
rds.name           <- "output/Contra_Costa_2020-06-10.Rds"

####
## Counterfactual: No SIP at all (continue pre-intervention beta through today)
####
counter.factual    <- TRUE
cf.type            <- "no_int"
int.movement       <- "pre"
int.type           <- "none"
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" ## Doesnt matter, ignored for all int.type != "tail"
sim_title          <- "No social distancing"
source("./COVID_simulate_cont.R")
SEIR.sim.f.ci.s1   <- SEIR.sim.f.ci %>% mutate(intervention = sim_title)
SEIR.sim.f.t.s1    <- SEIR.sim.f.t %>% mutate(intervention = sim_title)

####
## Counterfactual: Delay SIP by one week
####
counter.factual    <- TRUE
cf.type            <- "delay"
delay.time         <- 7
int.movement       <- "pre"
int.type           <- "none"
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" ## Doesnt matter, ignored for all int.type != "tail"
sim_title          <- "Social distancing delayed by one week"
source("./COVID_simulate_cont.R")
SEIR.sim.f.ci.s2   <- SEIR.sim.f.ci %>% mutate(intervention = sim_title)
SEIR.sim.f.t.s2    <- SEIR.sim.f.t %>% mutate(intervention = sim_title)

####
## Counterfactual: Lift SIP completely on May 1 (resume pre-intervention beta)
####
counter.factual    <- TRUE
cf.type            <- "may1"
int.movement       <- "pre"
int.type           <- "none"
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" ## Doesnt matter, ignored for all int.type != "tail"
sim_title          <- "Social distancing ended on May 1"
source("./COVID_simulate_cont.R")
SEIR.sim.f.ci.s3   <- SEIR.sim.f.ci %>% mutate(intervention = sim_title)
SEIR.sim.f.t.s3    <- SEIR.sim.f.t %>% mutate(intervention = sim_title)

####
## Counterfactual: Reality for comparison
####
counter.factual    <- FALSE
int.movement       <- "pre"
int.type           <- "none"
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" ## Doesnt matter, ignored for all int.type != "tail"
sim_title          <- "Reality"
source("./COVID_simulate_cont.R")
SEIR.sim.f.ci.s4   <- SEIR.sim.f.ci %>% mutate(intervention = sim_title)
SEIR.sim.f.t.s4    <- SEIR.sim.f.t %>% mutate(intervention = sim_title)

### Combine and plot all of these
SEIR.sim.f.t <- rbind(SEIR.sim.f.t.s4,
                      SEIR.sim.f.t.s1, 
                      SEIR.sim.f.t.s2, 
                      SEIR.sim.f.t.s3)
SEIR.sim.f.ci <- rbind(SEIR.sim.f.ci.s4,
                      SEIR.sim.f.ci.s1, 
                      SEIR.sim.f.ci.s2, 
                      SEIR.sim.f.ci.s3)

SEIR.sim.f.t  <- SEIR.sim.f.t %>% filter(date < max(county.data$date))
SEIR.sim.f.ci <- SEIR.sim.f.ci %>% filter(date < max(county.data$date))

fig4_colors <- c("#D67236", "grey30", "#0b775e", "magenta4")

{
(
gg <- ggplot() + 
  
  geom_ribbon(data = SEIR.sim.f.ci
  , aes(
    x = date
  , ymin = lwr
  , ymax = upr
  , fill = intervention
    ), alpha = 0.50) +
  
  geom_line(data = (SEIR.sim.f.t %>% filter(.id == "median"))
  , aes(
      x = date
    , y = value 
    , group  = interaction(.id, intervention)
    , colour = intervention
    ), alpha = 1, lwd = 1.5) + 

#  geom_ribbon(data = SEIR.sim.f.ci
#  , aes(
#    x = date
#  , ymin = lwr
#  , ymax = upr
#    ), fill = "grey30", alpha = 0.50) +
  
#  geom_line(data = (SEIR.sim.f.t %>% 
#                      {filter(., date <= as.Date("2020-06-08"), 
#                             .id == "median", 
#                             intervention == pull(., intervention) %>% 
#                               first)})
#  , aes(
#      x = date
#    , y = value 
#    , group  = .id
#    ), colour = "grey30", alpha = 1, lwd = 1.5) + 

    geom_point(data = county.data %>% 
                 mutate(D = cumsum(ifelse(is.na(deaths), 0, deaths))*
                          ifelse(is.na(deaths), NA, 1)) %>% 
                 select(date, any_of(plot_vars)) %>% 
                 pivot_longer(any_of(plot_vars))
    , aes(
      x = date
    , y = value)
      , colour = "black", lwd = 2) + 
  
#  geom_vline(xintercept = as.Date("2020-06-08"), linetype = "dashed", lwd = 0.65) +
  
  scale_colour_manual(
      values = fig4_colors
    , name   = "Scenario"
    # , labels = sim_names
    ) +
  scale_fill_manual(
      values = fig4_colors
    , name   = "Scenario"
    # , labels = sim_names
    ) +
  ylab(NULL) +
  scale_y_continuous(trans = "log10", 
                     breaks = c(1, 10, 100, 1000, 10000)) + 
  scale_x_date(labels = date_format("%b"), date_breaks = "2 month") +
  xlab("") +
  theme(
    axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 12)
  , legend.title     = element_text(size = 12)
  , legend.text      = element_text(margin = margin(t = 0.3, unit = 'cm'), size = 11)
  , plot.title       = element_text(size = 12)
  , legend.position  = c(0.245, 0.925)
  , axis.ticks.x     = element_blank()
  , strip.background = element_blank()
  , strip.placement = "outside"
  , strip.text = element_text(size = 18)) +
  facet_wrap(~ name, nrow = 2, scale = "free_y", 
             strip.position = "left", 
             labeller = as_labeller(c(cases = "Daily Cases", D = "Cumulative Deaths")))
)
}

## check point est for lives saved
SEIR.sim.f.ci %>% group_by(intervention, name) %>% filter(name == "D") %>% summarize(max(mid))
