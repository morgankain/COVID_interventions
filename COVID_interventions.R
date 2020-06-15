# set the base params that we will then update
source("./COVID_simulate_cont_params.R")

focal.county       <- "Contra Costa"
focal.state_abbr   <- "CA"
rds.name           <- "output/Contra Costa_TRUE_FALSE_0_2020-06-10_cont_temp_NZ_sig_exp_gamma_full_try_NA.Rds"

####
## Simulation: Lift all interventions June 8 (resume pre-intervention beta)
####
counter.factual    <- FALSE
int.movement       <- "pre"
int.type           <- "none"
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" ## Doesnt matter, ignored for all int.type != "tail"
sim_title          <- "Lift all"
# source sim script here
source("./COVID_simulate_cont.R")
## Save CIS and trajectories of predictions with some reasonable title 
SEIR.sim.f.ci.s1   <- SEIR.sim.f.ci %>% mutate(intervention = sim_title)
SEIR.sim.f.t.s1    <- SEIR.sim.f.t %>% mutate(intervention = sim_title)

####
## Simulation: Decrease SIP by 50% on June 8
####
counter.factual    <- FALSE
int.movement       <- "mid"
int.type           <- "none"
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" ## Doesnt matter, ignored for all int.type != "tail"
sim_title          <- "Reduce to 50% of shelter-in-place"
## Save predictions with some reasonable title 
source("./COVID_simulate_cont.R")
SEIR.sim.f.ci.s2   <- SEIR.sim.f.ci %>% mutate(intervention = sim_title)
SEIR.sim.f.t.s2    <- SEIR.sim.f.t %>% mutate(intervention = sim_title)

####
## Simulation: Decrease SIP by 50% + test & isolate on June 8
####
counter.factual    <- FALSE
int.movement       <- "mid"
int.type           <- "inf_iso"
iso_mild_level     <- 0.2
iso_severe_level   <- 0.2
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" ## Doesnt matter, ignored for all int.type != "tail"
sim_title          <- "Reduce to 50% of shelter-in-place\nwith testing and isolation"
source("./COVID_simulate_cont.R")
## Save predictions with some reasonable title 
SEIR.sim.f.ci.s3   <- SEIR.sim.f.ci %>% mutate(intervention = sim_title)
SEIR.sim.f.t.s3    <- SEIR.sim.f.t %>% mutate(intervention = sim_title)

####
## Simulation: Decrease SIP by 50% + truncate on June 8
####
counter.factual    <- FALSE
int.movement       <- "mid"
int.type           <- "tail"
int.beta_catch     <- 0.05            ## beta0 values caught by intervention; alternatively, specify by top percent of distribution to trim
int.beta_catch_type<- "pct"           ## if pct, treated as percentile, otherwise, as absoulte value
int.catch_eff      <-  1              ## effectiveness at catching beta0 values above the beta_catch (0 - 1)
int.init           <- "2020-06-08"
int.end            <- "2020-10-01" 
sim_title          <- "Reduce to 50% of shelter-in-place\nwith superspreading averted"

source("./COVID_simulate_cont.R")
## Save predictions with some reasonable title 
SEIR.sim.f.ci.s4   <- SEIR.sim.f.ci %>% mutate(intervention = sim_title)
SEIR.sim.f.t.s4    <- SEIR.sim.f.t %>% mutate(intervention = sim_title)

### Combine and plot all of these
SEIR.sim.f.t <- rbind(SEIR.sim.f.t.s1, 
                      SEIR.sim.f.t.s2, 
                      SEIR.sim.f.t.s3, 
                      SEIR.sim.f.t.s4)
SEIR.sim.f.ci <- rbind(SEIR.sim.f.ci.s1, 
                      SEIR.sim.f.ci.s2, 
                      SEIR.sim.f.ci.s3, 
                      SEIR.sim.f.ci.s4)

fig2_colors <- c("#B40F20", "#E58601", "#046C9A", "darkgreen", "grey30") # old colors

{
(
gg <- ggplot() + 
  
  geom_ribbon(data = (SEIR.sim.f.ci %>% filter(date >= as.Date("2020-06-08")))
  , aes(
    x = date
  , ymin = lwr
  , ymax = upr
  , fill = intervention
    ), alpha = 0.50) +
  
  geom_line(data = (SEIR.sim.f.t %>% filter(date >= as.Date("2020-06-08"), .id == "median"))
  , aes(
      x = date
    , y = value # mid_c
    , group  = interaction(.id, intervention)
    , colour = intervention
    ), alpha = 1, lwd = 1.5) + 

  geom_ribbon(data = (SEIR.sim.f.ci %>% filter(date <= as.Date("2020-06-08")))
  , aes(
    x = date
  , ymin = lwr
  , ymax = upr
    ), fill = "grey30", alpha = 0.50) +
  
  geom_line(data = (SEIR.sim.f.t %>% 
                      {filter(., date <= as.Date("2020-06-08"), 
                             .id == "median", 
                             intervention == pull(., intervention) %>% 
                               first)})
  , aes(
      x = date
    , y = value # mid_c
    , group  = .id
    ), colour = "grey30", alpha = 1, lwd = 1.5) + 

    geom_point(data = county.data %>% 
                 mutate(D = cumsum(ifelse(is.na(deaths), 0, deaths))*
                          ifelse(is.na(deaths), NA, 1)) %>% 
                 select(date, any_of(plot_vars)) %>% 
                 pivot_longer(any_of(plot_vars))
    , aes(
      x = date
    , y = value)
      , colour = "black", lwd = 2) + 
  
  geom_vline(xintercept = as.Date("2020-06-08"), linetype = "dashed", lwd = 0.65) +
  
  scale_colour_manual(
      values = fig2_colors
    , name   = "Scenario"
    # , labels = sim_names
    ) +
  scale_fill_manual(
      values = fig2_colors
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
  , legend.position  = c(0.245, 0.875)
  , axis.ticks.x     = element_blank()
  , strip.background = element_blank()
  , strip.placement = "outside"
  , strip.text = element_text(size = 18)) +
  facet_wrap(~ name, nrow = 2, scale = "free_y", 
             strip.position = "left", 
             labeller = as_labeller(c(cases = "Daily Cases", D = "Cumulative Deaths")))
)
}
