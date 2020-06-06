library(scales)

# select vars to be plotting 
plot_vars         <- c("cases", "D")

# set the base params that we will then update
source("./COVID_simulate_cont_params.R")

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
SEIR.sim.f.ci.s1    <- SEIR.sim.f.ci %>% mutate(intervention = sim_title)
SEIR.sim.f.t.s1     <- SEIR.sim.f.t %>% mutate(intervention = sim_title)

####
## Simulation: Decrease SIP by 50% on June 8
####
counter.factual    <- FALSE
int.movement       <- "mid"
int.type           <- "none"
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" ## Doesnt matter, ignored for all int.type != "tail"
sim_title  <- "Reduce to 50% of shelter-in-place"
## Save predictions with some reasonable title 
source("./COVID_simulate_cont.R")
SEIR.sim.f.ci.s2    <- SEIR.sim.f.ci %>% mutate(intervention = sim_title)
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
sim_title  <- "Reduce to 50% of shelter-in-place\nwith testing and isolation"
source("./COVID_simulate_cont.R")
## Save predictions with some reasonable title 
SEIR.sim.f.ci.s3    <- SEIR.sim.f.ci %>% mutate(intervention = sim_title)
SEIR.sim.f.t.s3   <- SEIR.sim.f.t %>% mutate(intervention = sim_title)

####
## Simulation: Decrease SIP by 50% + truncate on June 8
####
counter.factual    <- FALSE
int.movement       <- "mid"
int.type           <- "tail"
int.beta_catch     <- 0.05
int.beta_catch_type<- "pct"
int.catch_eff      <-  1
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" 
sim_title          <- "Reduce to 50% of shelter-in-place\nwith superspreading averted"

source("./COVID_simulate_cont.R")
## Save predictions with some reasonable title 
SEIR.sim.f.ci.s4    <- SEIR.sim.f.ci %>% mutate(intervention = sim_title)
SEIR.sim.f.t.s4   <- SEIR.sim.f.t %>% mutate(intervention = sim_title)

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


gg <- ggplot() + 
  
  geom_ribbon(data = (SEIR.sim.f.ci %>% filter(date >= as.Date("2020-06-08")))
  , aes(
    x = date
  , ymin = lwr
  , ymax = upr
  , fill = intervention
    ), alpha = 0.50) +
  
#  geom_line(data = (SEIR.sim.f.c.a %>% filter(date >= as.Date("2020-06-08"), .id != "median"))
#  , aes(
#      x = date
#    , y = cases # mid_c
#    , group  = interaction(.id, intervention)
#    , colour = intervention
#    ), alpha = 0.075, lwd = 0.20) + 
  
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
  
#  geom_line(data = (SEIR.sim.f.c.a %>% filter(date <= as.Date("2020-06-08"), .id != "median", intervention == "Lift all on June 8"))
#  , aes(
#      x = date
#    , y = cases # mid_c
#    , group  = .id
#    ), colour = "grey30", alpha = 0.075, lwd = 0.20) + 
  
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

#  geom_line(data = (SEIR.sim.f.c %>% filter(date <= as.Date("2020-06-08"), intervention == "Lift all on June 8"))
#  , aes(
#      x = date
#    , y = mid_c
#    ), colour = "grey30", alpha = 1, lwd = 1) + 

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
  , legend.position  = c(0.175, 0.875)
  , axis.ticks.x     = element_blank()
  , strip.background = element_blank()
  , strip.placement = "outside"
  , strip.text = element_text(size = 18)) +
  facet_wrap(~ name, nrow = 2, scale = "free_y", 
             strip.position = "left", 
             labeller = as_labeller(c(cases = "Daily Cases", D = "Cumulative Deaths")))
  # annotate("text", x = as.Date("2020-02-01"), y = 4700, label = "a", size = 8, fontface = 2)

######################################################################################################
######################################################################################################
######################################################################################################

####
## Counterfactual: No SIP at all (continue pre-intervention beta through today)
####
counter.factual    <- TRUE
cf.type            <- "no_int"
int.movement       <- "pre"
int.type           <- "none"
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" ## Doesnt matter, ignored for all int.type != "tail"

## Save predictions with some reasonable title 
SEIR.sim.f.c.c1    <- SEIR.sim.f.c %>% mutate(intervention = "No social distancing")
SEIR.sim.f.d.c1    <- SEIR.sim.f.d %>% mutate(intervention = "No social distancing")
SEIR.sim.f.D.c1    <- SEIR.sim.f.D %>% mutate(intervention = "No social distancing")
SEIR.sim.f.c.a.c1  <- SEIR.sim.f.c.a %>% mutate(intervention = "No social distancing")
SEIR.sim.f.d.a.c1  <- SEIR.sim.f.d.a %>% mutate(intervention = "No social distancing")
SEIR.sim.f.D.a.c1  <- SEIR.sim.f.D.a %>% mutate(intervention = "No social distancing")

####
## Counterfactual: Delay SIP by one week
####
counter.factual    <- TRUE
cf.type            <- "delay"
int.movement       <- "pre"
int.type           <- "none"
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" ## Doesnt matter, ignored for all int.type != "tail"

## Save predictions with some reasonable title 
SEIR.sim.f.c.c2    <- SEIR.sim.f.c %>% mutate(intervention = "Social distancing delayed by one week")
SEIR.sim.f.d.c2    <- SEIR.sim.f.d %>% mutate(intervention = "Social distancing delayed by one week")
SEIR.sim.f.D.c2    <- SEIR.sim.f.D %>% mutate(intervention = "Social distancing delayed by one week")
SEIR.sim.f.c.a.c2  <- SEIR.sim.f.c.a %>% mutate(intervention = "Social distancing delayed by one week")
SEIR.sim.f.d.a.c2  <- SEIR.sim.f.d.a %>% mutate(intervention = "Social distancing delayed by one week")
SEIR.sim.f.D.a.c2  <- SEIR.sim.f.D.a %>% mutate(intervention = "Social distancing delayed by one week")

####
## Counterfactual: Lift SIP completely on May 1 (resume pre-intervention beta)
####
counter.factual    <- TRUE
cf.type            <- "may1"
int.movement       <- "pre"
int.type           <- "none"
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" ## Doesnt matter, ignored for all int.type != "tail"

## Save predictions with some reasonable title 
SEIR.sim.f.c.c3    <- SEIR.sim.f.c %>% mutate(intervention = "Social distancing ended on May 1")
SEIR.sim.f.d.c3    <- SEIR.sim.f.d %>% mutate(intervention = "Social distancing ended on May 1")
SEIR.sim.f.D.c3    <- SEIR.sim.f.D %>% mutate(intervention = "Social distancing ended on May 1")
SEIR.sim.f.c.a.c3  <- SEIR.sim.f.c.a %>% mutate(intervention = "Social distancing ended on May 1")
SEIR.sim.f.d.a.c3  <- SEIR.sim.f.d.a %>% mutate(intervention = "Social distancing ended on May 1")
SEIR.sim.f.D.a.c3  <- SEIR.sim.f.D.a %>% mutate(intervention = "Social distancing ended on May 1")

####
## Also need reality. For this just use SEIR.sim.f.c.s1, d, etc.
####

SEIR.sim.f.c.s1.c   <- SEIR.sim.f.c.s1; SEIR.sim.f.c.s1.c$intervention <- "Reality"
SEIR.sim.f.d.s1.c   <- SEIR.sim.f.d.s1; SEIR.sim.f.d.s1.c$intervention <- "Reality"
SEIR.sim.f.D.s1.c   <- SEIR.sim.f.D.s1; SEIR.sim.f.D.s1.c$intervention <- "Reality"
SEIR.sim.f.c.a.s1.c <- SEIR.sim.f.c.a.s1; SEIR.sim.f.c.a.s1.c$intervention <- "Reality"
SEIR.sim.f.d.a.s1.c <- SEIR.sim.f.d.a.s1; SEIR.sim.f.d.a.s1.c$intervention <- "Reality"
SEIR.sim.f.D.a.s1.c <- SEIR.sim.f.D.a.s1; SEIR.sim.f.D.a.s1.c$intervention <- "Reality"

### Combine and plot all of these
SEIR.sim.f.c.c <- rbind(SEIR.sim.f.c.c1, SEIR.sim.f.c.c2, SEIR.sim.f.c.c3, SEIR.sim.f.c.s1.c)
SEIR.sim.f.d.c <- rbind(SEIR.sim.f.d.c1, SEIR.sim.f.d.c2, SEIR.sim.f.d.c3, SEIR.sim.f.d.s1.c)
SEIR.sim.f.D.c <- rbind(SEIR.sim.f.D.c1, SEIR.sim.f.D.c2, SEIR.sim.f.D.c3, SEIR.sim.f.D.s1.c)

SEIR.sim.f.c.a.c <- rbind(SEIR.sim.f.c.a.c1, SEIR.sim.f.c.a.c2, SEIR.sim.f.c.a.c3, SEIR.sim.f.c.a.s1.c)
SEIR.sim.f.d.a.c <- rbind(SEIR.sim.f.d.a.c1, SEIR.sim.f.d.a.c2, SEIR.sim.f.d.a.c3, SEIR.sim.f.d.a.s1.c)
SEIR.sim.f.D.a.c <- rbind(SEIR.sim.f.D.a.c1, SEIR.sim.f.D.a.c2, SEIR.sim.f.D.a.c3, SEIR.sim.f.D.a.s1.c)

SEIR.sim.f.c.a.c.gg <- SEIR.sim.f.c.a.c %>% filter(date <= "2020-06-08")
SEIR.sim.f.c.c.gg   <- SEIR.sim.f.c.c   %>% filter(date <= "2020-06-08")
SEIR.sim.f.D.a.c.gg <- SEIR.sim.f.D.a.c %>% filter(date <= "2020-06-08")
SEIR.sim.f.D.c.gg   <- SEIR.sim.f.D.c   %>% filter(date <= "2020-06-08")

SEIR.sim.f.c.a.c.gg$intervention <- factor(SEIR.sim.f.c.a.c.gg$intervention)
SEIR.sim.f.c.a.c.gg$intervention <- factor(SEIR.sim.f.c.a.c.gg$intervention
  , levels = c("Reality", "Social distancing delayed by one week", "Social distancing ended on May 1", "No social distancing"))
  
SEIR.sim.f.c.c.gg$intervention <- factor(SEIR.sim.f.c.c.gg$intervention)
SEIR.sim.f.c.c.gg$intervention <- factor(SEIR.sim.f.c.c.gg$intervention
  , levels = c("Reality", "Social distancing delayed by one week", "Social distancing ended on May 1", "No social distancing"))

SEIR.sim.f.D.a.c.gg$intervention <- factor(SEIR.sim.f.D.a.c.gg$intervention)
SEIR.sim.f.D.a.c.gg$intervention <- factor(SEIR.sim.f.D.a.c.gg$intervention
  , levels = c("Reality", "Social distancing delayed by one week", "Social distancing ended on May 1", "No social distancing"))
  
SEIR.sim.f.D.c.gg$intervention <- factor(SEIR.sim.f.D.c.gg$intervention)
SEIR.sim.f.D.c.gg$intervention <- factor(SEIR.sim.f.D.c.gg$intervention
  , levels = c("Reality", "Social distancing delayed by one week", "Social distancing ended on May 1", "No social distancing"))

fig4_colors <- c("grey30", "#D67236", "#0b775e", "magenta4")

gg.c <- ggplot(SEIR.sim.f.c.a.c.gg) + 
  
  geom_ribbon(data = (SEIR.sim.f.c.c.gg)
  , aes(
    x = date
  , ymin = lwr_c
  , ymax = upr_c
  , fill = intervention
    ), alpha = 0.50) +
  
  geom_line(data = (SEIR.sim.f.c.a.c.gg %>% filter(.id == "median"))
  , aes(
      x = date
    , y = cases
    , group  = interaction(.id, intervention)
    , colour = intervention
    ), alpha = 1, lwd = 1.5) + 

    geom_point(data = county.data
    , aes(
      x = date
    , y = cases)
      , colour = "black", lwd = 2) + 
  
  scale_colour_manual(
      values = fig4_colors
    , name   = "Counterfactual"
    , labels = c(
      "Reality"
    , "Social distancing delayed
by one week"
    , "Social distancing 
ended on May 1"
    , "No social distancing"
    )) +
  scale_fill_manual(
      values = fig4_colors
    , name   = "Counterfactual"
    , labels = c(
      "Reality"
    , "Social distancing delayed
by one week"
    , "Social distancing 
ended on May 1"
    , "No social distancing"
    )) +
  ylab("Daily Cases") +
  scale_y_continuous(trans = "log10", 
                     breaks = c(1, 10, 100, 1000, 10000)) + 
  scale_x_date(labels = date_format("%b"), date_breaks = "1 month") +
  xlab("") +
  theme(
    axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 12)
  , legend.title     = element_text(size = 12)
  , legend.text      = element_text(margin = margin(t = 0.3, unit = 'cm'), size = 11)
  , plot.title       = element_text(size = 12)
  , legend.position  = c(0.25, 0.80)
  , axis.ticks.x     = element_blank())


gg.D <- ggplot(SEIR.sim.f.D.a.c.gg) + 
  
  geom_ribbon(data = (SEIR.sim.f.D.c.gg)
  , aes(
    x = date
  , ymin = lwr_D
  , ymax = upr_D
  , fill = intervention
    ), alpha = 0.50) +
  
  geom_line(data = (SEIR.sim.f.D.a.c.gg %>% filter(.id == "median"))
  , aes(
      x = date
    , y = D
    , group  = interaction(.id, intervention)
    , colour = intervention
    ), alpha = 1, lwd = 1.5) + 

    geom_point(data = deaths
    , aes(
      x = date
    , y = deaths)
      , colour = "black", lwd = 2) + 
  
  scale_colour_manual(
      values = fig4_colors
    , name   = "Counterfactual"
    , labels = c(
      "Reality"
    , "Social distancing delayed
by one week"
    , "Social distancing 
ended on May 1"
    , "No social distancing"
    )) +
  scale_fill_manual(
      values = fig4_colors
    , name   = "Counterfactual"
    , labels = c(
      "Reality"
    , "Social distancing delayed
by one week"
    , "Social distancing 
ended on May 1"
    , "No social distancing"
    )) +
  ylab("Cumulative Deaths") +
  scale_y_continuous(trans = "log10", 
                     breaks = c(1, 10, 100, 1000, 10000)) + 
  scale_x_date(labels = date_format("%b"), date_breaks = "1 month") +
  xlab("") +
  theme(
    axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 12)
  , legend.title     = element_text(size = 12)
  , legend.text      = element_text(margin = margin(t = 0.3, unit = 'cm'), size = 11)
  , plot.title       = element_text(size = 12)
  , legend.position  = c(0.25, 0.80)
  , axis.ticks.x     = element_blank()) +
   guides(color = FALSE) +
   guides(fill = FALSE) 

gridExtra::grid.arrange(gg.c, gg.D)
