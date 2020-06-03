####
## Simulation: Lift all interventions June 8 (resume pre-intervention beta)
####
counter.factual    <- FALSE
int.movement       <- "pre"
int.type           <- "none"
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" ## Doesnt matter, ignored for all int.type != "tail"

## Save predictions with some reasonable title 
SEIR.sim.f.c.s1    <- SEIR.sim.f.c %>% mutate(intervention = "Lift all on June 8")
SEIR.sim.f.d.s1    <- SEIR.sim.f.d %>% mutate(intervention = "Lift all on June 8")
SEIR.sim.f.D.s1    <- SEIR.sim.f.D %>% mutate(intervention = "Lift all on June 8")
SEIR.sim.f.c.a.s1  <- SEIR.sim.f.c.a %>% mutate(intervention = "Lift all on June 8")
SEIR.sim.f.d.a.s1  <- SEIR.sim.f.d.a %>% mutate(intervention = "Lift all on June 8")
SEIR.sim.f.D.a.s1  <- SEIR.sim.f.D.a %>% mutate(intervention = "Lift all on June 8")

####
## Simulation: Decrease SIP by 50% on June 8
####
counter.factual    <- FALSE
int.movement       <- "mid"
int.type           <- "none"
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" ## Doesnt matter, ignored for all int.type != "tail"

## Save predictions with some reasonable title 
SEIR.sim.f.c.s2    <- SEIR.sim.f.c %>% mutate(intervention = "Reduce to 50% of shelter in place on June 8")
SEIR.sim.f.d.s2    <- SEIR.sim.f.d %>% mutate(intervention = "Reduce to 50% of shelter in place on June 8")
SEIR.sim.f.D.s2    <- SEIR.sim.f.D %>% mutate(intervention = "Reduce to 50% of shelter in place on June 8")
SEIR.sim.f.c.a.s2  <- SEIR.sim.f.c.a %>% mutate(intervention = "Reduce to 50% of shelter in place on June 8")
SEIR.sim.f.d.a.s2  <- SEIR.sim.f.d.a %>% mutate(intervention = "Reduce to 50% of shelter in place on June 8")
SEIR.sim.f.D.a.s2  <- SEIR.sim.f.D.a %>% mutate(intervention = "Reduce to 50% of shelter in place on June 8")

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

## Save predictions with some reasonable title 
SEIR.sim.f.c.s3    <- SEIR.sim.f.c %>% mutate(intervention = "Reduce to 50% of shelter in place with test and isolate June 8")
SEIR.sim.f.d.s3    <- SEIR.sim.f.d %>% mutate(intervention = "Reduce to 50% of shelter in place with test and isolate June 8")
SEIR.sim.f.D.s3    <- SEIR.sim.f.D %>% mutate(intervention = "Reduce to 50% of shelter in place with test and isolate June 8")
SEIR.sim.f.c.a.s3  <- SEIR.sim.f.c.a %>% mutate(intervention = "Reduce to 50% of shelter in place with test and isolate June 8")
SEIR.sim.f.d.a.s3  <- SEIR.sim.f.d.a %>% mutate(intervention = "Reduce to 50% of shelter in place with test and isolate June 8")
SEIR.sim.f.D.a.s3  <- SEIR.sim.f.D.a %>% mutate(intervention = "Reduce to 50% of shelter in place with test and isolate June 8")

### Combine and plot all of these
SEIR.sim.f.c <- rbind(SEIR.sim.f.c.s1, SEIR.sim.f.c.s2, SEIR.sim.f.c.s3)
SEIR.sim.f.d <- rbind(SEIR.sim.f.d.s1, SEIR.sim.f.d.s2, SEIR.sim.f.d.s3)
SEIR.sim.f.D <- rbind(SEIR.sim.f.D.s1, SEIR.sim.f.D.s2, SEIR.sim.f.D.s3)

SEIR.sim.f.c.a <- rbind(SEIR.sim.f.c.a.s1, SEIR.sim.f.c.a.s2, SEIR.sim.f.c.a.s3)
SEIR.sim.f.d.a <- rbind(SEIR.sim.f.d.a.s1, SEIR.sim.f.d.a.s2, SEIR.sim.f.d.a.s3)
SEIR.sim.f.D.a <- rbind(SEIR.sim.f.D.a.s1, SEIR.sim.f.D.a.s2, SEIR.sim.f.D.a.s3)

fig2_colors <- c("#B40F20", "#E58601", "#046C9A", "grey30") # old colors

library(scales)

# SEIR.sim.f.c.a <- SEIR.sim.f.c.a %>% filter(date > as.Date(sir_init.mid.t) + 1)
# SEIR.sim.f.D.a <- SEIR.sim.f.D.a %>% filter(date > as.Date(sir_init.mid.t) + 1)

gg.c <- ggplot(SEIR.sim.f.c.a) + 
  
  geom_ribbon(data = (SEIR.sim.f.c %>% filter(date >= as.Date("2020-06-08")))
  , aes(
    x = date
  , ymin = lwr_c
  , ymax = upr_c
  , fill = intervention
    ), alpha = 0.50) +
  
#  geom_line(data = (SEIR.sim.f.c.a %>% filter(date >= as.Date("2020-06-08"), .id != "median"))
#  , aes(
#      x = date
#    , y = cases # mid_c
#    , group  = interaction(.id, intervention)
#    , colour = intervention
#    ), alpha = 0.075, lwd = 0.20) + 
  
  geom_line(data = (SEIR.sim.f.c.a %>% filter(date >= as.Date("2020-06-08"), .id == "median"))
  , aes(
      x = date
    , y = cases # mid_c
    , group  = interaction(.id, intervention)
    , colour = intervention
    ), alpha = 1, lwd = 1.5) + 

  geom_ribbon(data = (SEIR.sim.f.c %>% filter(date <= as.Date("2020-06-08")))
  , aes(
    x = date
  , ymin = lwr_c
  , ymax = upr_c
    ), fill = "grey30", alpha = 0.50) +
  
#  geom_line(data = (SEIR.sim.f.c.a %>% filter(date <= as.Date("2020-06-08"), .id != "median", intervention == "Lift all on June 8"))
#  , aes(
#      x = date
#    , y = cases # mid_c
#    , group  = .id
#    ), colour = "grey30", alpha = 0.075, lwd = 0.20) + 
  
  geom_line(data = (SEIR.sim.f.c.a %>% filter(date <= as.Date("2020-06-08"), .id == "median", intervention == "Lift all on June 8"))
  , aes(
      x = date
    , y = cases # mid_c
    , group  = .id
    ), colour = "grey30", alpha = 1, lwd = 1.5) + 

#  geom_line(data = (SEIR.sim.f.c %>% filter(date <= as.Date("2020-06-08"), intervention == "Lift all on June 8"))
#  , aes(
#      x = date
#    , y = mid_c
#    ), colour = "grey30", alpha = 1, lwd = 1) + 

    geom_point(data = county.data
    , aes(
      x = date
    , y = cases)
      , colour = "black", lwd = 2) + 
  
  geom_vline(xintercept = as.Date("2020-06-08"), linetype = "dashed", lwd = 0.65) +
  
  scale_colour_manual(
      values = fig2_colors
    , name   = "Intervention"
    , labels = c(
      "Lift all on June 8"
    , "Reduce to 50% of shelter
in place on June 8"
    , "Reduce to 50% of shelter 
in place with 
test and isolate on June 8"
    )) +
  scale_fill_manual(
      values = fig2_colors
    , name   = "Intervention"
    , labels = c(
      "Lift all on June 8"
    , "Reduce to 50% of shelter
in place on June 8"
    , "Reduce to 50% of shelter 
in place with 
test and isolate on June 8"
    )) +
  ylab("Daily Cases") +
  scale_y_continuous(trans = "log10", 
                     breaks = c(1, 10, 100, 1000, 10000)) + 
  scale_x_date(labels = date_format("%b"), date_breaks = "2 month") +
  xlab("") +
  theme(
    axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 12)
  , legend.title     = element_text(size = 12)
  , legend.text      = element_text(margin = margin(t = 0.3, unit = 'cm'), size = 11)
  , plot.title       = element_text(size = 12)
  , legend.position  = c(0.25, 0.80)
  , axis.ticks.x     = element_blank())
  # annotate("text", x = as.Date("2020-02-01"), y = 4700, label = "a", size = 8, fontface = 2)

gg.D <- ggplot(SEIR.sim.f.D.a) + 
  
  geom_ribbon(data = (SEIR.sim.f.D %>% filter(date >= as.Date("2020-06-08")))
  , aes(
    x = date
  , ymin = lwr_D
  , ymax = upr_D
  , fill = intervention
    ), alpha = 0.50) +
  
## Lines instead of 
#  geom_line(data = (SEIR.sim.f.D.a %>% filter(date >= as.Date("2020-06-08"), .id != "median"))
#  , aes(
#      x = date
#    , y = D # mid_c
#    , group  = interaction(.id, intervention)
#    , colour = intervention
#    ), alpha = 0.25, lwd = 0.25) + 
  
  geom_line(data = (SEIR.sim.f.D.a %>% filter(date >= as.Date("2020-06-08"), .id == "median"))
  , aes(
      x = date
    , y = D # mid_c
    , group  = interaction(.id, intervention)
    , colour = intervention
    ), alpha = 1, lwd = 1.5) + 

  geom_ribbon(data = (SEIR.sim.f.D %>% filter(date <= as.Date("2020-06-08")))
  , aes(
    x = date
  , ymin = lwr_D
  , ymax = upr_D
    ), fill = "grey30", alpha = 0.50) +
  
#  geom_line(data = (SEIR.sim.f.D.a %>% filter(date <= as.Date("2020-06-08"), .id != "median", intervention == "Lift all on June 8"))
#  , aes(
#      x = date
#    , y = D # mid_c
#    , group  = .id
#    ), colour = "grey30", alpha = 0.25, lwd = 0.25) + 
  
  geom_line(data = (SEIR.sim.f.D.a %>% filter(date <= as.Date("2020-06-08"), .id == "median", intervention == "Lift all on June 8"))
  , aes(
      x = date
    , y = D # mid_c
    , group  = .id
    ), colour = "grey30", alpha = 1, lwd = 1.5) + 

#  geom_line(data = (SEIR.sim.f.c %>% filter(date <= as.Date("2020-06-08"), intervention == "Lift all on June 8"))
#  , aes(
#      x = date
#    , y = mid_c
#    ), colour = "grey30", alpha = 1, lwd = 1) + 

    geom_point(data = deaths
    , aes(
      x = date
    , y = deaths)
      , colour = "black", lwd = 2) + 
  
  geom_vline(xintercept = as.Date("2020-06-08"), linetype = "dashed", lwd = 0.65) +
  
  scale_colour_manual(
      values = fig2_colors
    , name   = "Intervention"
    , labels = c(
      "Lift all on June 8"
    , "Reduce to 50% of shelter
in place on June 8"
    , "Reduce to 50% of shelter 
in place with 
test and isolate June 8"
    )) +
  scale_fill_manual(
      values = fig2_colors
    , name   = "Intervention"
    , labels = c(
      "Lift all on June 8"
    , "Reduce to 50% of shelter
in place on June 8"
    , "Reduce to 50% of shelter 
in place with 
test and isolate June 8"
    )) +
  guides(fill = FALSE) +
  ylab("Cumulative Deaths") +
  scale_y_continuous(trans = "log10", 
                     breaks = c(1, 10, 100, 1000, 10000)) + 
  scale_x_date(labels = date_format("%b"), date_breaks = "2 month") +
  xlab("Date") +
  theme(
    axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 12)
  , legend.title     = element_text(size = 14)
  , legend.text      = element_text(margin = margin(t = 0.5, unit = 'cm'), size = 11)
  , plot.title       = element_text(size = 12)
  , legend.position  = c(0.25, 0.80)
  , axis.ticks.x     = element_blank()) +
   guides(color = FALSE) 
  # annotate("text", x = as.Date("2020-02-01"), y = 4700, label = "a", size = 8, fontface = 2)

gridExtra::grid.arrange(gg.c, gg.D)

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
