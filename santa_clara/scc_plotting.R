###########################################
## Plot output of sobol parameter search ##
###########################################

variable_params.l <- melt(variable_params[, c("E0", "sim_start", "int_length2", "sd_m2", "sd_m1", "int_length1", "paramset")], "paramset")
SEIR.out          <- left_join(variable_params.l, SEIR.sim.ss.t.ci, "paramset")
names(SEIR.out)   <- c("param_num", "Param.Name", "Param.Value", "Out.Name", "lwr", "est", "upr")

## Hideous hupercube plot
ggres1 <- ggplot(SEIR.out
  , aes(Param.Value, est)) +
    geom_point(lwd = 1) +
    facet_grid(Out.Name ~ Param.Name, scale = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(trans = "pseudo_log") +
    xlab("Parameter Values") +
    ylab("COVID Outcome Summary")

SEIR.sim.ss.t.ci.gg <- left_join(
  SEIR.sim.ss.t.ci
, variable_params
, by = "paramset")

SEIR.sim.ss.t.ci.gg.r <- SEIR.sim.ss.t.ci.gg %>%
  filter(name != "total_R")

SEIR.sim.ss.t.ci.gg.r[SEIR.sim.ss.t.ci.gg.r$name == "max_H", ]$name       <- "Maximum Simultaneous Hospitalizations"
SEIR.sim.ss.t.ci.gg.r[SEIR.sim.ss.t.ci.gg.r$name == "prop_ext", ]$name    <- "Proportion of Stochastic Runs that End"
SEIR.sim.ss.t.ci.gg.r[SEIR.sim.ss.t.ci.gg.r$name == "total_D", ]$name     <- "Total Deaths"
SEIR.sim.ss.t.ci.gg.r[SEIR.sim.ss.t.ci.gg.r$name == "when_max_H", ]$name  <- "Week of Maximum Hospitalizations"
SEIR.sim.ss.t.ci.gg.r[SEIR.sim.ss.t.ci.gg.r$name == "when_red_H", ]$name  <- "Week of First Decrease in Hospitalizations"

SEIR.sim.ss.t.ci.gg.r   <- SEIR.sim.ss.t.ci.gg.r %>% filter(SEIR.sim.ss.t.ci.gg.r$name != "Proportion of Stochastic Runs that End")

SEIR.sim.ss.t.ci.gg.r <- SEIR.sim.ss.t.ci.gg.r %>% mutate(iso_inf = 90)

SEIR.sim.ss.t.ci.gg.r.r <- rbind(
  SEIR.sim.ss.t.ci.gg.r.F
, SEIR.sim.ss.t.ci.gg.r.T
, SEIR.sim.ss.t.ci.gg.r.90
  )

SEIR.sim.ss.t.ci.gg.r.r <- SEIR.sim.ss.t.ci.gg.r.r %>%
  filter(paramset < 299) 

gg1 <- ggplot(SEIR.sim.ss.t.ci.gg.r.r[SEIR.sim.ss.t.ci.gg.r.r$name != "Week of First Decrease in Hospitalizations", ]
  , aes(int_length2, est)) +
  geom_point(aes(colour = sd_m2, shape = as.factor(iso_inf)), lwd = 2) +
#  geom_errorbar(aes(x = int_length2, ymin = lwr, ymax = upr), width = 0.1, alpha = 0.2) +
  scale_color_gradient(low = "steelblue3", high = "firebrick3", name = "Proportion of 
Contacts Remaining") +
  scale_shape_manual(values = c(19, 2, 4), name = "Infected Isolation", labels = c("No", "Yes")) +
  # scale_y_log10() + 
  facet_wrap(~name, scale = "free") +
  xlab("Length of Social Distancing Intervention") +
  ylab("COVID Summary Estimate") +
  theme(strip.text.x = element_text(size = 12, colour = "black", face = "bold"))

SEIR.sim.ss.t.ci.gg.r.r.d <- SEIR.sim.ss.t.ci.gg.r.r %>% 
  filter(name == "Week of Maximum Hospitalizations" | name == "Week of First Decrease in Hospitalizations") %>% 
  mutate(est.d = as.Date(est * 7, origin = sim_start))

gg2 <- ggplot(SEIR.sim.ss.t.ci.gg.r.r.d
  , aes(int_length2, est.d)) +
  geom_point(aes(colour = sd_m2, shape = iso_inf), lwd = 2.5) +
#  geom_errorbar(aes(x = int_length2, ymin = lwr, ymax = upr), width = 0.1, alpha = 0.2) +
  scale_color_gradient(low = "steelblue3", high = "firebrick3", name = "Proportion of 
Contacts Remaining") +
  scale_shape_manual(values = c(19, 2)) +
  facet_wrap(~name, scale = "free") +
  xlab("Length of Social Distancing Intervention") +
  ylab("COVID Summary Estimate") 

## R0
# beta0est * (1/3) * (2/3) * (7) + beta0est * (2/3) * (2 + 5.76) * (1 - 0.956) + beta0est * (2/3) * (2 + 7) * 0.956


SEIR.sim.ss.t.ci.gg.r.r <- SEIR.sim.ss.t.ci.gg.r.r %>%
  mutate(
    Baseline      = beta0est * (1/3) * (2/3) * (7) + beta0est * (2/3) * (0.5 + 5.76) * (1 - 0.956) + beta0est * (2/3) * (0.5 + 7) * 0.956
      ) %>%
  mutate(
    Intervention  = sd_m2 * beta0est * (1/3) * (2/3) * (7) + sd_m2 * beta0est * (2/3) * (0.5 + 5.76) * (1 - 0.956) + sd_m2 * beta0est * (2/3) * (0.5 + 7) * 0.956
  )

R0_hist     <- SEIR.sim.ss.t.ci.gg.r.r[!duplicated(SEIR.sim.ss.t.ci.gg.r.r[, c("beta0est", "Baseline", "Intervention", "beta0est")]), ]
R0_hist     <- melt(R0_hist[, c("Baseline", "Intervention", )])

ggplot(R0_hist, aes(x = value)) + 
 geom_histogram(aes( y = ..density.., fill = variable), colour = "black", bins = 50, lwd = 0.2, alpha = 0.5) +
 geom_density(aes(fill = variable), alpha = .2) +
 scale_fill_manual(name = "", values = c("firebrick3", "steelblue3")) +
 xlab("R0 Estimate") + ylab("Density") +
 geom_vline(xintercept =
     0.285 * 0.5422 * (1/3) * (2/3) * (7) + 
     0.285 * 0.5422 * (2/3) * (0.5 + 5.76) * (1 - 0.956) + 
     0.285 * 0.5422 * (2/3) * (0.5 + 7) * 0.956
   , lwd = 1, linetype = "dashed") +
  annotate("text", x = 2.75, y = 1.2, label = "72% reduction in contacts at 
average estimated base transmission rate
for R0 = 1")


