#########################################
## Summarize fits from scc_forecasts.R ##
#########################################

variable_params.l <- melt(variable_params[, c("E0", "sim_start", "int_length2", "sd_m2", "sd_m1", "int_length1", "paramset")], "paramset")
SEIR.out          <- left_join(variable_params.l, SEIR.sim.ss.t.ci, "paramset")
names(SEIR.out)   <- c("param_num", "Param.Name", "Param.Value", "Out.Name", "lwr", "est", "upr")

SEIR.sim.ss.t.ci.gg <- left_join(
  SEIR.sim.ss.t.ci
, variable_params
, by = "paramset")

SEIR.sim.ss.t.ci.gg.r <- SEIR.sim.ss.t.ci.gg %>%
  dplyr::filter(name != "total_R")

SEIR.sim.ss.t.ci.gg.r[SEIR.sim.ss.t.ci.gg.r$name == "max_H", ]$name       <- "Maximum Simultaneous Hospitalizations"
SEIR.sim.ss.t.ci.gg.r[SEIR.sim.ss.t.ci.gg.r$name == "prop_ext", ]$name    <- "Proportion of Stochastic Runs that End"
SEIR.sim.ss.t.ci.gg.r[SEIR.sim.ss.t.ci.gg.r$name == "total_D", ]$name     <- "Total Deaths"
SEIR.sim.ss.t.ci.gg.r[SEIR.sim.ss.t.ci.gg.r$name == "when_max_H", ]$name  <- "Week of Maximum Hospitalizations"
SEIR.sim.ss.t.ci.gg.r[SEIR.sim.ss.t.ci.gg.r$name == "when_red_H", ]$name  <- "Week of First Decrease in Hospitalizations"

SEIR.sim.ss.t.ci.gg.r   <- SEIR.sim.ss.t.ci.gg.r %>% 
  dplyr::filter(SEIR.sim.ss.t.ci.gg.r$name != "Proportion of Stochastic Runs that End")

#SEIR.sim.ss.t.ci.gg.r.r <- rbind(
#  SEIR.sim.ss.t.ci.gg.r.0
#, SEIR.sim.ss.t.ci.gg.r.60
#, SEIR.sim.ss.t.ci.gg.r.90
#  )

## stupid dplyr error with custom function... :(
SEIR.sim.ss.t.ci.gg.r <- SEIR.sim.ss.t.ci.gg.r %>%
  mutate(Baseline = 0, Intervention = 0)
for (i in 1:nrow(SEIR.sim.ss.t.ci.gg.r)) {
SEIR.sim.ss.t.ci.gg.r[i, ]$Baseline     <- with(SEIR.sim.ss.t.ci.gg.r[i, ], covid_R0(beta0est, fixed_params, 1))
SEIR.sim.ss.t.ci.gg.r[i, ]$Intervention <- with(SEIR.sim.ss.t.ci.gg.r[i, ], covid_R0(beta0est, fixed_params, sd_m2))
}

R0_hist     <- SEIR.sim.ss.t.ci.gg.r[!duplicated(SEIR.sim.ss.t.ci.gg.r[, c("beta0est", "Baseline", "Intervention", "beta0est")]), ]
R0_hist     <- melt(R0_hist[, c("Baseline", "Intervention")])
