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

ggplot(SEIR.sim.ss.t.ci.gg.r
  , aes(int_length2, est)) +
  geom_point(aes(colour = sd_m2)) +
  geom_errorbar(aes(x = int_length2, ymin = lwr, ymax = upr), width = 0.1) +
  facet_wrap(~name, scale = "free") +
  xlab("Length of Social Distancing Intervention") +
  ylab("COVID Summary Estimate") 

  