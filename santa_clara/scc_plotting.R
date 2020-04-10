###########################################
## Plot output of sobol parameter search ##
###########################################

## Hideous hupercube plot
ggres1 <- ggplot(SEIR.out
  , aes(Param.Value, est)) +
    geom_point(lwd = 1) +
    facet_grid(Out.Name ~ Param.Name, scale = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(trans = "pseudo_log") +
    xlab("Parameter Values") +
    ylab("COVID Outcome Summary")

## Plot with only runs where strength of social distancing is lowered
gg1_only_red <- ggplot(SEIR.sim.ss.t.ci.gg.r[SEIR.sim.ss.t.ci.gg.r$name != "Week of First Decrease in Hospitalizations", ]
  , aes(int_length2, est)) +
  geom_point(aes(colour = sd_m2), lwd = 2) +
#  geom_errorbar(aes(x = int_length2, ymin = lwr, ymax = upr), width = 0.1, alpha = 0.2) +
  scale_color_gradient(low = "steelblue3", high = "firebrick3", name = "Proportion of 
Contacts Remaining") +
  facet_wrap(~name, scale = "free") +
  xlab("Length of Social Distancing Intervention") +
  ylab("COVID Summary Estimate") +
  theme(strip.text.x = element_text(size = 12, colour = "black", face = "bold"))

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

gghist <- ggplot(R0_hist, aes(x = value)) + 
 geom_histogram(aes( y = ..density.., fill = variable), colour = "black", bins = 50, lwd = 0.2, alpha = 0.5) +
 geom_density(aes(fill = variable), alpha = .2) +
 scale_fill_manual(name = "", values = c("firebrick3", "steelblue3")) +
 xlab("R0 Estimate") + ylab("Density") 

ggfilename <- paste("output/", paste("R0_hist", focal.county, max(deaths$date), sep = "_"), sep = "")

ggsave(filename = paste(ggfilename, ".pdf", sep = ""), gghist, device = "pdf")
