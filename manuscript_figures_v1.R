# running this script creates figures but requires pre-exisisting fits and simulation to run

# set up ---
needed_packages <- c(
  "pomp"
  , "plyr"
  , "dplyr"
  , "ggplot2"
  , "magrittr"
  , "scales"
  , "lubridate"
  , "tidyr")

lapply(needed_packages, require, character.only = TRUE)

# colors 
fig2_colors = c("#B40F20", "#E58601", "#046C9A", "grey30") # old colors
fig3_colors = c("#273046",
                "#366d8a",
                "#78B7C5")
fig4_colors = c("grey30",
                "#D67236",
                "#0b775e")

source("ggplot_theme.R")
fitting        <- F
fit.E0         <- F
source("COVID_pomp.R")

# read deaths data used for ditting 
deaths <- read.csv("us-counties.txt") %>% 
  mutate(date = as.Date(date)) %>% 
  filter(county == "Santa Clara") %>% 
  mutate(deaths_cum = deaths) %>% 
  mutate(deaths = deaths_cum - lag(deaths_cum)) %>% 
  replace_na(list(deaths = 0)) 


####
## Figure 1 ----
####
# load files
all_fits = alply(list.files("output", full.names = T, pattern = "Santa Clara"),
                 1, 
                 function(file_path){
                   print(file_path)
                   fit <- readRDS(file_path)
                   full_Reff <- matrix(nrow = dim(fit$param_array)[1], 
                                       ncol = dim(fit$param_array)[2], 
                                       data = 0)
                   full_R0 <- matrix(nrow = dim(fit$param_array)[1], 
                                     ncol = dim(fit$param_array)[2], 
                                     data = 0)
                   # still using double forloop...                     
                   for (i in 1:dim(full_Reff)[1]) {
                     temp_dyn <- fit$dynamics_summary %>% filter(paramset == i)
                     for (j in 1:dim(full_Reff)[2]) {
                       full_Reff[i, j] <- with(fit$variable_params[i, ], covid_R0(
                         beta0est     = fit$param_array[i, j, 1]
                         , fixed_params = c(fit$fixed_params, unlist(fit$variable_params[i, ]))
                         , sd_strength  = fit$param_array[i, j, 2]
                         , prop_S = unlist(temp_dyn[temp_dyn$name == "S_now", 3]) /
                           (fit$fixed_params["N"] - temp_dyn[temp_dyn$name == "total_D", 3]))
                       )
                       full_R0[i, j] <- with(fit$variable_params[i, ], covid_R0(
                         beta0est     = fit$param_array[i, j, 1]
                         , fixed_params = c(fit$fixed_params, unlist(fit$variable_params[i, ]))
                         , sd_strength  = 1
                         , prop_S = 1)
                       )
                     }
                   }
                   fit$variable_params$county = gsub("./output/", "", 
                                                     strsplit(file_path, split = "_")[[1]][1])
                   fit$fit_sip_eff = strsplit(file_path, split = "_")[[1]][2]
                   fit$days_minus = as.numeric(strsplit(file_path, split = "_")[[1]][4])
                   fit$county = gsub("./output/", "", 
                                     strsplit(file_path, split = "_")[[1]][1])
                   fit$variable_params$fit_sip_eff = strsplit(file_path, split = "_")[[1]][2]
                   fit$variable_params$days_minus = as.numeric(strsplit(file_path, split = "_")[[1]][4])
                   fit$variable_params %<>% mutate(date = as.Date("2020-04-22") - days_minus)
                   fit$full_Reff = full_Reff
                   fit$full_R0 = full_R0
                   return(fit)
                 })

names(all_fits) = gsub(".rds", "",
                       sapply(strsplit(list.files("./output", full.names = T, pattern = "Santa Clara"), 
                                       split = "/"), "[", 3))

# pull out all Reff and R0 estimates from the most recent fit in Santa Clara County
R_hist <- rbind(data.frame(name = "Reff", 
                            value = as.vector(all_fits[["Santa Clara_TRUE_FALSE_0_2020-04-25_final"]]$full_Reff),
                            county = "Santa Clara"),
                 data.frame(name = "R0", 
                            value = as.vector(all_fits[["Santa Clara_TRUE_FALSE_0_2020-04-25_final"]]$full_R0),
                            county = "Santa Clara")) 

# make plot with Reff and R0 histogram
gg1 <- R_hist %>% mutate(name = factor(name, levels = c("R0", "Reff"))) %>%
  ggplot(aes(x = value, fill = name)) + 
  geom_histogram(aes(y = ..density..,),
                 position = "identity",
                 colour = "black", bins = 50,
                 lwd = 0.2, alpha = 0.5) +
  geom_density(alpha = 0.7) +
  geom_vline(xintercept = 1, lty = "dashed", alpha = 0.8) + 
  scale_fill_manual(name = "", 
                    values = fig2_colors[1:2]) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.text=element_text(size=15),
        axis.title=element_text(size=18)) +
  annotate("text", x = 1.5, y = 3.3, 
           hjust = 0.5,
           size = 6,
           parse = TRUE, 
           label = expression(atop(paste("R" [e]*" under"), "shelter-in-place"))
  ) +
  annotate("text", x = 3, y = 2.1, 
           hjust = 0.5,
           size = 6,
           parse = TRUE, label = as.character(expression("R"[0]))) +
  xlab(expression("R"[0]^{}*" estimate")) + 
  ylab("Density")

# save the figure
ggsave("manuscript_materials/Figure1.pdf", plot = gg1, 
       device = "pdf", width = 8, height = 6,
       units = "in")
####
## Figure 2 ----
####
format.labels = function(label, char_max = 6){stringr::str_pad(label, max(char_max), pad = " ")}

fig2_traj <- readRDS("manuscript_materials/fig2data.Rds")
SEIR.sim.f.a    <- gg.data[[1]]
SEIR.sim.f.cf.a <- gg.data[[2]]
gg.data <- readRDS("output/counterfactual_hist.Rds")
SEIR.sim.ss.m   <- gg.data[[1]]
SEIR.sim.ss.a   <- gg.data[[2]]



# 2A: deaths ----

gg2.a <- ggplot(SEIR.sim.f.a) + 
  geom_line(data = (SEIR.sim.f.a %>% filter(.id != "median") %>% filter(date >= as.Date("2020-06-01")))
  , aes(x     = date
      , y     = D
      , group = interaction(.id, Intervention)
      , colour = Intervention
    ), alpha = 0.1) + 
  geom_line(data = (SEIR.sim.f.a %>% filter(.id == "median")%>% filter(date >= as.Date("2020-06-01")))
  , aes(x      = date
      , y      = D
      , colour = Intervention
    ), lwd = 1.5) +
  scale_colour_manual(values = fig2_colors) +
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12)
  , legend.title = element_text(size = 12)
  , plot.title = element_text(size = 12)) +
  ylab("Cumulative Deaths") +
  # scale_y_log10(breaks = c(1, 10, 100, 1000, 5000), labels = format.labels(c(1, 10, 100, 1000, 5000))) +
  scale_y_continuous(trans = "log10", 
                     breaks = c(1, 10, 100, 1000, 5000),
                     labels = format.labels) + 
  geom_line(data = (SEIR.sim.f.a %>% filter(.id != "median") %>% filter(date <= as.Date("2020-06-01")))
  , aes(x     = date
      , y     = D
      , group = .id
    ), alpha = 0.1, colour = "grey30") +
  geom_line(data = (SEIR.sim.f.a %>% filter(.id == "median") %>% filter(date <= as.Date("2020-06-01")))
            , aes(x     = date
                  , y     = D
                  , group = .id
            ), lwd = 1.5, colour = "grey30") +
  # geom_point(data = deaths, aes(date, deaths_cum), lwd = 2, colour = "forestgreen") +
  geom_point(data = deaths, aes(date, deaths_cum), lwd = 2, colour = "black") +
  geom_vline(xintercept = as.Date("2020-06-01"), linetype = "dashed", lwd = 0.65) +
  theme(legend.title = element_text(size = 14),
        legend.position = c(0.8, 0.25),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
        ) +
  scale_x_date(labels = date_format("%b"), date_breaks = "2 month") +
  xlab("") +
  annotate("text", x = as.Date("2020-02-01"), y = 4700, label = "a",
           size = 8, fontface = 2)
# gg2.a 

## 2B: cases ----

gg2.b <- ggplot(SEIR.sim.f.a) + 
  geom_line(data = (SEIR.sim.f.a %>% filter(.id != "median")  %>% filter(date >= as.Date("2020-06-01")))
  , aes(x     = date
      , y     = I 
      , group = interaction(.id, Intervention)
      , colour = Intervention
    ), alpha = 0.1) + 
  geom_line(data = (SEIR.sim.f.a %>% filter(.id == "median") %>% filter(date >= as.Date("2020-06-01")))
  , aes(x      = date
      , y      = I 
      , colour = Intervention
    ), lwd = 1.5) +
  scale_x_date(labels = date_format("%b"), date_breaks = "2 month") +
  scale_colour_manual(values = fig2_colors) + 
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(angle = 0.5, hjust = 0.5, size = 12)
  , legend.title = element_text(size = 12)
  , plot.title = element_text(size = 12)) +
  xlab("Date") + 
  ylab("Concurrent Infections") +
  # scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000), labels = format.labels(c(1, 10, 100, 1000, 10000, 100000))) +
  scale_y_continuous(trans = "log10", 
                     breaks = c(1, 10, 100, 1000, 10000, 100000),
                     labels = c(1, 10, 100, 1000, 10000, 100000)) + 
  geom_line(data = (SEIR.sim.f.a %>% filter(.id != "median") %>% filter(date <= as.Date("2020-06-01")))
            , aes(x     = date
                  , y     = I # Ia + Ip + Is + Im#+ I
                  , group = .id
            ), alpha = 0.1, colour = "grey30") +
  geom_line(data = (SEIR.sim.f.a %>% filter(.id == "median") %>% filter(date <= as.Date("2020-06-01")))
            , aes(x     = date
                  , y     = I # Ia + Ip + Is + Im#+ I
                  , group = .id
            ), colour = "grey30", lwd = 1.5) +
  geom_vline(xintercept = as.Date("2020-06-01"), linetype = "dashed", lwd = 0.65) +
  theme(legend.title = element_text(size = 14),
        legend.position = "none") +
  annotate("text", x = as.Date("2020-02-01"), y = 3e5, label = "b",
           size = 8, fontface = 2)

gg2 <- gridExtra::grid.arrange(gg2.a, gg2.b, 
                        ncol = 1)
ggsave("figures/Figure2.pdf", plot = gg2, 
       device = "pdf", width = 9, height = 7,
       units = "in")
####
## Figure 3 ----
####
## Simulation details:
## 1) 100 parameter sets from the 200 fits chosen at random.
## 2) 300 stochastic epidemic simulations for 12 test-and-isolate interventions simulated (see figure for strenghts)
## 3) 95% CI for cumulative deaths taken across all parameter sets and stochastic simulations for each intervention scenario

fig3          <- readRDS("manuscript_materials/fig3data.Rds")
ttd           <- fig3[[2]]
SEIR.sim.f.b1 <- fig3[[1]]

(gg.3 <- ggplot(ttd, aes(red_shelt, est, colour = test_and_isolate)) + 
    geom_hline(data = SEIR.sim.f.b1, aes(yintercept = est), linetype = "solid", lwd = 0.5, alpha = 0.5, colour = "black") + 
    geom_hline(data = SEIR.sim.f.b1, aes(yintercept = lwr), linetype = "dotted", lwd = 0.5, alpha = 0.5, colour = "black") + 
    geom_hline(data = SEIR.sim.f.b1, aes(yintercept = upr), linetype = "dotted", lwd = 0.5, alpha = 0.5, colour = "black") + 
    geom_hline(data = SEIR.sim.f.b1, aes(yintercept = lwr2), linetype = "dashed", lwd = 0.5, alpha = 0.5, colour = "black") + 
    geom_hline(data = SEIR.sim.f.b1, aes(yintercept = upr2), linetype = "dashed", lwd = 0.5, alpha = 0.5, colour = "black") + 
    geom_point(lwd = 4) +
    scale_x_continuous(labels = c("40%", "50%", "60%", "70%")) +
    geom_errorbar(aes(ymin = lwr, ymax = upr, colour = test_and_isolate), width = 0.01, lwd = 1) +
    scale_color_manual(
      values = c("lightblue", "royalblue3", "black")
      , labels = c("70%", "80%", "90%")
      , name   = "Test-and-Isolate 
Effectiveness") +
    #  xlab("Proportional of Baseline Infective Contacts in the General Population
    #Under Social Distancing") +
    xlab("Baseline Social Distancing Effectiveness") +    
    ylab("Total Deaths") +
    theme(
      legend.text = element_text(size = 10)
      , legend.title = element_text(size = 12)
      , legend.position = c(0.85, 0.75)))

ggsave("figures/figure3.pdf", plot = gg3,
       device = "pdf", width = 6, height = 5,
       units = "in")
####
## Figure 4 ----
####

# greens? 
# grey for what has already happened?
# fig4_colors = c("dodgerblue4", "firebrick4", "darkgoldenrod4")

## 4A: deaths ----
gg4.a <- ggplot(SEIR.sim.f.cf.a) + 
  geom_line(data = (SEIR.sim.f.cf.a %>% filter(.id != "median" & Scenario != "Reality" & date >= as.Date("2020-03-17")))
            , aes(x     = date
                  , y     = D
                  , group = interaction(.id, Scenario)
                  , colour = Scenario
            ), alpha = 0.1) + 
  geom_line(data = (SEIR.sim.f.cf.a %>% filter(.id == "median" & Scenario != "Reality" & date >= as.Date("2020-03-17")))
            , aes(x      = date
                  , y      = D
                  , colour = Scenario
            ), lwd = 1.5) +
  geom_line(data = (SEIR.sim.f.cf.a %>% filter(.id != "median" & Scenario == "Reality"))
            , aes(x     = date
                  , y     = D
                  , group = interaction(.id, Scenario)
            ), color = "grey30", alpha = 0.1) + 
  geom_line(data = (SEIR.sim.f.cf.a %>% filter(.id == "median" & Scenario == "Reality"))
            , aes(x      = date
                  , y      = D
            ), color = "grey30", lwd = 1.5) +
  # scale_x_date(labels = date_format("%Y-%b"), date_breaks = "1 month") +
  scale_x_date(labels = date_format("%b"), date_breaks = "1 month") +
  scale_colour_manual(
    values = fig4_colors[2:3]
  , labels = c(
   "Shelter in Place
Starting One Week Later"
   , "Test and Isolate 
Starting on March 17, 2020"
   , "Reality"
    )) +
  geom_vline(xintercept = as.Date("2020-03-17"), lty = "dashed") + 
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12)
  , legend.title = element_text(size = 12)
  , plot.title = element_text(size = 12)) +
  ylab("Cumulative Deaths") +
  xlab("") +
  geom_point(data = deaths, aes(date, deaths_cum), lwd = 2, color = "black") +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 5000), labels = c(1, 10, 100, 1000, 5000)) +
    # geom_point(data = deaths, aes(date, deaths), lwd = 2, colour = "forestgreen") +
  theme(legend.title = element_text(size = 12)
        , legend.text = element_text(margin = margin(b = 10, unit = "pt"), size = 10)
        , axis.text.x = element_blank()
        , axis.ticks.x = element_blank()) +
  guides(colour = FALSE) +
  annotate("text", x = as.Date("2020-02-01"), y = 1000, label = "a",
           size = 8, fontface = 2)
# gg4.a

## 4B: cases ----
gg4.b <- SEIR.sim.f.cf.a %>% 
  mutate(Scenario  = factor(Scenario, 
                            levels = c("Shelter in Place One Week Later",
                                       "Test and Isolate starting on March 17",
                                       "Reality"))) %>%
  arrange(Scenario) %>%
  ggplot() + 
  geom_line(data = (SEIR.sim.f.cf.a %>% filter(.id != "median" & Scenario != "Reality" & date >= as.Date("2020-03-17")))
  , aes(x     = date
      , y     = I
      , group = interaction(.id, Scenario)
      , colour = Scenario
    ), alpha = 0.1) + 
  geom_line(data = (SEIR.sim.f.cf.a %>% filter(.id == "median" & Scenario != "Reality" & date >= as.Date("2020-03-17")))
  , aes(x      = date
      , y      = I
      , colour = Scenario
    ), lwd = 1.5) +
  geom_line(data = (SEIR.sim.f.cf.a %>% filter(.id != "median" & Scenario == "Reality"))
            , aes(x     = date
                  , y     = I
                  , group = interaction(.id, Scenario)
            ), color = "grey30", alpha = 0.1) + 
  geom_line(data = (SEIR.sim.f.cf.a %>% filter(.id == "median" & Scenario == "Reality"))
            , aes(x      = date
                  , y      = I
            ), color = "grey30", lwd = 1.5) +
  geom_vline(xintercept = as.Date("2020-03-17"), lty = "dashed") + 
  # scale_x_date(labels = date_format("%Y-%b"), date_breaks = "1 month") +
  scale_x_date(labels = date_format("%b"), date_breaks = "1 month") +
  xlab("Date") + 
  scale_colour_manual(
    values = fig4_colors[2:3]
  , labels = c(
    "Shelter in Place
Starting One Week Later"
   , "Test and Isolate 
Starting on March 17, 2020"
   , "Reality"
    )) +
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12)
  , legend.title = element_text(size = 12)
  , plot.title = element_text(size = 12)) +
  ylab("Cases") +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000), labels = c(1, 10, 100, 1000, 10000)) +
   theme(legend.title = element_text(size = 12)
     , legend.text = element_text(margin = margin(b = 10, unit = "pt"), size = 10)) +
    guides(colour = FALSE) +
  annotate("text", x = as.Date("2020-02-01"), y = 0.4e5, label = "b",
           size = 8, fontface = 2)

## Super ugly way to subtract for the difference histograms. Use
# SEIR.sim.ss.m and SEIR.sim.ss.a for raw

SEIR.sim.ss.m.ow <- SEIR.sim.ss.m[SEIR.sim.ss.m$Scenario == "Shelter in Place One Week Later", ]$D - 
  SEIR.sim.ss.m[SEIR.sim.ss.m$Scenario == "Reality", ]$D
SEIR.sim.ss.m.ow <- data.frame(
  Scenario = "OW"
, D        = SEIR.sim.ss.m.ow)
SEIR.sim.ss.m.ti <- SEIR.sim.ss.m[SEIR.sim.ss.m$Scenario == "Test and Isolate starting on March 17", ]$D - 
  SEIR.sim.ss.m[SEIR.sim.ss.m$Scenario == "Reality", ]$D 
SEIR.sim.ss.m.ti <- data.frame(
  Scenario = "TI"
, D        = SEIR.sim.ss.m.ti)
SEIR.sim.ss.a.ow <- SEIR.sim.ss.a[SEIR.sim.ss.a$Scenario == "Shelter in Place One Week Later", ]$D - 
  SEIR.sim.ss.a[SEIR.sim.ss.a$Scenario == "Reality", ]$D
SEIR.sim.ss.a.ow <- data.frame(
  Scenario = "OW"
, D        = SEIR.sim.ss.a.ow)
SEIR.sim.ss.a.ti <- SEIR.sim.ss.a[SEIR.sim.ss.a$Scenario == "Test and Isolate starting on March 17", ]$D - 
  SEIR.sim.ss.a[SEIR.sim.ss.a$Scenario == "Reality", ]$D 
SEIR.sim.ss.a.ti <- data.frame(
  Scenario = "TI"
, D        = SEIR.sim.ss.a.ti)

SEIR.sim.ss.m2 <- rbind(SEIR.sim.ss.m.ow, SEIR.sim.ss.m.ti)
SEIR.sim.ss.a2 <- rbind(SEIR.sim.ss.a.ow, SEIR.sim.ss.a.ti)

# 4C: change in mortality April 22 ----
gg4.c <- SEIR.sim.ss.a2 %>%
  mutate(Scenario = factor(Scenario, levels = c("R", "OW", "TI"))) %>%
  ggplot(aes(x = D, fill = Scenario)) + 
  geom_histogram(aes(y = ..density..,),
                 position = "identity",
                 colour = "black", bins = 50,
                 lwd = 0.2, alpha = 0.5) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(
    name = "Scenario"
    , values = fig4_colors
      , labels = c(
        "Reality\n"
        ,"Shelter in Place\nStarting One Week Later"
        , "Test and Isolate\nStarting on March 17, 2020"
    ),
   drop = FALSE) +
  theme(legend.title = element_text(size = 12)
        , legend.position = c(0.67, 0.8)
        , legend.text = element_text(margin = margin(b = 10, unit = "pt"), size = 10)) + 
  ylab("Density") +
  xlab("Difference in Total Deaths") +
  ggtitle("Estimated Mortality as of April 22, 2020") +
  annotate("text", x = -220, y = 0.022, label = "c",
           size = 8, fontface = 2)


gg4 <- gridExtra::grid.arrange(gg4.a, gg4.b, gg4.c, 
                        layout_matrix = matrix(c(1,2,3,3), nrow = 2, ncol = 2)) 
ggsave("figures/figure4.pdf", plot = gg4, 
              device = "pdf", width = 10, height = 6,
              units = "in")
                        
## Simple simulation that shows lightswitch working. Requires different parameters

# supplementary figures -----

# Figure S2: violins over time ----
Rt_time <- ldply(all_fits, function(fit){
  return(rbind(data.frame(name = "Reff", 
                          value = as.vector(fit$full_Reff),
                          county = fit$county,
                          fit_sip_eff = fit$fit_sip_eff,
                          days_minus = fit$days_minus,
                          stringsAsFactors = F),
               data.frame(name = "R0", 
                          value = as.vector(fit$full_R0),
                          county = fit$county,
                          fit_sip_eff = fit$fit_sip_eff,
                          days_minus = fit$days_minus,
                          stringsAsFactors = F)))
})
gg.S3 <- Rt_time %>% 
  mutate(date = as.Date("2020-04-22") - as.numeric(days_minus),
         date = format.Date(date, "%b %d")) %>%
  filter(county == "Santa Clara" & fit_sip_eff == "TRUE") %>% 
  ggplot(aes(x = as.factor(date), y = value, fill = name, color = name)) +
  geom_violin(trim = FALSE, color = "black", position = "identity") + 
  geom_hline(yintercept = 1, lty = "dashed") + 
  scale_color_manual(values = fig2_colors[1:2]) + 
  scale_fill_manual(values = fig2_colors[1:2]) + 
  ylab("R0 and Reff estimates") +
  xlab("Last day of data used in fit") 
ggsave("figures/FigureS3.pdf", plot = gg.S3,
       device = "pdf", width = 10, height = 6,
       units = "in")

# Figure S3 recovered/seroprevalence distribution ----
gg.S4 <- SEIR.sim.ss.a %>% filter(Scenario == "Reality") %>% 
  mutate(recov_pct = R/(S + E + I + H + R)) %>%
  ggplot(aes(x = recov_pct, y = ..density..)) + 
  geom_histogram(bins = 50) + 
  geom_density() +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
  xlab("Estimated percent of population recovered on April 22, 2020")
ggsave("figures/figureS4.pdf", plot = gg.S4, 
       device = "pdf", width = 8, height = 6,
       units = "in")

