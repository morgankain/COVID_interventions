hospit     <- read.csv("ccc_data.csv")
hospit     <- hospit %>% 
  mutate(date = as.Date(REPORT_DATE)) %>% 
  filter(CURRENT_HOSPITALIZED != "NULL") %>% 
  mutate(ch = as.numeric(as.character(CURRENT_HOSPITALIZED))) %>% 
  dplyr::select(date, ch)

SEIR.sim.f.s.f90 <- SEIR.sim.f.s.f90 %>%
  mutate(
    end = 90
  , iso = F
    )

SEIR.sim.f.s.F <- rbind(
  SEIR.sim.f.s.f240
#, SEIR.sim.f.s.f120
, SEIR.sim.f.s.f90
)

SEIR.sim.f.s.F <- SEIR.sim.f.s.F %>% 
  mutate(end = as.factor(end))

SEIR.sim.f.s.t90.50 <- SEIR.sim.f.s.t90.50 %>%
  mutate(
    end   = 90
  , iso   = T
  , sd_m3 = 50
    )

SEIR.sim.f.s.T50 <- rbind(
#  SEIR.sim.f.s.t120.50
 SEIR.sim.f.s.t90.50
, SEIR.sim.f.s.t60.50
)

SEIR.sim.f.s.T50 <- SEIR.sim.f.s.T50 %>% 
  mutate(end = as.factor(end))

gg.T35 <- ggplot(SEIR.sim.f.s.T35[SEIR.sim.f.s.T35$date < "2021-05-15", ]) + 
  geom_line(aes(x = date, y = est, colour = as.factor(end))) + 
  geom_ribbon(aes(x = date, ymin = lwr, ymax = upr, fill = as.factor(end), colour = NA), alpha = 0.2) +
  scale_x_date(labels = date_format("%Y-%b"), date_breaks = "1 month") +
  scale_y_log10() +
  xlab("Date") + 
  ylab("Currently Hospitalized") +
  guides(color = FALSE) +
  scale_color_manual(values = c("dodgerblue4", "brown4", "springgreen4")
  , labels = c("iso_start", "2020-06-15", "2020-07-15")
  , name = "Date shelter in
place relaxes if
infected isolation
possible") +
  scale_fill_manual(
    values = c("dodgerblue4", "brown4", "springgreen4")
  , labels = c("2020-05-16", "2020-06-15", "2020-07-15")
  , name = "Date shelter in place 
relaxes if infected 
isolation is possible") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)
  , legend.title = element_text(size = 12)) +
  ggtitle("Relax shelter in place to 35% of baseline contacts, given that
infected isolation is possible") +
  geom_point(data = hospit , aes(date, ch))

gg.T50 <- ggplot(SEIR.sim.f.s.T50[SEIR.sim.f.s.T50$date < "2021-05-15", ]) + 
  geom_line(aes(x = date, y = est, colour = as.factor(end))) + 
  geom_ribbon(aes(x = date, ymin = lwr, ymax = upr, fill = as.factor(end), colour = NA), alpha = 0.2) +
  scale_x_date(labels = date_format("%Y-%b"), date_breaks = "1 month") +
  scale_y_log10() +
  xlab("Date") + 
  ylab("Currently Hospitalized") +
  guides(color = FALSE) +
  scale_color_manual(values = c("dodgerblue4", "brown4", "springgreen4")
  , labels = c("iso_start", "2020-06-15", "2020-07-15")
  , name = "Date shelter in
place relaxes if
infected isolation
possible") +
  scale_fill_manual(
    values = c("dodgerblue4", "brown4", "springgreen4")
  , labels = c("2020-05-16", "2020-06-15", "2020-07-15")
  , name = "Date shelter in place 
relaxes if infected 
isolation is possible") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)
  , legend.title = element_text(size = 12)) +
  ggtitle("Relax shelter in place to 50% of baseline contacts, given that
infected isolation is possible") +
  geom_point(data = hospit , aes(date, ch))
  

gg.T65 <- ggplot(SEIR.sim.f.s.T65) + 
  geom_line(aes(x = date, y = est, colour = as.factor(end))) + 
  geom_ribbon(aes(x = date, ymin = lwr, ymax = upr, fill = as.factor(end), colour = NA), alpha = 0.2) +
  scale_x_date(labels = date_format("%Y-%b"), date_breaks = "1 month") +
  scale_y_log10() +
  xlab("Date") + 
  ylab("Currently Hospitalized") +
  guides(color = FALSE) +
  scale_color_manual(values = c("dodgerblue4", "brown4", "springgreen4")
  , labels = c("iso_start", "2020-06-15", "2020-07-15")
  , name = "Date shelter in
place relaxes if
infected isolation
possible") +
  scale_fill_manual(
    values = c("dodgerblue4", "brown4", "springgreen4")
  , labels = c("2020-05-16", "2020-06-15", "2020-07-15")
  , name = "Date shelter in place 
relaxes if infected 
isolation is possible") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)
  , legend.title = element_text(size = 12)) +
  ggtitle("Relax shelter in place to 65% of baseline contacts, given that
infected isolation is possible")

gg.F <- ggplot(SEIR.sim.f.s.F[SEIR.sim.f.s.F$date < "2021-05-15", ]) + 
  geom_line(aes(x = date, y = est, colour = as.factor(end))) + 
  geom_ribbon(aes(x = date, ymin = lwr, ymax = upr, fill = as.factor(end), colour = NA), alpha = 0.2) +
  scale_x_date(labels = date_format("%Y-%b"), date_breaks = "1 month") +
  scale_y_log10() +
  xlab("Date") + 
  ylab("Currently Hospitalized") +
  guides(color = FALSE) +
  scale_color_manual(values = c("brown4", "springgreen4", "dodgerblue4")
  , labels = c("2020-05-16", "2020-07-15", "2020-11-12")
  , name = "Date shelter in place removed") +
  scale_fill_manual(
    values = c("dodgerblue4", "brown4", "springgreen4")
  , labels = c("2020-05-16", "2020-07-15", "2020-11-12")
  , name = "Date shelter in place removed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)
  , legend.title = element_text(size = 12)) +
  ggtitle("Shelter in place lifted with no stepwise exit strategy") +
  geom_point(data = hospit , aes(date, ch))

library(gridExtra)
grid.arrange(gg.F, gg.T35, gg.T50, ncol = 1)
