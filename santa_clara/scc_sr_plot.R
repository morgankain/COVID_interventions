### plot of a single run

SEIR.sim <- SEIR.sim %>%
  dplyr::mutate(
      date    = sim_start + day - 1
    , total_I = Is + Im + Ia + Ip
    )

SEIR.sim %>%
  ggplot() +
  
    geom_line(data = (SEIR.sim %>% filter(.id == "median"))
  , aes(  
    x     = date
  , y     = H_new
  , colour = as.factor(inf_iso))
  , lwd = 1.5) +
  
  geom_line(
    data = (SEIR.sim %>% filter(.id != "median"))
  , aes(
    x      = date
  , y      = H
  , group  = .id
  , colour = inf_iso
    ), alpha = 0.5, colour = "grey90") + 
  
  scale_x_date(labels = date_format("%Y-%b")) +
  scale_y_log10() + 
  xlab("Date") + 
  ylab("Hospitalized") +
  guides(color = FALSE)
