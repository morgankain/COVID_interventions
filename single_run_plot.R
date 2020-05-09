### plot of a single run

SEIR.sim <- SEIR.sim %>%
  dplyr::mutate(
      date    = variable_params[i, ]$sim_start + day - 1
    , total_I = Is + Im + Ia + Ip
    )

SEIR.sim %>%
  ggplot() +
  geom_line(
    data = (SEIR.sim %>% filter(.id != "median"))
  , aes(
    x      = date
  , y      = H
  , group  = .id
    ), alpha = 0.5, colour = "grey90") + 
  geom_line(data = (SEIR.sim %>% filter(.id == "median"))
  , aes(  
    x     = date
  , y     = H
    )
  , lwd = 1.5) +
  scale_x_date(labels = date_format("%Y-%b")) +
  scale_y_log10() + 
  xlab("Date") + 
  ylab("Hospitalized") +
  guides(color = FALSE)
