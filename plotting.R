source("ggplot_theme.R")

sim %>% 
  mutate(date = sim_start + day - 1) %>%
  # filter(.id == 3) %>%
ggplot() +
  geom_line(aes(x=date, 
                # y = Is + Im + Ia + Ip,
                y = H,
                group=.id, 
                color = .id == "median")) + 
  scale_x_date(labels = date_format("%Y-%b")) +
  guides(color=FALSE)+
  scale_color_manual(values=c("#D5D5D3", "#24281A")) 
