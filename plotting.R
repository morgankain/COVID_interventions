source("ggplot_theme.R")

sim %>% 
  mutate(date = sim_start + day - 1) %>%
ggplot() +
  geom_line(aes(x=date, 
               # y = Is + Im + Ia + Ip,
                y = get("H"),
                group=.id, 
                color = .id == "median")) + 
  scale_x_date(labels = date_format("%Y-%b")) +
  xlab("Date") + ylab("Cases") +
  guides(color=FALSE)+
  scale_color_manual(values=c("#D5D5D3", "#24281A")) 


sim %>% dplyr::group_by(day) %>% summarize(thresh_c = mean(thresh_crossed))

