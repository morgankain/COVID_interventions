##### !!
## Not currently used, but want to fit pmcmc with uncertainty in beta0 following the likelihood profile for beta0
##### !!

## super ugly way of retreiving likelihood plot, come back to this when adding pmcmc
ggout <- mifs_local %>%
  traces() %>%
  melt() %>%
  ggplot(aes(x = iteration, y = value)) +
  geom_line() +
  guides(color = FALSE) +
  facet_wrap(~variable, scales = "free_y") +
  theme_bw()

mif.l  <- ggout$data %>% filter(variable == "loglik" | variable == "beta0") %>% droplevels()
mif.l  <- pivot_wider(mif.l, names_from = variable, values_from = value)
ggout2 <- ggplot(mif.l, aes(beta0, loglik)) + geom_point()

