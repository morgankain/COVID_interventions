## Could simply just wrap this into another loop over pre-sampled parameter values
# for (j in seq_along(param_set)) ...

poss_beta0 <- data.frame(
  beta0  = seq(0.2, 0.7, by = 0.02)
, loglik = 0
)


## Crude brute force creation of the likelihood surface to get a baseline time check -- compare to using mif2 below just for a sanity check
timeoneparam <- system.time({
  
for (i in 1:nrow(poss_beta0)) {
 
  pf <- pfilter(covid, 
        params = c(fixed_params, c(beta0 = poss_beta0[i, "beta0"], N = 1.938e6, E0 = 3)),
        Np = 10000) %>% logLik(.) 
  
  poss_beta0[i, "loglik"] <- pf
  
  if (((i / 4) %% 1) == 0) {print(i)}
}

})

ggplot(poss_beta0, aes(beta0, loglik)) + geom_point()

## with say, 200 samples from the vectors of parameters we want to vary this will only take
paste(round(timeoneparam[3] * 200 / 3600, 3), "hours", sep = " ")

## Compare to using mif2 from pomp
timeoneparam <- system.time({
mifs_local <- covid %>%
        mif2(
          t0 = 15,
          params=c(fixed_params, c(beta0 = 2.5/7, N = 1.938e6, E0 = 3)),
          Np=3000,
          Nmif=50,
          cooling.fraction.50=0.5,
          rw.sd=rw.sd(beta0=0.02)
        )
})

ggout <- mifs_local %>%
  traces() %>%
  melt() %>%
  ggplot(aes(x=iteration,y=value))+
  geom_line()+
  guides(color=FALSE)+
  facet_wrap(~variable,scales="free_y")+
  theme_bw()

mif.l <- ggout$data %>% filter(variable == "loglik" | variable == "beta0") %>% droplevels()
mif.l <- pivot_wider(mif.l, names_from = variable, values_from = value)
ggplot(mif.l, aes(beta0, loglik)) + geom_point()

## Pleasant to see this is actually faster
paste(round(timeoneparam[3] * 200 / 3600, 3), "hours", sep = " ")

## What am i missing here.. why don't these agree on beta0?
plot(mifs_local)


## Regardless of option, from here simulate runs, store results, and calc samples... 

## can give the CI but also a hypercube plot that best illustrates the relationship between each of the variables 
 ## and the summary statistics of interest
