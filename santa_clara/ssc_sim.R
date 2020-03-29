## Could simply just wrap this into another loop over pre-sampled parameter values
# for (j in seq_along(param_set)) ...

poss_beta0 <- data.frame(
  beta0  = seq(0.2, 0.7, by = 0.02)
, loglik = 0
)

timeoneparam <- system.time({
  
for (i in 1:nrow(poss_beta0)) {
 
  pf <- pfilter(covid, 
        params = c(fixed_params, c(beta0 = poss_beta0[i, "beta0"], N = 100000, E0 = 10)),
        Np = 5000) %>% logLik(.) 
  
  poss_beta0[i, "loglik"] <- pf
  
  if (((i / 4) %% 1) == 0) {print(i)}
}

})

ggplot(poss_beta0, aes(beta0, loglik)) + geom_point()

## with say, 200 samples from the vectors of parameters we want to vary this will only take
paste(round(timeoneparam[3] * 200 / 3600, 3), "hours", sep = " ")

## then from here simulate runs, store results, and calc samples... 

## can give the CI but also a hypercube plot that best illustrates the relationship between each of the variables 
 ## and the summary statistics of interest
