Modeling effort for COVID to predict responses to various intervention strategies

###
NOTE: unfolded_Jun15.Rds is unmodified data from SafeGraph that we can't distribute. You can request an account with them to model COVID to get access to the data.
Here we provide SCC_fuzzed_mobility.rds which is mobility data for Santa Clara with a bit of random noise added so that you can run the model.
Just change: mobility.file <- "unfolded_Jun15.rds" to "SCC_fuzzed_mobility.rds"
###

For cleanest scripts to fit and simulate allowing for all options use:
COVID_fit_cont.R
COVID_simulate_cont.R

COVID_fit.R and COVID_simulate.R are now for an old model version

As of June 19, 2020 the model has moved to multiple boxes per class to move from exponential to gamma distributed rates
