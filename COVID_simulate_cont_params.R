####
## Parameters
####

mobility.file      <- "mobility/unfolded_Jun15.rds"  ## Mobility data

fitting            <- FALSE   ## Small change in pomp objects if fitting or simulating
use.rds            <- TRUE    ## Run from a previously stored RDS
usable.cores       <- 1       ## Number of cores to use to fit, for simulation irrelevant

# date_origin        <- as.Date("2019-12-31")

sir_init.mid       <- FALSE         ## Starts the epidemic from some non-zero timepoint
sir_init.mid.t     <- "2020-05-28"  ## Date to simulate forward from

## Sim and plotting details
nsim               <- 100     ## Number of epidemic simulations for each parameter set
plot.log10         <- TRUE    ## Plot on a log10 scale or not
plot.median        <- TRUE    ## Use the median? or the mean?

loglik.max         <- TRUE    ## Use just the maximum likelihood fit for one trajectory (TRUE) or many as specified by loglik.thresh
loglik.thresh      <- 2       ## Keep parameter sets with a likelihood within top X loglik units
loglik.num         <- NA
params.all         <- TRUE    ## Keep all fitted parameters above loglik thresh?...
nparams            <- 50      ## ...if FALSE, pick a random subset for speed

fit.with           <- "D_C"   ## Fit with D (deaths) or H (hospitalizations) -- Need your own data for H -- or both deaths and cases (D_C).
fit_to_sip         <- TRUE    ## Fit beta0 and shelter in place strength simultaneously?
meas.nb            <- TRUE    ## Negative binomial measurement process?
import_cases       <- FALSE   ## Use importation of cases?
fit.minus          <- 0       ## Use data until X days prior to the present

ci.stoc            <- 0.1     ## Size of the CI to use (0.1 means 80% CI)
ci.epidemic        <- TRUE    ## whether to limit to simulations where epidmeics occurs. epidemic defined by > 2 * E_init eventually are in the recovered class
ci.epidemic_cut    <- 100     ## Criteria of throwing away a stochastic realization as not resulting in an epidemic (total # infected)
plot_vars          <- c("cases", "D")

## Intervention details
counter.factual     <-  FALSE
sim_end             <- "2020-12-31"     ## How many days to run the simulation
int.init            <- c("2020-07-01", "2020-08-01") # need to be chronologically in order, any after sim_end get cut
int.type            <- c("none", "none") # options are none, tail, inf_iso; starts at corresponding int.init
int.movement        <- c("post", "post", "post") # movement returned to after data ends, then movement for each intervention period

iso_mild_level      <- 1 #0.2
iso_severe_level    <- 1 #0.2
int.beta_catch      <- 0.0            ## beta0 values caught by intervention; alternatively, specify by top percent of distribution to trim
int.beta_catch_type <- "pct"           ## if pct, treated as percentile, otherwise, as absoulte value
int.catch_eff       <-  0              ## effectiveness at catching beta0 values above the beta_catch (0 - 1)
int.beta0_k         <- 0.16    ## Heterogeneity k for the Neg Bin that explains variation in individual infectiveness
