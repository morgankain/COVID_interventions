####
## Parameters
####

mobility.file      <- "mobility/unfolded_Jun03.rds" ## Mobility data

fitting            <- FALSE   ## Small change in pomp objects if fitting or simulating
use.rds            <- TRUE    ## Run from a previously stored RDS
more.params.uncer  <- FALSE   ## Fit with more (FALSE) or fewer (TRUE) point estimates for a number of parameters
nsim               <- 200     ## Number of epidemic simulations for each parameter set
fit.E0             <- TRUE    ## Was E0 also fit?
fixed.E0           <- !fit.E0 ## Accidental extra param...
usable.cores       <- 1       ## Number of cores to use to fit, for simulation irrelevant

int.beta0_sigma    <- 0.16    ## Heterogeneity k for the Neg Bin that explains variation in individual infectiveness


sir_init.mid       <- FALSE         ## Starts the epidemic from some non-zero timepoint
sir_init.mid.t     <- "2020-05-28"  ## Date to simulate forward from

## Sim and plotting details
sim_length         <- 200     ## How many days to run the simulation
state.plot         <- "D"     ## State variable for plotting (Hospit [H], Death [D], or Cases [C])
plot.log10         <- TRUE    ## Plot on a log10 scale or not
print.plot         <- FALSE

loglik.max         <- TRUE    ## Use just the maximum likelihood fit for one trajectory (TRUE) or many as specified by loglik.thresh
loglik.thresh      <- 2       ## Keep parameter sets with a likelihood within top X loglik units
params.all         <- TRUE    ## Keep all fitted parameters above loglik thresh?...
nparams            <- 50      ## ...if FALSE, pick a random subset for speed

fit.with           <- "D_C"   ## Fit with D (deaths) or H (hospitalizations) -- Need your own data for H -- or both deaths and cases (D_C).
fit_to_sip         <- TRUE    ## Fit beta0 and shelter in place strength simultaneously?
meas.nb            <- TRUE    ## Negative binomial measurement process?
import_cases       <- FALSE   ## Use importation of cases?
fit.minus          <- 0       ## Use data until X days prior to the present

detect.logis       <- TRUE    ## Use logistic detection function. Don't change this, the other option is crap

ci.stoc            <- 0.1     ## Size of the CI to use (0.1 means 80% CI)
plot_vars          <- c("cases", "D")

## Need to specify, but ignored for all scenarios that isn't the tail chopping intervention
int.beta_catch     <- 0.05            ## beta0 values caught by intervention; alternatively, specify by top percent of distribution to trim
int.beta_catch_type<- "pct"           ## if pct, treated as percentile, otherwise, as absoulte value
int.catch_eff      <-  1              ## effectiveness at catching beta0 values above the beta_catch (0 - 1)

