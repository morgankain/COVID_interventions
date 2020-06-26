####
## Parameters
####

set.seed(10001)

mobility.file      <- "mobility/unfolded_Jun15.rds"  ## Mobility data


fitting            <- FALSE   ## Small change in pomp objects if fitting or simulating
## TRUE if COVID_fit previously run, FALSE if COVID_fit was just run and global environment is still full
use.rds            <- TRUE    

con_theta          <- FALSE

int.beta0_k        <- 0.16       ## heterogeneity value

## Where to simulate from
sir_init.mid       <- FALSE         ## Starts the epidemic from some non-zero timepoint
sir_init.mid.t     <- "2020-05-28"  ## Date to simulate forward from

## Sim and plotting details
state.plot         <- "D"     ## State variable for plotting (Hospit [H], Death [D], or Cases [C])
plot.log10         <- TRUE    ## Plot on a log10 scale or not
print.plot         <- FALSE

counter.factual    <- FALSE    ## If true do a special analysis that ignores a lot of these other parameters
cf.type            <- "may1"   ## Specifically modeled counterfactual analyses: no_int, delay, may1 coded for now

## Test and Isolate parameters
iso_mild_level        <- 0.2
iso_severe_level      <- 0.2
iso_mild_level_post   <- 0.6
iso_severe_level_post <- 0.6

loglik.thresh      <- 2       ## Keep parameter sets with a likelihood within top X loglik units
loglik.max         <- TRUE
params.all         <- TRUE    ## Keep all fitted parameters above loglik thresh?...
nparams            <- 50      ## ...if FALSE, pick a random subset for speed

fit.with           <- "D_C"   ## Fit with D (deaths) or H (hospitalizations) -- Need your own data for H -- or both deaths and cases (D_C).
fit_to_sip         <- TRUE    ## Fit beta0 and shelter in place strength simultaneously?
meas.nb            <- TRUE    ## Negative binomial measurement process?
import_cases       <- FALSE   ## Use importation of cases?
fit.minus          <- 0       ## Use data until X days prior to the present

detect.logis       <- T

ci.epidemic        <- TRUE    ## Remove all epidemic simulations that dont take off

ci.epidemic_cut    <- 100      ## Criteria of throwing away a stochastic realization as not resulting in an epidemic (total # infected)

nE  <- 3        # must be >= 2
nIa <- 7
nIp <- 2        # must be >= 2
nIm <- nIs <- 5 # must be >= 2
desired.R          <- 2

