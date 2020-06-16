####
## Simulation: Continue with current workplace closures (closed)
####
counter.factual    <- FALSE
int.movement       <- "post"
int.type           <- "none"
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" ## Doesnt matter, ignored for all int.type != "tail"

####
## Simulation: Immediately relax all current restrictions on workplaces (open)
####
counter.factual    <- FALSE
int.movement       <- "pre"
int.type           <- "none"
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" ## Doesnt matter, ignored for all int.type != "tail"

####
## Simulation: relax current social distancing 2 weeks after peak. Still need to code peak14 into the pomp object like hospital threshold was
####
counter.factual    <- FALSE
int.movement       <- "post"
int.type           <- "days_post_peak"
int.post_peak_thresh<- 14
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" ## Doesnt matter, ignored for all int.type != "tail"

####
## Simulation: relax social distancing when the number of new daily cases is at 1% of peak. Also need to track for this int as well 
####
counter.factual    <- FALSE
int.movement       <- "post"
int.type           <- "one_perc"
int.init           <- "2020-06-08"
int.end            <- "2020-08-01" ## Doesnt matter, ignored for all int.type != "tail"
