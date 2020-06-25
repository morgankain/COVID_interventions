
int.init <- as.Date(int.init)
sim_end  <- as.Date(sim_end)

    # some catches for problematic interventions
    if(min(int.init) < (max(pomp_data$day) + date_origin)){ # catch if interventions start before the data end
      stop("first intervention must start after data ends")
    } else if(is.unsorted(int.init)){
      stop("interventions must be specified in chronological order")
    } else if(length(int.init) != length(int.type)){
      stop("int.type and int.init must have same length")
    } else if(length(int.init) + 1 != length(int.movement)){
      stop("length(int.init) != length(int.movement[-1])\nfirst element of int.movement specifies movement type after data ends and before interventions begin")
    }

    keep_ints <- int.init < sim_end
    # warning for dropped interventions
    if(sum(keep_ints) < length(int.init)){ 
      print(paste("dropping last ", length(int.init) - sum(keep_ints), " interventions that begin after sim_end"))
    }
    
    int.init     <- int.init[keep_ints] # trim to interventions before the sim_end date
    int.type     <- int.type[keep_ints] # apply same trim to intervention types and movenent
    int.movement <- int.movement[c(TRUE, keep_ints)] # need to add a T for the movement that occurs between data and first intervention
    
    int.phases <- c( 
      nrow(pomp_data), # first phase is the data
      # second phase is until first intervention, rest are between interventions, last is from last intervention to sim_end date
      c(int.init, sim_end) - c(max(pomp_data$day) + date_origin, int.init) 
    )
    
    post_mob <- mean(tail(pomp_data$sip_prop, 3))
    pre_mob  <- mean(head(pomp_data$sip_prop, 7))
    mid_mob  <- (pre_mob + post_mob)/2
    
    mob.covtab <- covariate_table(
      sip_prop           = c(pomp_data$sip_prop, 
                             rep(as.numeric(mapvalues(int.movement, 
                                                      from = c("pre", "mid", "post"),
                                                      to = c(pre_mob, mid_mob, post_mob))), 
                                 int.phases[-1]))       
      , order            = "constant"
      , times            = seq(min(pomp_data$day), as.numeric(sim_end - date_origin), by = 1)
      , iso_mild_level   = rep((c("none", "none", int.type) == "inf_iso")*iso_mild_level, int.phases)
      , iso_severe_level = rep((c("none", "none", int.type) == "inf_iso")*iso_severe_level, int.phases)
      , intervention     = rep(as.numeric(mapvalues(c("none", "none", int.type),
                                         from = c("none", "tail", "inf_iso"),
                                         to   = c(0, 1, 2))), int.phases)   
      , iso_mild_level_post   = iso_mild_level_post
      , iso_severe_level_post = iso_severe_level_post
      , check_int             = c(rep(1, int.phases[1]), rep(0, int.phases[2] + int.phases[3]))
      , thresh_inf            = thresh_inf.val
      )
    
    