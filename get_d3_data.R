  sim %>% 
    select(day,.id,E) %>% 
    pivot_wider(names_from= .id,values_from = E)
  
  write.csv(sim,"sim_data.csv")
