#SimData includes the likelihood function and the compiled sampling algorithm  
  
  params<-list(
    theta = c(10, 0.01),
    s.val = c(1, 4),
    pi.s = .25,
    beta = .99,
    epsilon = 10^-3
    #value for determining x.M
  )

  
  likelihood <- function(df){
    group_by(df, engine_id) %>%
      mutate(likelihood = max(10^-12, prob_1*cumprod(prob_2))) %>%
      ungroup()
  }
  
  sampling<-function(samples = 100000){
    x_m <- (params$theta[1]*max(params$s.val)/params$theta[2]) + (log(params$epsilon)/log(1/2))
    
    true_ccp_df<-data.frame(engine_id=c(1,4),
                            prob_1 = 1/2,
                            prob_2 = 1/2)%>%
      expand(mileage = seq(0, x_m), engine_id, prob_1, prob_2) %>%
      mutate(
        prob_1p = 1/2,
        prob_2p = 1-prob_1p)%>%
      iterate(.,params)%>%
      likelihood(.)%>%
      mutate(pi.s = ifelse(engine_id == 1, params$pi.s, 1-params$pi.s))%>%
      group_by(mileage)%>%
      summarise(likelihood_s = sum(pi.s*likelihood))%>%
      sample_n(samples, weight = likelihood_s, replace = T)%>%
      group_by(mileage)%>%
      summarise(count = n())%>%
      saveRDS(paste0('C:/Users/yakov/Downloads/structural estimation/OPNS_Arcidiacono_lab/code/', 'solution_partIII.rds'))
    
  }
  
