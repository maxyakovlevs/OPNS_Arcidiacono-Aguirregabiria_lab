  params<-list(
    theta = c(10, 0.01),
    s.val = c(1, 4),
    pi.s = .25,
    beta = .99,
    epsilon = 10^-3
    #value for determining x.M
  )
  
  
  
  bellman<-function(df, parameters){
    df%>%mutate(prob_1p = prob_1,
                prob_2p = prob_2,
                value = parameters$theta[1]*engine_id-parameters$theta[2]*mileage+
                  parameters$beta*(log(prob_1p)-log(prob_2p)),
                prob_1 = 1/(1+exp(value)),
                prob_2 = 1-prob_1)
    
  }
  
  iterate<-function(df, parameters){
    temp_delta = 1
    while(temp_delta>10^-3){
      df_upd<-bellman(df,parameters)
      temp_delta = max(abs(df_upd$prob_1 - df_upd$prob_1p))
      print(temp_delta)
      df = df_upd
    }
    return(df)
  }
  
  
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
    saveRDS(paste0('C:/Users/yakov/Downloads/', 'solution_partIII.rds'))
  
  }
  
  
  
