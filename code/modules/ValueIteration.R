#ValueIteration - includes a function for Bellman contractor and iteration algorithm


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
  