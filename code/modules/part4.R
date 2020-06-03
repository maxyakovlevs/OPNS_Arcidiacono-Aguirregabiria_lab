

  
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
                value = parameters$theta[1]*s-parameters$theta[2]*x+
                  parameters$beta*(log(prob_1p)-log(prob_2p)),
                value = ifelse(value< -(10^8),10^-6,value),
                prob_1 = 1/(1+exp(value)),
                prob_2 = 1-prob_1)
    
  }
  
  iterate<-function(df, parameters){
    temp_delta = 1
    while(temp_delta>10^-6){
      df_upd<-bellman(df,parameters)
      temp_delta = max(abs(df_upd$prob_1 - df_upd$prob_1p))
      print(temp_delta)
      df = df_upd
    }
    return(df)
  }


sim_data<-readRDS('C:/Users/yakov/Downloads/sim_data.RDS')

  test_fn <- function(df = sim_data, parameters = params) {
    sim_data %>%
      rbind(c(1, max(filter(., s == 1)$x) + 1, 1)) %>%
      rbind(c(2, max(filter(., s == 2)$x) + 1, 1)) %>%
        group_by(s) %>%
          mutate(prob_1 = n / sum(n)) %>%
          mutate(prob_1 = cumsum(lag(prob_1, default = 0, order_by = x))) %>%
        ungroup() %>%
          mutate(prob_2 = 1 - prob_1) %>%
          arrange(s, x) %>%
          iterate(., parameters) %>%
          mutate(log_prob_1 = log(prob_1)) %>%
          mutate(log_prob_2 = log(prob_1)) %>%
          select(s, x, prob_1, prob_2, log_prob_1, log_prob_2) %>%
        group_by(s) %>%
          mutate(log.like = log_prob_2 + cumsum(lag(log_prob_1, default = 0, order_by = x))) %>%
          summarise(., sum_like = sum(log.like)) %>%
          select(s, sum_like)
  }
  


  start <- list(
    theta = c(1000, 1000),
    s.val = c(1, 4),
    pi.s = .25,
    beta = .99,
    epsilon = 10^-6
    # value for determining x.M
  )


  optim(
    par = start, 
      fn = test_fn, 
        method = "L-BFGS-B", 
          control = list(fnscale = -1))

  num.tries <- 5
  
  num.tries %>%
    seq() %>%
    plyr::llply(
      function(l) {
        optim(
          par = start, 
            fn = test_fn, 
            method = "L-BFGS-B", 
              control = list(fnscale = -1))
      },
      .parallel = TRUE
    )
