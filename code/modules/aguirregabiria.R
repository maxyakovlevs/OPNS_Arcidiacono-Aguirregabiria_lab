params<-list(
  Q = 20,
  mu = 1,
  sigma = 1,
  p.r = 5,
  p.w = 1, 
  alpha = 0.3,
  eta  = 4,
  beta = 0.95,
  num_grid =100,
  periods = 10,
  epsilon = 10^-3
)

sampling <- function(params){


  V_hat_value <- rep(0, params$num_grid)
  
  inc_val <- function(v.1 ,v.2)
  {
    log(exp(v.1)+exp(v.2))-digamma(1)
  }
  
  g <- function(i,s){
    v.1 <- params$p.r*s - params$p.w*(params$Q-i+s) - params$alpha*(i-s) - params$eta
    v.2 <- params$p.r*s - params$alpha*(i-s)
    inc_val(v.1, v.2) * dlnormTrunc(s, params$mu, params$sigma, min=0, max=i)
  }
  
  delta <- 1
  while(delta > params$epsilon){
    
    f <- function(i){
      integrate(g, i=i, lower=0, upper=i) %>%
        .$value
    }
    
    V_hat <<- Vectorize(chebappxf(f, dims=params$num_grid, intervals=list(c(0,params$Q))))
    
    g <- function(i,s){ 
      v.1 <- params$p.r*s - params$p.w*(params$Q-i+s) - params$alpha*(i-s) - params$eta + params$beta*V_hat(params$Q)
      v.2 <- params$p.r*s - params$alpha*(i-s) + params$beta*V_hat(i-s)
      inc_val(v.1, v.2) * dlnormTrunc(s, params$mu, params$sigma, min=0, max=i)
    }
    
    V_hat_upd <- V_hat(seq(from=0, to=params$Q, length.out=params$num_grid))
    delta <- max(abs(V_hat_value-V_hat_upd))
    V_hat_value <- V_hat_upd
  }
  
  
  CCP <- function(df){   
    df %>%
      mutate(
        v.1 = params$p.r*s - params$p.w*(params$Q-i+s) - params$alpha*(i-s) - params$eta + params$beta *V_hat(params$Q),
        v.2 = params$p.r*s - params$alpha*(i-s) + params$beta * V_hat(i-s),                    
        f = v.2 - v.1,
        prob_order = pmax(.Machine$double.xmin, 1/(1+exp(f))),
        decision = ifelse(runif(1) > prob_order, 1, 0)
      )  
  }
  

  sim.data <-
    data.frame(
      d = rlnorm(params$periods, params$mu, params$sigma),
      i = 10) %>%
    mutate(
      s = min(d, i),
      prob = 0.5,
      decision = NA
    )%>%
      mutate(
        i = ifelse(decision==1, params$Q, i-s),
        s = min(i, d),
        t = seq(params$periods))%>%
      rowwise()%>%
      do(CCP(.))%>%
      saveRDS(paste0('C:/Users/yakov/Downloads/sim_data1.RDS'))
}


df<-readRDS('sim_data.RDS', refhook=NULL)
  

estimation <- function(df,
                      params) {
  
  
   df %>%
    mutate(
      i.s = i - s,
      CCP.est = predict(npreg(decision ~ i.s))
      )
  
  MLE <- function(params){
    
    mu <- params$mu
    sigma <- params$sigma
    
    df %>%
      group_by(i,s) %>%
      mutate(inc.val = ifelse(i!=0, dlnormTrunc(x = s, mu, sigma, min = 0, max = i), NA)) %>%
      ungroup() %>%
      mutate(p_hat = ifelse(decision == 1, CCP.est, 1 - CCP.est)) %>%
      filter(!is.na(inc.val)) %>%
      summarise(sum(log(p_hat * inc.val))) %>%
      prod(-1) %>%
      as.numeric
  }
  
  
  start = data.frame(1,0)
  est.ms <- optim(start[1:2], likelihood, method = "L-BFGS-B", lower = c(-Inf, 0), control = list(trace = TRUE))$par
  

  estimate_f <- function(params){
    
    f.hat <- function(i,s, params){
      decision <- ifelse(params$Q - i + s > 0, 1, 0)
      i.s <- i - s
      func <- (params$p.r*s - params$p.w*(params$Q-i+s) - params$alpha*(i-s) - decision * params$eta - log(df$CCP.est)) * dlnormTrunc(s, params_n[1], params_n[2], min = 0, max = i)
    }
    
    f <- function(i){
      integrate(
        f.hat,
        i = i,
        pr = params$p.r,
        alpha = params$alpha,
        eta = params$eta,
        lower = 0,
        upper = i
      ) %>%
        .$value
    }
    
    f.tilde <- Vectorize(chebappxf(f(i), dims = num_grid, intervals = list(c(0,params$Q))))
    
    v <- function(x){
      p.w*(params$Q-x) + params$eta + params$beta*(f.tilde(x) - f.tilde(params$Q))
    }
    
    loglikelihood <- df %>%
      mutate(
        p.v0 = 1/(1 + exp(-v(i.s))),
        pv.1 = 1/(1 + exp(v(i.s))),
        like = ifelse(decision == 1, p.v0, p.v1),
        log.like <- pv.1 + cumsum(lag(p.v0, default=0, order_by=d))
      ) %>%
      summarise(sum(log.like)) %>%
      prod(-1) %>%
      as.numeric()
    
  }
  start = data.frame(1,1,1)
  est.rest <- optim(start, estimate_f, method = "L-BFGS-B", lower = c(0,0, -Inf), control = list(trace = TRUE))
  
  return(est.ms, est.rest )
}
