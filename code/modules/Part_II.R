  #Part II
  
  
  sim_values <- function(num.trials=1000,lambda=3,df.min=.03,df.max=3)
  {tibble(
    agent_i=seq(num.trials)) %>%
      mutate(x = rpois(num.trials,lambda)+1,
             value = runif(num.trials,df.min,df.max)) %>%
      group_by(agent_i,x)%>%
      expand(number_choices = 1:x, value) %>%
      mutate(student=rt(number_choices,value),
             prob = exp(student),
             prob = prob / sum(prob),
             check_probs = near(sum(prob),1),
             check_NA = sum(is.na(prob)),
             lse = log(sum(exp(value/x))),
             lse = exp(lse)/sum(exp(lse)),
             check_lse = near(sum(lse),1),
             check_NA_lse = sum(is.na(lse)),
      )
  }
  
  sim_values()
  saveRDS(paste0('C:/Users/yakov/Downloads/', 'solution_partII.rds'))
  
