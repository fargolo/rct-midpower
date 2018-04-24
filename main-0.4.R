################################################################
# title      : conditional_power
# description: Procedure for calculating statistical power in ongoing
#              clinical trials with dichotomous outcomes
# date       : April - 10, 2018
# version    : 0.4
# r-version  : 3.4.1
# authors    : Guilherme M. Magnavita, Felipe C. Argolo
# usage      : Specify parameters of (1) non-inferiority margin,
#              (2) current size of treatment group, (3) current
#              number of successes in treatment group, (4) current
#              size of control group, (5) current number of 
#              successes in treatment group, (6) number of patients 
#              remaining for randomization, (7) randomization ratio,
#              (8) probability of success in treatment group,
#              (9) probability of success in control group,
#              (10) attrition rate - defaulted to zero
# notes      : probability of success are defaulted to the current
#              ratio in each group
# example    : conditional_power(margin=0.2, trat_npr=24, trat_sucpr=13,
#                  cont_npr=21, cont_sucpr=10, n_to_rand=51,
#                  rand_ratio=0.5)


conditional_power <- function(margin, trat_npr, trat_sucpr,
                              cont_npr, cont_sucpr, n_to_rand,
                              rand_ratio=0.5, 
                              p_trat=trat_sucpr/trat_npr, 
                              p_cont=cont_sucpr/cont_npr, attrition_rate=0){
  
  n_rand_followed <- floor(n_to_rand*(1-attrition_rate))
  power_prob <- 0
  for (trat_n_add in 0:n_rand_followed){ 
    cont_n_add <- n_rand_followed - trat_n_add 
    for (trat_succ_add in 0:trat_n_add){
      for(cont_succ_add in 0:cont_n_add){
        a <- trat_sucpr + trat_succ_add
        b <- cont_sucpr + cont_succ_add
        c <- trat_npr + trat_n_add
        d <- cont_npr + cont_n_add
        margobs_iter <- -(a/c - b/d - 1.64*sqrt((a*(c-a))/((c)^3) + (b*(d-b))/(d)^3))
        if(margobs_iter < margin){
          # Probability of randomizingpatients to treatment
          prob_one <- dbinom(trat_n_add, n_rand_followed, rand_ratio)
          # Probability of successes in treatment
          prob_two <- dbinom(trat_succ_add, trat_n_add, p_trat)
          # Probability of successes in control
          prob_three <- dbinom(cont_succ_add, cont_n_add, p_cont)
          prob_final <- prob_one*prob_two*prob_three
          power_prob <- power_prob + prob_final
        }
      }
    }
  }
  return(power_prob)
}