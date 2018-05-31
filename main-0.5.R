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

require(magrittr)
require(purrr)
require(dplyr)

margobs <- function(a,b,c,d){
  -(a/c - b/d - 1.64*sqrt((a*(c-a))/((c)^3) + (b*(d-b))/(d)^3))}

joint_binom <- function(trat_succ_add,cont_succ_add,
                        trat_n_add,cont_n_add,
                        rand_ratio,n_rand_followed,
                        p_trat,p_cont){
  # Probability of randomizingpatients to treatment
  prob_one <- dbinom(trat_n_add, n_rand_followed, rand_ratio)
  # Probability of successes in treatment
  prob_two <- dbinom(trat_succ_add, trat_n_add, p_trat)
  # Probability of successes in control
  prob_three <- dbinom(cont_succ_add, cont_n_add, p_cont)
  joint_prob <- prob_one*prob_two*prob_three
  return(joint_prob)}

conditional_power <- function(margin, trat_npr, trat_sucpr,
                              cont_npr, cont_sucpr, n_to_rand,
                              rand_ratio=0.5, 
                              p_trat=trat_sucpr/trat_npr, 
                              p_cont=cont_sucpr/cont_npr, attrition_rate=0){
  n_rand_followed <- floor(n_to_rand*(1-attrition_rate))
  
  c <- trat_npr:(trat_npr+n_rand_followed)
  d <- cont_npr:(cont_npr+n_rand_followed)
  a <- trat_sucpr:(trat_sucpr+n_rand_followed)
  b <- cont_sucpr:(cont_sucpr+n_rand_followed)
  tot_rand <- trat_npr+cont_npr+n_rand_followed
  
  scenario_sims <- data.frame(expand.grid(a,b,c,d)) %>% setNames(c("a","b","c","d")) %>%
    mutate(marg_obs = margobs(a,b,c,d)) %>% # calculate margins
    subset(marg_obs < margin & c + d == tot_rand) %>% # non-inf scenarios
    mutate(trat_s=a - trat_sucpr, cont_s=b - cont_sucpr,
           trat_n=c - trat_npr, cont_n=d - cont_npr) %>% # retrieving args for joint binom 
    mutate(p=joint_binom(trat_s, cont_s, trat_n, cont_n, #joint probability
                         rand_ratio,n_rand_followed,p_trat,p_cont))
  sum(scenario_sims$p)} 