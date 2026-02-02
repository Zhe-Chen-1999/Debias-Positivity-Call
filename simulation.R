########### Simulation Studies ########### 

# Evaluate type-I error rates and power for the following three p-values:

# (1) Unadjusted p-value: based on a one-sided test using the primary samples, ignoring the paired control data.
# (2) Maximally adjusted p-value: applies the most conservative adjustment, ensuring validity of the p-value over all plausible batch effects consistent with the paired control samples’ data.
# (3) Minimally adjusted p-value: imposes the minimal adjustment such that the adjusted p-value cannot be falsified by the paired control samples’ data.

###########################################

source("debias_positivity_functions.R")
library(dplyr)
library(xtable)

# Parameters:
## N00: the control sample size at baseline
## N10: the control sample size post-vaccination
## N0: the primary sample size at baseline
## N1: the primary sample size post-vaccination  

simulate_trial <- function(N1, N0, N10, N00, 
                           tp_beta = 500, fp_beta = 2000, C.tp_beta = 0.00029,
                           effect_size = 3, par_ci_alpha = 0.01) {
  
  # Simulate true responder status
  true_responder <- rbinom(1, 1, 0.5)
  
  # Simulate true proportions of positive cells in the primary sample at baseline and post-vaccination
  p_B.tp <- rbeta(1, 1, tp_beta) 
  p_V.tp <- ifelse(true_responder == 0, p_B.tp, min(effect_size*p_B.tp, 1)) # effect size: p_B.tp + rbeta(1, 1, 1000)
  
  # Simulate true proportions of positive cells in the control sample at baseline and post-vaccination
  p_C.B.tp <-  p_C.V.tp <- C.tp_beta # qbeta(0.25, 1, tp_beta)
 
  # Simulate assay false positive rates (moderate batch effect)
  p_B.fp <- rbeta(1, 3, fp_beta)
  p_V.fp <- rbeta(1, 6, fp_beta) 
  
  # Simulate assay false positive rates (no batch effect)
  # p_B.fp <- p_V.fp <- rbeta(1, 1, fp_beta) 
  
  # Simulate assay false positive rates (small batch effect)
  # p_B.fp <- rbeta(1, 1, fp_beta)
  # p_V.fp <- rbeta(1, 2, fp_beta) 
  
  # Simulate assay false positive rates (large batch effect)
  # p_B.fp <- rbeta(1, 1, fp_beta)
  # p_V.fp <- rbeta(1, 5, fp_beta) 
  
  
  # Simulate assay false negative rates (same for both time points)
  p_B.fn <- p_V.fn <- rbeta(1, 1, 5)  
  
  # Expected proportions of observed positive cells in the primary sample
  p_B.obs.p = p_B.tp*(1-p_B.fn) + (1-p_B.tp)*p_B.fp
  p_V.obs.p = p_V.tp*(1-p_V.fn) + (1-p_V.tp)*p_V.fp
  
  # Observed positive count in the primary sample
  n0 = rbinom(1, N0, p_B.obs.p)
  n1 = rbinom(1, N1, p_V.obs.p)
  
  # Expected proportions of observed positive cells in the control sample
  p_C.B.obs.p = p_C.B.tp*(1-p_B.fn) + (1-p_C.B.tp)*p_B.fp
  p_C.V.obs.p = p_C.V.tp*(1-p_V.fn) + (1-p_C.V.tp)*p_V.fp
 
  # Observed positive count in the control sample
  n00 = rbinom(1, N00, p_C.B.obs.p)
  n10 = rbinom(1, N10, p_C.V.obs.p) 
  
  # p-value calculated ignoring misclassification/adjusting for control samples
  res = adjust_pval(n00, N00, n10, N10,
                    n0, N0, n1, N1,
                    alpha = par_ci_alpha, negative = FALSE) 
  
  # p-value calculated using the true misclassification rates
  true_p_value = 1 - pnorm(test.stat(n0, N0, n1, N1, 
                                  p_B.fp = p_B.fp, p_B.fn = p_B.fn, 
                                  p_V.fp = p_V.fp, p_V.fn = p_V.fn))
  
  return(data.frame(responder_status = true_responder,
           true_p_value = true_p_value, 
           unadj_p_value = res$p_unadj, 
           berger_boos_p_value = res$p_adj_BB, 
           min_adj_p_value = res$p_adj,
           n0 = n0, n1 = n1, n00 = n00, n10 = n10,
           p_B.tp = p_B.tp, p_V.tp = p_V.tp, p_C.B.tp = p_C.B.tp, p_C.V.tp = p_C.V.tp,
           p_B.fp = p_B.fp, p_V.fp = p_V.fp, p_B.fn = p_B.fn, p_V.fn = p_V.fn,
           effect_size = effect_size))
}

set.seed(123)
n_sim <- 2000  # number of simulations

res = data.frame()
for (N00 in c(100000)) {
  for (effect_size in c(2,4,10)) {
    for (i in 1:n_sim) {
      res = rbind(res, simulate_trial(N1 = 50000, N0 = 50000, N10 = N00, N00 = N00, 
                                      effect_size = effect_size, C.tp_beta = 0.001))
      cat(i, '\n')
    }
  }
}
write.csv(res, "N00_100000.csv")

res = data.frame()
# primary sample size and control sample size, can vary
for (N00 in c(50000)) {
  for (effect_size in c(2,4,10)) {
    for (i in 1:n_sim) {
      res = rbind(res, simulate_trial(N1 = 50000, N0 = 50000, N10 = N00, N00 = N00, 
                                      effect_size = effect_size, C.tp_beta = 0.001))
      cat(i, '\n')
    }
  }
}
write.csv(res, "N00_50000.csv")

res = data.frame()
for (N00 in c(10000)) {
  for (effect_size in c(2,4,10)) {
    for (i in 1:n_sim) {
      res = rbind(res, simulate_trial(N1 = 50000, N0 = 50000, N10 = N00, N00 = N00, 
                                      effect_size = effect_size, C.tp_beta = 0.005))
      cat(i, '\n')
    }
  }
}
write.csv(res, "N00_10000.csv")

res = data.frame()
for (N00 in c(1000)) {
  for (effect_size in c(2,4,10)) {
    for (i in 1:n_sim) {
      res = rbind(res, simulate_trial(N1 = 50000, N0 = 50000, N10 = N00, N00 = N00, 
                                      effect_size = effect_size, C.tp_beta = 0.03))
      cat(i, '\n')
    }
  }
}
write.csv(res, "N00_1000.csv")

################## Results Analysis ################## 

# Read the simulated data 
res_all = read.csv("N00_100000.csv")

# Delete rows containing NA
sum(apply(res_all, 1, function(row) any(is.na(row))))
rows_with_na = apply(res_all, 1, function(row) any(is.na(row)))
res_all = res_all[!rows_with_na, ]

# Type-I error rates and power associated with each method across different data generating processes.
# Type I error rate is the probability of rejecting the null hypothesis given that it is true. 
res_summary = res_all %>% 
  group_by(effect_size) %>%
  summarize(
    # True p-value
    typeIerror_true = mean((responder_status==0) & (true_p_value < 0.05)),
    power_true = mean((responder_status==1) & (true_p_value < 0.05)),
    
    # Unadjusted p-value (inflated level)
    typeIerror_unadj = mean((responder_status==0) & (unadj_p_value < 0.05)),
    power_unadj = mean((responder_status==1) & (unadj_p_value < 0.05)),
    
    # Maximally adjusted p-value (valid: type I error <= alpha)
    typeIerror_bb = mean((responder_status==0) & (berger_boos_p_value < 0.05)),
    power_bb = mean((responder_status==1) & (berger_boos_p_value < 0.05)),
    
    # Minimally adjusted p-value
    typeIerror_minadj = mean((responder_status==0) & (min_adj_p_value < 0.05)),
    power_minadj = mean((responder_status==1) & (min_adj_p_value < 0.05))
  )

xtable(100*res_summary, digits = 1)
