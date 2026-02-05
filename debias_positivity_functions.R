### Helper functions for determining vaccine responders using single-cell assays and paired control samples ###

# Function to calculate the value of the test statistic  

## Parameters:  
## n0: Observed positive counts at baseline  
## N0: Sample size at baseline  
## n1: Observed positive counts post-vaccination  
## N1: Sample size post-vaccination  
## p_B.fp: False positive rate at baseline  
## p_V.fp: False positive rate post-vaccination  
## p_B.fn: False negative rate at baseline  
## p_V.fn: False negative rate post-vaccination 

test.stat <- function(n0 = 8, N0 = 93883, n1 = 43, N1 = 212650,
                      p_B.fp, p_B.fn = 0, p_V.fp, p_V.fn = 0){
  
  # observed positive rate at baseline
  p_B.obs = n0/N0
  
  # estimated true rate of positive outcome at baseline 
  p_B.est = (p_B.obs - p_B.fp)/(1- p_B.fn - p_B.fp)
  
  # observed positive rate at post-vaccination
  p_V.obs = n1/N1
  
  # estimated true rate of positive outcome at post-vaccination 
  p_V.est = (p_V.obs - p_V.fp)/(1- p_V.fn - p_V.fp)
  
  # test statistics
  p_pool.est = (N1*p_V.est + N0*p_B.est)/(N1+N0)
  pooled.var.est = p_pool.est*(1-p_pool.est)*(1/N1 + 1/N0)
  test_stat = (p_V.est - p_B.est)/sqrt(pooled.var.est)
  
  return(test_stat)
}


# (1 - alpha) two-sided confidence set for misclassification rates (p_B.fp, p_B.fn, p_V.fp, p_V.fn)

# Input parameters:
## alpha: Significance level for the two-sided confidence interval (default is 0.05).
## nstep: Number of steps to discretize the range of misclassification rates (default is 21).

## Generic control: no additional constraints on any of the 4 parameters  
CI <- function(alpha = 0.05, n0 = 8, N0 = 93883, n1 = 43, N1 = 212650, nstep = 21){
  
  ci = data.frame()
  
  for (p_B.fp in seq(0, n0/N0-1e-6, length.out = nstep)) { # p_B.obs - p_B.fp > 0
    for (p_V.fp in seq(0, n1/N1-1e-6, length.out = nstep)) { # p_V.obs - p_V.fp > 0
      for (p_B.fn in seq(0, 1-n0/N0 - 1e-6, length.out = nstep)){  
        for (p_V.fn in seq(0, 1-n1/N1 - 1e-6, length.out = nstep)){ 
          test_stat = test.stat(n0, N0, n1, N1, p_B.fp, p_B.fn, p_V.fp, p_V.fn)
          
          if(is.nan(test_stat)){
            ci = rbind(ci, data.frame(p_B.fp = NA, p_V.fp = NA,
                                      p_B.fn = NA, p_V.fn = NA))
          }else if(qnorm(alpha/2) <= test_stat &  test_stat <= qnorm(1-alpha/2)){
              ci <- rbind(ci, data.frame(p_B.fp = p_B.fp, p_V.fp = p_V.fp,
                                         p_B.fn = p_B.fn, p_V.fn = p_V.fn))
          }
        }
      }
    }
  }
  return(ci)
}

## 3 parameters: assume p_B.fn = p_V.fn
CI_three_param <- function(alpha = 0.05, n0 = 8, N0 = 93883, n1 = 43, 
                           N1 = 212650, nstep = 21){
  
  ci = data.frame()
  
  # Define safe upper bounds to prevent negative sequence limits
  # If n0 or n1 is 0, the upper bound becomes 0.
  get_seq <- function(limit, steps) {
    if (limit <= 0) return(0) 
    return(seq(0, limit, length.out = steps))
  }

  # Define the grids
  p_B.fp_seq <- get_seq(n0/N0 - 1e-6, nstep)
  p_V.fp_seq <- get_seq(n1/N1 - 1e-6, nstep)
  p_B.fn_seq <- get_seq(min(1 - n0/N0, 1 - n1/N1) - 1e-6, nstep)
  
  for (p_B.fp in p_B.fp_seq) { # p_B.obs - p_B.fp > 0
    for (p_V.fp in p_V.fp_seq) { # p_V.obs - p_V.fp > 0
      for (p_B.fn in p_B.fn_seq){
        
          test_stat = test.stat(n0, N0, n1, N1, p_B.fp, p_B.fn, p_V.fp, p_B.fn)
          
          if(is.nan(test_stat)){
            ci = rbind(ci, data.frame(p_B.fp = NA, p_V.fp = NA,
                                      p_B.fn = NA, p_V.fn = NA))
          }else if((qnorm(alpha/2) <= test_stat) &  (test_stat <= qnorm(1-alpha/2))){
            ci <- rbind(ci, data.frame(p_B.fp = p_B.fp, p_V.fp = p_V.fp,
                                       p_B.fn = p_B.fn, p_V.fn = p_B.fn))
          }
      }
    }
  }
  return(ci)
}

## 2 parameters: assume p_B.fn = p_V.fn = 0 
CI_two_param <- function(alpha = 0.05, n0 = 8, N0 = 93883, n1 = 43, N1 = 212650,
                side = "Two-sided", nstep = 21){
  
  ci = data.frame()
  
  for (p_B.fp in seq(0, n0/N0-1e-6, length.out = nstep)) {
    for (p_V.fp in seq(0, n1/N1-1e-6, length.out = nstep)) {
      test_stat = test.stat(n0, N0, n1, N1, p_B.fp, p_B.fn = 0, 
                            p_V.fp, p_V.fn = 0)
      
      if(is.nan(test_stat)){
        ci = rbind(ci, data.frame(p_B.fp = NA, p_V.fp = NA,
                                  p_B.fn = NA, p_V.fn = NA))
      }else if(side == "Two-sided"){
        if(qnorm(alpha/2) <= test_stat &  test_stat <= qnorm(1-alpha/2)){
          ci <- rbind(ci, data.frame(p_B.fp = p_B.fp, p_V.fp = p_V.fp,
                                     p_B.fn = 0, p_V.fn = 0) )
        }
      }else{
        if(test_stat <= qnorm(1-alpha)){ ## favoring large p_V.fp and small p_B.fp
          ci <- rbind(ci, data.frame(p_B.fp = p_B.fp, p_V.fp = p_V.fp,
                                     p_B.fn = 0, p_V.fn = 0) )
        }
      }
    }
  }
  return(ci)
}

## Negative control: non-informative confidence interval for p_B.fn and p_V.fn  
CI_two_param_sep <- function(alpha = 0.05, n0 = 8, N0 = 93883, n1 = 43, N1 = 212650,
                         side = "Two-sided", nstep = 21){
  
  CI = tidyr::crossing(p_B.fp = seq(binom.test(x = n0, n = N0, conf.level = 1-alpha/2)$conf.int[1], 
                             binom.test(x = n0, n = N0, conf.level = 1-alpha/2)$conf.int[2], 
                             length.out = 11),
                       p_V.fp = seq(binom.test(x = n1, n = N1, conf.level = 1-alpha/2)$conf.int[1], 
                                 binom.test(x = n1, n = N1, conf.level = 1-alpha/2)$conf.int[2], 
                                 length.out = 11),
                       p_B.fn = seq(0, 1, length.out = nstep),
                       p_V.fn = seq(0, 1, length.out = nstep))
  
  return(CI)
}


# One-sided p-values for all combinations of misclassification rates within the confidence set
calc_pregion <- function(cis, n0 = 31, N0 = 69540, n1 = 85, N1 = 93562){
  pvals = data.frame()
  for (i in 1:nrow(cis)) {
    if(any(is.na(cis[i,]))){
      pval = NA
    }else{
      test_stat = test.stat(n0, N0, n1, N1, 
                          p_B.fp = cis$p_B.fp[i], p_B.fn = cis$p_B.fn[i], 
                          p_V.fp = cis$p_V.fp[i], p_V.fn = cis$p_V.fn[i])
      if(is.nan(test_stat)){
        pval = NA
      }else{
        pval = 1-pnorm(test_stat)
      }
    }
    pvals = rbind(pvals,
                  data.frame(p_B.fp = cis$p_B.fp[i], p_B.fn = cis$p_B.fn[i], 
                             p_V.fp = cis$p_V.fp[i], p_V.fn = cis$p_V.fn[i],
                             pval = pval))
  }
  return(pvals)
}


# Unadjusted, minimally adjusted, and maximally adjusted p-values (all one-sided)
# for identifying vaccine responders, assuming p_B.fn = p_V.fn 

# Parameters:
## n00: Observed positive counts for the control sample at baseline  
## n10: Observed positive counts for the control sample post-vaccination  
## N00: The control sample size at baseline
## N10: The control sample size post-vaccination
## n0: Observed positive counts for the primary sample at baseline  
## n1: Observed positive counts for the primary sample post-vaccination  
## N0: The primary sample size at baseline
## N1: The primary sample size post-vaccination  

adjust_pval <- function(n00 = 8, N00 = 93883, n10 = 43, N10 = 212650,
                        n0 = 31, N0 = 69540, n1 = 85, N1 = 93562,
                        alpha = 0.05, alpha_prime = 0.01) {
  
  if(n0 == 0 & n1 == 0){ # test stat = c'/sqrt(-c) or 0/0 = NaN -> pval = NA ## This person does not show any response at T1 (compared to the paired control)
    primary_p = 1
    p_adj_maximum = 1
    p_adj_minimal = 1
  }else if(n00 == 0 & n10 == 0){ # test stats inverting for nuisance CI = 0/0 = NaN -> CI = NA -> pvals = NA ## perfect assay with zero misclassification rates (no observed positive -> false positive rates = 0 & assume perfect negative control:no true positive in the negative control -> false negative rates = 0 )
    # Unadjusted p-value ignoring misclassification
    primary_p = 1 - pnorm(test.stat(n0, N0, n1, N1, 
                                    p_B.fp = 0, p_B.fn = 0, 
                                    p_V.fp = 0, p_V.fn=0))
    p_adj_maximum = primary_p
    p_adj_minimal = primary_p
  }else{
    # Unadjusted p-value ignoring misclassification
    primary_p = 1 - pnorm(test.stat(n0, N0, n1, N1, 
                                    p_B.fp = 0, p_B.fn = 0, 
                                    p_V.fp = 0, p_V.fn=0))
    
    # Maximally adjusted p-value: usual Berger and Boos approach
    ci_alpha_prime = CI_three_param(alpha = alpha_prime, n00, N00, n10, N10) # Assuming p_B.fn = p_V.fn 
    pval_set_alpha_prime = calc_pregion(ci_alpha_prime, n0, N0, n1, N1)
    pval_set_alpha_prime = pval_set_alpha_prime %>%
      filter(!is.nan(pval)) %>%
      filter(!is.na(pval))
    
    if(nrow(pval_set_alpha_prime)==0){
      p_adj_maximum = NA
    }else{
      p_adj_maximum = max(pval_set_alpha_prime$pval) + alpha_prime
      p_adj_maximum = min(p_adj_maximum, 1)
    }
    
    # Minimally adjusted p-value
    ci_alpha = CI_three_param(alpha = alpha, n00, N00, n10, N10)
    pval_set_alpha = calc_pregion(ci_alpha, n0, N0, n1, N1)
    pval_set_alpha = pval_set_alpha %>%
      filter(!is.nan(pval)) %>%
      filter(!is.na(pval))
    
    if(nrow(pval_set_alpha)==0){
      p_adj_minimal = NA
    }else{
      max_p = max(pval_set_alpha$pval)
      min_p = min(pval_set_alpha$pval)
      p_adj_minimal = ifelse((primary_p >= min_p) & (primary_p <= max_p), primary_p,
                             ifelse(primary_p < min_p, min_p, max_p))
    }
  }
  
  return(list(p_unadj = primary_p,
              p_adj_minimal = p_adj_minimal,
              p_adj_maximum = p_adj_maximum))
}






