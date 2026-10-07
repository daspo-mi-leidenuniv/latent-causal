##########################################################
# Step 1: Selection of the Latent Markov Model (LMM)     #
##########################################################
# Set the working directory to the root of this repository before running the script
rm(list = ls())
library(dplyr)
library(tidyr)
library(LMest)

source("functions/lmm_functions.R")
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

long_dat <- read.csv2("data/fake_data.csv")
J <- max(long_dat$cycle)

# Balanced long dataset (responses Y are NA after dropout)
lmm_data <- prepare_lmm_data(long_dat, J)


#################################################
# (I) Number of latent states k                 #
#################################################
# LM models with treatment A affecting both initial and transition probabilities,
# without baseline covariates, for k = 1,...,K.
K <- 6
out_k <- fit_lmm(lmm_data, k = 1:K)
info_k <- lmm_info(out_k, step = "Selection of k")
info_k

# Selected number of latent states: minimum BIC
k_final <- info_k$k[which.min(info_k$BIC)]
k_final


#################################################
# (II) Forward selection of baseline covariates #
#################################################
# At each step, each remaining baseline covariate X1-X4 is added either on the initial
# or on the transition probabilities; the model with the lowest BIC is retained.
# The procedure stops when no candidate model reduces the BIC.
init_cov <- NULL
trans_cov <- NULL
best_bic <- info_k$BIC[info_k$k == k_final]
info_fs <- NULL
step <- 0
repeat{
  step <- step + 1
  candidates <- c(lapply(setdiff(X_VARS, init_cov),  function(x) list(init = c(init_cov, x), trans = trans_cov)),
                  lapply(setdiff(X_VARS, trans_cov), function(x) list(init = init_cov, trans = c(trans_cov, x))))
  if(length(candidates) == 0) break
  info_step <- bind_rows(lapply(candidates, function(cand){
    out <- fit_lmm(lmm_data, k = k_final, init_cov = cand$init, trans_cov = cand$trans)
    lmm_info(out, step = paste("Forward selection - step", step), cand$init, cand$trans)
  }))
  print(info_step)
  info_fs <- bind_rows(info_fs, info_step)
  best <- which.min(info_step$BIC)
  if(info_step$BIC[best] >= best_bic) break
  best_bic <- info_step$BIC[best]
  init_cov <- candidates[[best]]$init
  trans_cov <- candidates[[best]]$trans
}

# Supplementary Table B.1: information criteria of the fitted LM models
table_B1 <- bind_rows(info_k, info_fs) %>%
  mutate(across(c(lk, AIC, BIC), ~ round(.x, 2)))
table_B1
write.csv(table_B1, "results/tables/table_B1.csv", row.names = FALSE)


#################################################
# (III) Final LM model                          #
#################################################
lmm_label(init_cov, trans_cov)

# Re-fit the final model with standard errors (out_se) and posterior probabilities V (output)
out_final <- fit_lmm(lmm_data, k = k_final, init_cov = init_cov, trans_cov = trans_cov,
                     out_se = TRUE, output = TRUE)
lmm_final <- out_final$out.single[[1]]

save(lmm_final, lmm_data, init_cov, trans_cov, file = "results/LMM_final.RData")
