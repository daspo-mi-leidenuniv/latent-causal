##########################################################
# Auxiliary functions for Step 2: Cox-MSM                #
##########################################################

# Function fit_cox_msm() fits the weighted marginal structural Cox model
#   h(t | cumA(t)) = h0(t) * exp(beta * cumA(t)),
# where cumA(t) is the cumulative number of treatment adjustments up to the current cycle.
# Stabilized weights are truncated at the given percentiles (trunc).
# Intervals with tstop <= tstart carry no follow-up time and are removed.
fit_cox_msm <- function(dat_w, trunc = c(0.05, 0.95)){
  dat_cox <- dat_w %>% filter(tstop > tstart)
  q <- quantile(dat_cox$sw, probs = trunc)
  dat_cox <- dat_cox %>% mutate(sw_trunc = pmin(pmax(sw, q[1]), q[2]))
  coxph(Surv(tstart, tstop, event) ~ cumA, data = dat_cox,
        weights = sw_trunc, cluster = id)
}


# Function rmst_strategy() computes the restricted mean survival time RMST(tau), tau = 1,...,tau_max
# (in days), under the treatment strategy a = (a_1,...,a_J), where cycle j starts at cycle_start[j]
# and the last cycle lasts until the end of follow-up:
#   S_a(t) = exp( - sum_{s <= t} dH0(s) * exp(beta * cumA_a(s)) ),   RMST_a(tau) = int_0^tau S_a(t) dt.
# Since event times are recorded in days, S_a(t) is constant on [d, d+1) and the integral is exact.
rmst_strategy <- function(fit, a, cycle_start, tau_max){
  bh <- basehaz(fit, centered = FALSE)
  dH0 <- diff(c(0, bh$hazard))
  j <- pmax(findInterval(bh$time, cycle_start, left.open = TRUE), 1)
  H <- cumsum(dH0 * exp(coef(fit)[["cumA"]] * cumsum(a)[j]))
  days <- 0:(tau_max - 1)
  S <- exp(-c(0, H)[findInterval(days, bh$time) + 1])
  cumsum(S)
}


# Function rmst_differences() returns the RMST differences RMST_a(tau) - RMST_a0(tau)
# between each strategy in the list strategies and the reference strategy a0.
rmst_differences <- function(fit, strategies, a0, cycle_start, tau_max){
  r0 <- rmst_strategy(fit, a0, cycle_start, tau_max)
  bind_rows(lapply(names(strategies), function(s)
    data.frame(strategy = s,
               tau = 1:tau_max,
               RD = rmst_strategy(fit, strategies[[s]], cycle_start, tau_max) - r0)))
}


# Function boot_sample() resamples patients with replacement and assigns new unique ids,
# keeping all the cycles of each sampled patient (in both dat_A and dat_D).
boot_sample <- function(step2_data){
  ids <- unique(step2_data$dat_A$id)
  map <- data.frame(id = sample(ids, length(ids), replace = TRUE), boot_id = seq_along(ids))
  resample <- function(d) inner_join(map, d, by = "id", relationship = "many-to-many") %>%
    mutate(id = boot_id) %>% select(-boot_id)
  list(dat_A = resample(step2_data$dat_A), dat_D = resample(step2_data$dat_D))
}


# Function boot_rmst_differences() computes B bootstrap replicates of the RMST differences.
# In each replicate, IPTW, IPDW and the Cox-MSM are re-estimated on the bootstrap sample,
# while the RR profiles from Step 1 are kept fixed.
boot_rmst_differences <- function(step2_data, spec_iptw, spec_ipdw, strategies, a0,
                                  cycle_start, tau_max, B = 200, trunc = c(0.05, 0.95)){
  bind_rows(lapply(1:B, function(b){
    dat_w_b <- add_weights(boot_sample(step2_data), spec_iptw, spec_ipdw)
    fit_b <- fit_cox_msm(dat_w_b, trunc = trunc)
    rmst_differences(fit_b, strategies, a0, cycle_start, tau_max) %>% mutate(b = b)
  }))
}
