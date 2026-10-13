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


# Function surv_strategy() computes the marginal survival probability under the treatment strategy
# a = (a_1,...,a_J), where cycle j starts at cycle_start[j] and the last cycle lasts until the end
# of follow-up:
#   S_a(t) = exp( - sum_{s <= t} dH0(s) * exp(beta * cumA_a(s)) ),
# with H0 the baseline cumulative hazard of the Cox-MSM and cumA_a(s) the cumulative number of
# treatment modifications under strategy a up to the cycle containing time s.
# The curve is computed up to the time horizon t_max.
surv_strategy <- function(fit, a, cycle_start, t_max = Inf){
  bh <- basehaz(fit, centered = FALSE)
  bh <- bh[bh$time <= t_max, ]
  dH0 <- diff(c(0, bh$hazard))
  j <- pmax(findInterval(bh$time, cycle_start, left.open = TRUE), 1)
  H <- cumsum(dH0 * exp(coef(fit)[["cumA"]] * cumsum(a)[j]))
  data.frame(time = c(0, bh$time), surv = exp(-c(0, H)))
}
