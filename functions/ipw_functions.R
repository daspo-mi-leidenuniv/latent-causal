##########################################################
# Auxiliary functions for Step 2: IPTW and IPDW          #
##########################################################

# Function build_step2_data() combines the observed data with the RR profiles and returns:
#   - dat_A: one row per patient and observed cycle j = 1,...,M (used for IPTW and the Cox-MSM);
#   - dat_D: one row per patient and cycle j = 2,...,min(M+1,J) with the dropout indicator
#            D = 1 if the patient is not observed at cycle j (dropout after cycle j-1), 0 otherwise.
# Time-varying covariates at cycle j:
#   - RR_s:     RR profiles at cycle j-1 (latent confounders preceding treatment at cycle j)
#   - A_lag:    treatment at cycle j-1
#   - cumA_lag: cumulative number of treatment adjustments up to cycle j-1
#   - cumA:     cumulative number of treatment adjustments up to cycle j (Cox-MSM exposure)
build_step2_data <- function(long_dat, rr, J = max(long_dat$cycle)){
  rr_lag <- rr %>% mutate(cycle = cycle + 1)

  dat_A <- long_dat %>%
    select(id, all_of(X_VARS), cycle, tstart, tstop, event, M, A) %>%
    arrange(id, cycle) %>%
    group_by(id) %>%
    mutate(A_lag = lag(A, default = 0),
           cumA = cumsum(A),
           cumA_lag = lag(cumA, default = 0)) %>%
    ungroup() %>%
    left_join(rr_lag, by = c("id", "cycle"))

  dropout_rows <- dat_A %>%
    filter(cycle == M, M < J) %>%
    select(id, all_of(X_VARS), cycle, M, A, cumA) %>%
    mutate(cycle = cycle + 1, A_lag = A, cumA_lag = cumA, A = NA, cumA = NA) %>%
    left_join(rr_lag, by = c("id", "cycle"))

  dat_D <- bind_rows(dat_A %>% filter(cycle >= 2) %>% mutate(D = 0),
                     dropout_rows %>% mutate(D = 1)) %>%
    select(id, all_of(X_VARS), cycle, D, A_lag, cumA_lag, starts_with("RR_")) %>%
    arrange(id, cycle)

  list(dat_A = dat_A, dat_D = dat_D)
}


# Function ipw_specs() defines the candidate specifications for IPTW and IPDW
# (den: covariates of the denominator model; num: covariates of the numerator model).
ipw_specs <- function(rr_vars, baseline = X_VARS){
  list(
    "IPW0" = list(den = rr_vars,                           num = NULL),
    "IPW1" = list(den = c(baseline, rr_vars),              num = NULL),
    "IPW2" = list(den = c(rr_vars, "A_lag"),               num = "A_lag"),
    "IPW3" = list(den = c(baseline, rr_vars, "A_lag"),     num = "A_lag"),
    "IPW4" = list(den = c(rr_vars, "cumA_lag"),            num = "cumA_lag"),
    "IPW5" = list(den = c(baseline, rr_vars, "cumA_lag"),  num = "cumA_lag")
  )
}


# Function ipw_formula() builds the right-hand side formula of a numerator/denominator model.
ipw_formula <- function(x) reformulate(if(length(x) == 0) "1" else x)


# Function compute_iptw() computes the stabilized IPTW for each patient and cycle (package ipw):
#   - cycle 1: ipwpoint(), treatment model conditional on baseline covariates only;
#   - cycles j >= 2: ipwtm() with type = "all", treatment model defined by spec;
#     ipwtm() already returns the cumulative product of the ratios over cycles 2,...,j;
# the weight at cycle j >= 2 is then multiplied by the weight at cycle 1.
compute_iptw <- function(dat_A, spec, baseline = X_VARS){
  dat_A <- as.data.frame(dat_A %>% arrange(id, cycle))
  c1 <- dat_A$cycle == 1
  
  w1 <- do.call(ipwpoint, list(exposure = as.name("A"),
                               family = "binomial",
                               link = "logit",
                               numerator = ~ 1,
                               denominator = ipw_formula(baseline),
                               data = dat_A[c1, ],
                               trace = FALSE))$ipw.weights
  
  w2 <- do.call(ipwtm, list(exposure = as.name("A"),
                            family = "binomial",
                            link = "logit",
                            numerator = ipw_formula(spec$num),
                            denominator = ipw_formula(spec$den),
                            id = as.name("id"),
                            timevar = as.name("cycle"),
                            type = "all",
                            data = dat_A[!c1, ],
                            trace = FALSE))$ipw.weights
  
  dat_A$w <- NA
  dat_A$w[c1] <- w1
  dat_A$w[!c1] <- w2
  dat_A %>%
    group_by(id) %>%
    mutate(sw_iptw = ifelse(cycle == 1, w, w * w[cycle == 1])) %>%
    ungroup() %>%
    pull(sw_iptw)
}


# Function compute_ipdw() computes the stabilized IPDW for each patient and observed cycle (package ipw):
# ipwtm() with type = "first" returns, at cycle j >= 2, the cumulative product up to cycle j of the ratios
# P(D = 0 | num covariates) / P(D = 0 | den covariates); the weight is equal to 1 at cycle 1.
compute_ipdw <- function(dat_A, dat_D, spec){
  dat_D <- as.data.frame(dat_D %>% arrange(id, cycle))
  dat_D$sw_ipdw <- do.call(ipwtm, list(exposure = as.name("D"),
                                       family = "binomial",
                                       link = "logit",
                                       numerator = ipw_formula(spec$num),
                                       denominator = ipw_formula(spec$den),
                                       id = as.name("id"),
                                       timevar = as.name("cycle"),
                                       type = "first",
                                       data = dat_D,
                                       trace = FALSE))$ipw.weights
  dat_A %>%
    arrange(id, cycle) %>%
    left_join(dat_D %>% filter(D == 0) %>% select(id, cycle, sw_ipdw), by = c("id", "cycle")) %>%
    mutate(sw_ipdw = replace_na(sw_ipdw, 1)) %>%
    pull(sw_ipdw)
}


# Function add_weights() adds the IPTW, IPDW and final stabilized weights sw = sw_iptw * sw_ipdw.
add_weights <- function(step2_data, spec_iptw, spec_ipdw){
  step2_data$dat_A %>%
    arrange(id, cycle) %>%
    mutate(sw_iptw = compute_iptw(step2_data$dat_A, spec_iptw),
           sw_ipdw = compute_ipdw(step2_data$dat_A, step2_data$dat_D, spec_ipdw),
           sw = sw_iptw * sw_ipdw)
}
