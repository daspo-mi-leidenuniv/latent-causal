##########################################################
# Auxiliary functions for Step 1: Latent Markov Models   #
##########################################################
# Reference for the LM models:
# - Bartolucci, F., Pandolfi, S., and Pennoni, F. (2017). LMest: An R Package for Latent Markov Models
#   for Longitudinal Categorical Data. Journal of Statistical Software, 81(4):1-38.

# Longitudinal responses (observed manifestations of the latent condition)
Y_VARS <- paste0("Y", 1:8)
# Baseline covariates
X_VARS <- paste0("X", 1:4)


# Function prepare_lmm_data() creates the balanced long dataset required by LMest:
#   - one row per patient and cycle j = 1,...,J;
#   - responses Y are NA for cycles after dropout (j > M);
#   - baseline covariates are carried over and A is set to 0 for cycles after dropout
#     (these values only enter the predicted, non-observed part of the latent process).
prepare_lmm_data <- function(long_dat, J = max(long_dat$cycle)){
  long_dat %>%
    select(id, cycle, all_of(X_VARS), A, all_of(Y_VARS)) %>%
    complete(id, cycle = 1:J) %>%
    group_by(id) %>%
    fill(all_of(X_VARS), .direction = "downup") %>%
    ungroup() %>%
    mutate(A = replace_na(A, 0L)) %>%
    arrange(id, cycle) %>%
    as.data.frame()
}


# Function fit_lmm() fits a LM model with k latent states, where the treatment A always
# affects both the initial and the transition probabilities, and
#   - init_cov:  baseline covariates affecting the initial probabilities
#   - trans_cov: baseline covariates affecting the transition probabilities
fit_lmm <- function(data, k, init_cov = NULL, trans_cov = NULL,
                    nrep = 20, tol1 = 10^-8, seed = 1234, ...){
  latent_formula <- as.formula(paste("~", paste(c(init_cov, "A"), collapse = " + "),
                                     "|", paste(c(trans_cov, "A"), collapse = " + ")))
  out <- lmestSearch(responsesFormula = as.formula(paste(paste(Y_VARS, collapse = " + "), "~ NULL")),
                     latentFormula = latent_formula,
                     version = "categorical",
                     index = c("id", "cycle"),
                     data = data,
                     k = k,
                     nrep = nrep,
                     tol1 = tol1,
                     seed = seed,
                     ...)
  out
}


# Function lmm_label() returns a readable label of the latent formula.
lmm_label <- function(init_cov = NULL, trans_cov = NULL){
  paste0("Initial: ", paste(c(init_cov, "A"), collapse = " + "),
         " | Transition: ", paste(c(trans_cov, "A"), collapse = " + "))
}


# Function lmm_info() extracts the information criteria of a fitted lmestSearch object.
lmm_info <- function(out, step, init_cov = NULL, trans_cov = NULL){
  data.frame(step = step,
             model = lmm_label(init_cov, trans_cov),
             k = out$k,
             np = sapply(out$out.single, function(x) x$np),
             lk = out$lkv,
             AIC = out$Aic,
             BIC = out$Bic,
             row.names = NULL)
}


# Function reference_state() returns the latent state with the lowest average probability
# of observing Y = 1 (i.e., the "healthiest" state), used as reference state R.
reference_state <- function(fit){
  which.min(rowMeans(matrix(fit$Psi[2, , ], nrow = fit$k)))
}


# Function rr_profiles() computes the longitudinal relative risk profiles:
#   RR_s(i,j) = log( P(U_ij = s | Y_i) / P(U_ij = R | Y_i) ),  s != R,
# i.e., the additive log-ratio transform of the posterior probabilities V (Equation (4))
# with respect to the reference latent state R.
# Probabilities are bounded away from 0 (eps) to avoid infinite log-ratios.
rr_profiles <- function(fit, ids, ref = reference_state(fit), eps = 1e-12){
  V <- pmax(fit$V, eps) # n x k x TT
  n <- dim(V)[1]; k <- dim(V)[2]; TT <- dim(V)[3]
  out <- data.frame(id = rep(ids, each = TT), cycle = rep(1:TT, times = n))
  for(s in setdiff(1:k, ref)){
    out[[paste0("RR_", s)]] <- as.vector(t(log(V[, s, ] / V[, ref, ])))
  }
  out
}


# Function state_prevalence() computes the marginal latent state probabilities over cycles
# averaged over all subjects, using the individual-specific initial and transition probabilities.
state_prevalence <- function(fit){
  n <- dim(fit$Piv)[1]; k <- dim(fit$Piv)[2]; TT <- dim(fit$PI)[4]
  piv_ind <- array(NA, dim = c(n, TT, k))
  piv_ind[, 1, ] <- fit$Piv
  for(i in 1:n){
    for(t in 2:TT){
      piv_ind[i, t, ] <- piv_ind[i, t - 1, ] %*% fit$PI[, , i, t]
    }
  }
  prev <- apply(piv_ind, c(2, 3), mean)
  dimnames(prev) <- list(cycle = 1:TT, state = 1:k)
  prev
}


# Function transition_matrix() returns the estimated transition probability matrix
# averaged over all subjects and observed cycles (j >= 2) with treatment A = a.
transition_matrix <- function(fit, lmm_data, a){
  TT <- dim(fit$PI)[4]
  A_obs <- ifelse(is.na(lmm_data[[Y_VARS[1]]]), NA, lmm_data$A) # cycles after dropout excluded
  A_mat <- matrix(A_obs, ncol = TT, byrow = TRUE) # n x TT (data sorted by id, cycle)
  idx <- which(A_mat[, 2:TT, drop = FALSE] == a, arr.ind = TRUE)
  Pi <- Reduce(`+`, lapply(seq_len(nrow(idx)), function(r)
    fit$PI[, , idx[r, 1], idx[r, 2] + 1])) / nrow(idx)
  dimnames(Pi) <- list(from = 1:fit$k, to = 1:fit$k)
  Pi
}


# Function term_labels() returns readable labels of the covariates affecting the initial probabilities.
term_labels <- function(terms, latex = FALSE){
  labels <- terms
  labels[terms == "(Intercept)"] <- "Intercept"
  labels[terms == "A"] <- if(latex) "Initial treatment $A_{\\bullet 1}$" else "Initial treatment A_.1"
  labels
}


# Function latex_table_1() returns the LaTeX code (booktabs) of Table 1, given
#   - Be: matrix of the regression parameters of the initial probabilities (covariates x states);
#   - Pi: list of the transition probability matrices when A = 0 and A = 1.
latex_table_1 <- function(Be, Pi){
  k <- ncol(Be)
  nc <- k + 2
  fmt <- function(x) sprintf("%.2f", x)
  row <- function(...) paste(paste(c(...), collapse = " & "), "\\\\")
  c("\\begin{table}[ht]",
    "\\centering",
    "\\caption{Estimated regression parameters affecting the distribution of the initial probabilities and estimated transition probabilities.}",
    paste0("\\begin{tabular}{ll", strrep("c", k), "}"),
    "\\toprule",
    paste0("\\multicolumn{", nc, "}{l}{\\textbf{Regression parameters for initial probabilities}} \\\\"),
    "\\midrule",
    row("", "$u$", 1:k),
    paste0("\\cmidrule{3-", nc, "}"),
    sapply(1:nrow(Be), function(p)
      row(term_labels(rownames(Be)[p], latex = TRUE), paste0("$\\hat{\\beta}_{", p - 1, "u}$"),
          ifelse(Be[p, ] == 0, "0", fmt(Be[p, ])))),
    "\\midrule",
    paste0("\\multicolumn{", nc, "}{l}{\\textbf{Transition probabilities from $\\tilde{u}$ to $u$}} \\\\"),
    unlist(lapply(names(Pi), function(a)
      c("\\midrule",
        row(paste0("when $A_{ij} = ", a, "$"), "$\\tilde{u} \\backslash u$", 1:k),
        paste0("\\cmidrule{3-", nc, "}"),
        sapply(1:k, function(h) row("", h, fmt(Pi[[a]][h, ])))))),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}")
}
