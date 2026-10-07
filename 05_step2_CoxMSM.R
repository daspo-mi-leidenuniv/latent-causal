##########################################################
# Step 2: Marginal structural Cox model (Cox-MSM)        #
##########################################################
# Set the working directory to the root of this repository before running the script
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(ipw)
library(survival)

source("functions/lmm_functions.R")
source("functions/ipw_functions.R")
source("functions/msm_functions.R")
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Weighted dataset and final IPTW/IPDW specifications from 04_step2_IPWselection.R
load("results/IPW_final.RData")
J <- max(dat_weighted$cycle)

# Intervals without follow-up time (tstop <= tstart) are removed when fitting the Cox-MSM
dat_weighted %>% filter(tstop <= tstart) %>% count(event)


###################
# Final Cox-MSM    #
###################
# Stabilized weights truncated at the 5th and 95th percentiles
trunc <- c(0.05, 0.95)
cox_msm <- fit_cox_msm(dat_weighted, trunc = trunc)
summary(cox_msm)

save(cox_msm, file = "results/CoxMSM_final.RData")


#########################################
# RMST differences between strategies   #
#########################################
# Planned schedule: one cycle every 14 days
cycle_start <- 14 * (0:(J - 1))
# Time horizon tau (days)
one_month <- 30.4375
tau_max <- round(72 * one_month)

# Reference strategy: no treatment adjustment
a0 <- rep(0, J)
# Treatment strategies of interest
strategies <- list(
  "One adjustment (cycle 3)"             = c(0, 0, 1, 0, 0),
  "Adjustment at cycles 2 and 4"         = c(0, 1, 0, 1, 0),
  "Adjustment at every cycle"            = rep(1, J)
)

res_RD <- rmst_differences(cox_msm, strategies, a0, cycle_start, tau_max)

# 95% bootstrap confidence intervals (IPTW, IPDW and Cox-MSM re-estimated in each replicate)
B <- 200
set.seed(1234)
boot_RD <- boot_rmst_differences(step2_data, spec_iptw, spec_ipdw, strategies, a0,
                                 cycle_start, tau_max, B = B, trunc = trunc)
CI_RD <- boot_RD %>%
  group_by(strategy, tau) %>%
  summarise(lower = quantile(RD, 0.025), upper = quantile(RD, 0.975), .groups = "drop")

res_RD <- res_RD %>%
  left_join(CI_RD, by = c("strategy", "tau")) %>%
  mutate(strategy = factor(strategy, levels = names(strategies)))

# RMST differences (days) every 12 months
res_RD %>% filter(tau %in% round(one_month * seq(12, 72, by = 12)))

save(res_RD, boot_RD, file = "results/CoxMSM_RMST.RData")


############
# Figure 7 #
############
figure_7 <- ggplot(res_RD, aes(x = tau / one_month, y = RD)) +
  geom_hline(yintercept = 0, color = "gray30", linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = strategy), alpha = 0.2) +
  geom_line(aes(color = strategy), linewidth = 1) +
  facet_wrap(~ strategy) +
  scale_color_manual(values = c("#1C3F94", "goldenrod", "#9e2a2b")) +
  scale_fill_manual(values = c("#1C3F94", "goldenrod", "#9e2a2b")) +
  scale_x_continuous(breaks = seq(0, 72, by = 12)) +
  labs(x = "Time since start of treatment [months]",
       y = "RMST difference vs no adjustment [days]",
       title = "Estimated RMST differences over time") +
  theme_light(base_size = 12) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))
figure_7
ggsave("results/figures/figure_7.pdf", plot = figure_7, width = 12, height = 5)
