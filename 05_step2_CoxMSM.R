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


###############################################
# Marginal survival under treatment strategies #
###############################################
# Planned schedule: one cycle every 14 days (the last cycle lasts until the end of follow-up)
cycle_start <- 14 * (0:(J - 1))

# Treatment strategies of interest
strategies <- list(
  "Never modified"  = rep(0, J),                    # (0,...,0)
  "Alternating"     = rep(c(0, 1), length.out = J), # (0,1,0,1,...)
  "Always modified" = rep(1, J)                     # (1,...,1)
)

# Time horizon: the curves are shown until at least 10% of the patients are still at risk,
# to avoid unstable estimates at the end of follow-up (few patients at risk)
at_risk <- function(t) n_distinct(dat_weighted$id[dat_weighted$tstart < t & dat_weighted$tstop >= t])
event_times <- sort(unique(dat_weighted$tstop[dat_weighted$event == 1]))
n_risk <- sapply(event_times, at_risk)
t_max <- max(event_times[n_risk >= 0.10 * n_distinct(dat_weighted$id)])
t_max / 365.25

surv_curves <- bind_rows(lapply(names(strategies), function(s)
  surv_strategy(cox_msm, strategies[[s]], cycle_start, t_max) %>% mutate(strategy = s))) %>%
  mutate(strategy = factor(strategy, levels = names(strategies)),
         years = time / 365.25)

# Marginal survival probabilities every year
surv_curves %>%
  group_by(strategy) %>%
  reframe(year = 1:floor(max(years)),
          surv = sapply(year, function(y) surv[max(which(years <= y))])) %>%
  pivot_wider(names_from = year, names_prefix = "year_", values_from = surv)

save(surv_curves, file = "results/CoxMSM_survival.RData")


############
# Figure 7 #
############
# Left panel: full range of survival probabilities; right panel: zoomed view
panels <- c("Marginal survival probability under different strategies", "(zoomed)")
df_fig7 <- bind_rows(surv_curves %>% mutate(panel = panels[1]),
                     surv_curves %>% mutate(panel = panels[2])) %>%
  mutate(panel = factor(panel, levels = panels))

figure_7 <- ggplot(df_fig7, aes(x = years, y = surv, color = strategy)) +
  geom_step(linewidth = 0.8) +
  # extend the y-axis of the left panel to 0
  geom_blank(data = data.frame(years = 0, surv = 0, strategy = names(strategies)[1],
                               panel = factor(panels[1], levels = panels))) +
  facet_wrap(~ panel, scales = "free_y") +
  scale_color_manual(values = c("Never modified" = "#3A6EA5", "Alternating" = "#F5A300",
                                "Always modified" = "#E3607A"),
                     breaks = c("Alternating", "Always modified", "Never modified"),
                     name = "Strategy") +
  scale_x_continuous(breaks = 0:ceiling(max(surv_curves$years))) +
  labs(x = "Time since treatment initiation (years)", y = "Survival probability") +
  theme_light(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(color = "black", size = 12, hjust = 0))
figure_7
ggsave("results/figures/figure_7.pdf", plot = figure_7, width = 12, height = 5.5)
