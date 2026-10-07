##########################################################
# Step 2: Specification of IPTW and IPDW                 #
##########################################################
# Set the working directory to the root of this repository before running the script
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(ipw)

source("functions/lmm_functions.R")
source("functions/ipw_functions.R")
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

long_dat <- read.csv2("data/fake_data.csv")

# RR profiles computed in 03_step1_LMMresults.R
load("results/RR_profiles.RData")
rr_vars <- grep("^RR_", names(RR_profiles), value = TRUE)

# Datasets for the treatment (dat_A) and dropout (dat_D) models
step2_data <- build_step2_data(long_dat, RR_profiles)

# Candidate specifications (same covariate sets for IPTW and IPDW):
#   IPW0: RR profiles
#   IPW1: baseline covariates + RR profiles
#   IPW2: RR profiles + previous treatment A_{j-1}                        (numerator: A_{j-1})
#   IPW3: baseline covariates + RR profiles + previous treatment A_{j-1}  (numerator: A_{j-1})
#   IPW4: RR profiles + cumulative treatment up to j-1                    (numerator: cumulative treatment)
#   IPW5: baseline covariates + RR profiles + cumulative treatment        (numerator: cumulative treatment)
specs <- ipw_specs(rr_vars)


# Stabilized weights under each specification
weights_all <- bind_rows(lapply(names(specs), function(s){
  step2_data$dat_A %>%
    arrange(id, cycle) %>%
    select(id, cycle) %>%
    mutate(spec = s,
           sw_iptw = compute_iptw(step2_data$dat_A, specs[[s]]),
           sw_ipdw = compute_ipdw(step2_data$dat_A, step2_data$dat_D, specs[[s]]))
}))

# Distribution of the stabilized weights by cycle
weights_all %>%
  pivot_longer(c(sw_iptw, sw_ipdw), names_to = "weight", values_to = "sw") %>%
  group_by(weight, spec, cycle) %>%
  summarise(mean = mean(sw), min = min(sw), max = max(sw), .groups = "drop") %>%
  print(n = Inf)

plot_weights <- function(dat, var, title){
  ggplot(dat, aes(x = factor(cycle), y = log(.data[[var]]))) +
    geom_hline(yintercept = 0, color = "gray30", linetype = "dashed") +
    geom_boxplot(fill = "#1C3F94", alpha = 0.3, outlier.size = 0.8) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 2.5, color = "#9e2a2b") +
    facet_wrap(~ spec, ncol = 2) +
    labs(x = "Cycle j", y = "Log stabilized weight", title = title) +
    theme_light(base_size = 12) +
    theme(strip.text = element_text(face = "bold"))
}


###########################
# Supplementary Figure C.1 #
###########################
figure_C1 <- plot_weights(weights_all %>% mutate(spec = sub("IPW", "IPTW", spec)), "sw_iptw",
                          "Stabilized IPTW")
figure_C1
ggsave("results/figures/figure_C1.pdf", plot = figure_C1, width = 8, height = 9)


###########################
# Supplementary Figure C.2 #
###########################
figure_C2 <- plot_weights(weights_all %>% mutate(spec = sub("IPW", "IPDW", spec)) %>% filter(cycle > 1),
                          "sw_ipdw", "Stabilized IPDW")
figure_C2
ggsave("results/figures/figure_C2.pdf", plot = figure_C2, width = 8, height = 9)


##########################
# Final IPTW and IPDW     #
##########################
spec_iptw <- specs[["IPW4"]]
spec_ipdw <- specs[["IPW4"]]

dat_weighted <- add_weights(step2_data, spec_iptw, spec_ipdw)
summary(dat_weighted$sw)

save(step2_data, dat_weighted, spec_iptw, spec_ipdw, file = "results/IPW_final.RData")
