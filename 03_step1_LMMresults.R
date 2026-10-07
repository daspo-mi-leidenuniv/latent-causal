##########################################################
# Step 1: Results of the final LMM and RR profiles        #
##########################################################
# Set the working directory to the root of this repository before running the script
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(LMest)

source("functions/lmm_functions.R")
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Final LM model selected in 02_step1_LMMselection.R
load("results/LMM_final.RData")
k <- lmm_final$k
TT <- lmm_final$TT
lmm_label(init_cov, trans_cov)


############
# Figure 5 #
############
# (a) Estimated conditional response probabilities P(Y_r = 1 | U = u)
df_psi <- expand.grid(state = 1:k, Y = Y_VARS) %>%
  mutate(prob = as.vector(lmm_final$Psi[2, , ]),
         state = factor(state),
         Y = factor(Y, levels = rev(Y_VARS)))

p_psi <- ggplot(df_psi, aes(x = state, y = Y, fill = prob)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", prob)), size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1), name = "Probability") +
  labs(x = "Latent state u", y = NULL, title = "(a) Conditional response probabilities P(Y=1|U=u)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

# (b) Latent states prevalence over cycles averaged over all the subjects
df_prev <- as.data.frame(as.table(state_prevalence(lmm_final)), responseName = "prob") %>%
  mutate(cycle = as.numeric(as.character(cycle)))

p_prev <- ggplot(df_prev, aes(x = cycle, y = 100 * prob, color = state)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:TT) +
  ylim(0, 100) +
  labs(x = "Cycle j", y = "Probability (%)", color = "Latent state",
       title = "(b) Latent states prevalence over cycles") +
  theme_light(base_size = 12)

figure_5 <- ggarrange(p_psi, p_prev, ncol = 2)
figure_5
ggsave("results/figures/figure_5.pdf", plot = figure_5, width = 13, height = 5)


###########
# Table 1 #
###########
# (a) Estimated regression parameters (and standard errors) of the initial probabilities
#     (multinomial logits with respect to latent state 1)
Be <- as.matrix(lmm_final$Be)
seBe <- as.matrix(lmm_final$seBe)
tab_initial <- data.frame(part = "Initial probabilities: logit parameters (SE)",
                          row = rownames(Be),
                          matrix(sprintf("%.3f (%.3f)", Be, seBe), nrow = nrow(Be),
                                 dimnames = list(NULL, paste0("state_", 2:k))))

# (b) Estimated transition probabilities under A = 0 and A = 1
tab_trans <- bind_rows(lapply(0:1, function(a){
  Pi <- transition_matrix(lmm_final, lmm_data, a)
  data.frame(part = paste0("Transition probabilities: A = ", a),
             row = paste0("from state ", 1:k),
             matrix(sprintf("%.3f", Pi), nrow = k, dimnames = list(NULL, paste0("state_", 1:k))))
}))

table_1 <- bind_rows(tab_initial, tab_trans) %>% relocate(any_of(paste0("state_", 1:k)), .after = row)
table_1
write.csv(table_1, "results/tables/table_1.csv", row.names = FALSE, na = "")


###########################################
# Longitudinal relative risk (RR) profiles #
###########################################
# Reference latent state R: state with the lowest probabilities of Y = 1
ref_state <- reference_state(lmm_final)
ref_state

ids <- unique(lmm_data$id)
RR_profiles <- rr_profiles(lmm_final, ids, ref = ref_state)
head(RR_profiles)

save(RR_profiles, ref_state, file = "results/RR_profiles.RData")


############
# Figure 6 #
############
# RR profiles of two example patients
ids_example <- c(1, 2)
df_rr <- RR_profiles %>%
  filter(id %in% ids_example) %>%
  pivot_longer(starts_with("RR_"), names_to = "state", values_to = "RR") %>%
  mutate(state = paste0("log(P(U=", sub("RR_", "", state), ")/P(U=", ref_state, "))"),
         id = paste("Patient", id))

figure_6 <- ggplot(df_rr, aes(x = cycle, y = RR, color = state)) +
  geom_hline(yintercept = 0, color = "gray30", linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ id) +
  scale_x_continuous(breaks = 1:TT) +
  labs(x = "Cycle j", y = "Relative risk profile", color = NULL,
       title = "Longitudinal relative risk profiles") +
  theme_light(base_size = 12) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
figure_6
ggsave("results/figures/figure_6.pdf", plot = figure_6, width = 10, height = 5)
