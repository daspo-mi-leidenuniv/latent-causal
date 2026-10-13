##########################################################
# Step 1: Results of the final LMM and RR profiles        #
##########################################################
# Set the working directory to the root of this repository before running the script
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
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
# Estimated conditional response probabilities P(Y_r = 1 | U = u)
df_psi <- expand.grid(state = 1:k, Y = Y_VARS) %>%
  mutate(prob = as.vector(lmm_final$Psi[2, , ]),
         state = factor(state),
         Y = factor(Y, levels = rev(Y_VARS)))

figure_5 <- ggplot(df_psi, aes(x = state, y = Y, fill = prob)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", prob)), size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1), name = "Probability") +
  labs(x = "Latent state u", y = NULL, title = "Conditional response probabilities P(high-grade toxicity | U = u)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

figure_5
ggsave("results/figures/figure_5.pdf", plot = figure_5, width = 7, height = 5)


###########################
# Supplementary Figure B.1 #
###########################
# Latent states prevalence over cycles averaged over all the subjects
df_prev <- as.data.frame(as.table(state_prevalence(lmm_final)), responseName = "prob") %>%
  mutate(cycle = as.numeric(as.character(cycle)))

# Colors of the latent states 1-5
state_colors <- c("1" = "#2A9D8F", "2" = "#F4A261", "3" = "#595959", "4" = "#E63946", "5" = "#457B9D")

figure_B1 <- ggplot(df_prev, aes(x = cycle, y = 100 * prob, color = state)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = state_colors, name = "Latent state (u)") +
  scale_x_continuous(breaks = 1:TT) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25)) +
  labs(x = "Cycle (j)", y = "Prevalence (%)",
       title = "Latent states prevalence over cycles") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 11),
        legend.key.height = unit(1.2, "cm"),
        legend.key.width = unit(1, "cm"))

figure_B1
ggsave("results/figures/figure_B1.pdf", plot = figure_B1, width = 7, height = 5)


###########
# Table 1 #
###########
# Estimated regression parameters affecting the distribution of the initial probabilities
# (multinomial logits with respect to latent state 1, whose parameters are fixed to 0)
Be <- cbind(0, as.matrix(lmm_final$Be))
dimnames(Be) <- list(rownames(lmm_final$Be), 1:k)
Be
# Estimated transition probabilities from state u~ (rows) to state u (columns) when A = 0 and A = 1
Pi <- lapply(c("0" = 0, "1" = 1), function(a) transition_matrix(lmm_final, lmm_data, a))
Pi

# Table in csv format
table_1 <- bind_rows(
  data.frame(section = "Regression parameters for initial probabilities",
             term = term_labels(rownames(Be)),
             symbol = paste0("beta_", 0:(nrow(Be) - 1), "u"),
             matrix(ifelse(Be == 0, "0", sprintf("%.2f", Be)), nrow = nrow(Be), dimnames = list(NULL, paste0("u=", 1:k))),
             check.names = FALSE),
  bind_rows(lapply(names(Pi), function(a)
    data.frame(section = paste0("Transition probabilities from u~ to u when A_ij = ", a),
               term = "",
               symbol = paste0("u~=", 1:k),
               matrix(sprintf("%.2f", Pi[[a]]), nrow = k, dimnames = list(NULL, paste0("u=", 1:k))),
               check.names = FALSE)))
)
table_1
write.csv(table_1, "results/tables/table_1.csv", row.names = FALSE)

# Table in LaTeX format (requires the booktabs package)
writeLines(latex_table_1(Be, Pi), "results/tables/table_1.tex")


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
