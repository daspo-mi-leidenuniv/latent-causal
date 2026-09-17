library(dplyr)
library(ggplot2)
library(tidyr)

setwd("~/github/lmm-ipw")

long_dat = read.csv2("data/fake_data.csv")
head(long_dat)


############
# Figure 3 #
############
df_cumA <- long_dat %>%
  group_by(id) %>%
  summarise(cumA = sum(A == 1, na.rm = TRUE), .groups = "drop")

figure_3 = ggplot(df_cumA, aes(x = factor(cumA, levels = 0:5))) +
  geom_bar(fill = "#1C3F94") +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Number of cycles with A=1", y = "Number of patients",
       title = "Cumulative treatment history") +
  theme_light()
figure_3

ggsave("results/figures/figure_3.pdf", plot = figure_3, width = 5, height = 6)


############
# Figure 4 #
############

# 1. Patient x cycle grid
J <- max(long_dat$cycle) 
Y_vars <- paste0("Y", 1:8)

dat_plot <- expand_grid(id = unique(long_dat$id), cycle = 1:J) %>%
  left_join(long_dat %>% 
              select(id, cycle, all_of(Y_vars)), 
            by = c("id", "cycle")) %>%
  left_join(long_dat %>%
              group_by(id) %>%
              summarise(M = max(cycle),.groups = "drop"),
            by = "id")

# 2. Long format
dat_plot_long <- dat_plot %>%
  pivot_longer(cols = all_of(Y_vars), names_to = "Outcome", 
               values_to = "Value") %>%
  mutate(Value = case_when(cycle > M ~ "Dropout",
                           Value == 0 ~ "0",
                           Value == 1 ~ "1",
                           TRUE ~ NA_character_)) %>%
  filter(!is.na(Value)) %>%
  mutate(Value = factor(Value, levels = c("0", "1", "Dropout")))

# 3. Percentages
dat_percent <- dat_plot_long %>%
  count(Outcome, cycle, Value) %>%
  group_by(Outcome, cycle) %>%
  mutate(Percentage = 100 * n / sum(n)) %>%
  ungroup()

# 4. Plot
figure_4 = ggplot(dat_percent, aes(x = factor(cycle), y = Percentage, fill = Value)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  facet_wrap(~ Outcome, ncol = 4) +
  scale_fill_manual(name = "Observed level",
                    values = c("0" = "#fce6c6", "1" = "#9e2a2b","Dropout" = "blue")) +
  labs(x = "Cycle j", y = "Percentage of patients (%)",
       title = "Observed responses over cycles") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 11), legend.position = "top")
figure_4

ggsave("results/figures/figure_4.pdf", plot = figure_4, width = 15, height = 8)

