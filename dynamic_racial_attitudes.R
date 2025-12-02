#Lauren Palladino Olson

library(dplyr)
library(tidyr)
library(broom)
library(survival)
library(ggplot2)
library(sandwich)
library(lmtest)
library(readr)
library(moments)
library(scales)

setwd("~/Documents/[RESEARCH] Dynamic Racial Attitudes")
outdir <- "replication_outputs"
if(!dir.exists(outdir)) dir.create(outdir)

anes <- read_csv("anes_mergedfile_2016-2020-2024panel_20251030.csv", na = c("", " ", "-9","-8","-7","-2","-1"))

map_vote <- function(x, dem_code, rep_code){
  xnum <- as.numeric(x)
  v <- rep(NA_real_, length(xnum))
  v[xnum == rep_code] <- 1
  v[xnum == dem_code] <- 0
  v
}
anes <- anes %>%
  mutate(
    rep2016 = map_vote(V162034a, dem_code = 1, rep_code = 2),
    rep2020 = map_vote(V202073, dem_code = 1, rep_code = 2),
    rep2024 = map_vote(V242067, dem_code = 1, rep_code = 2)
  )

# figure 2 - vote switching counts by category
anes_switch <- anes %>%
  mutate(
    switch_2016_2020 = case_when(
      rep2016 == 0 & rep2020 == 0 ~ "Stayed Dem",
      rep2016 == 1 & rep2020 == 1 ~ "Stayed Rep",
      rep2016 == 0 & rep2020 == 1 ~ "Dem→Rep",
      rep2016 == 1 & rep2020 == 0 ~ "Rep→Dem",
      TRUE ~ NA_character_
    ),
    switch_2020_2024 = case_when(
      rep2020 == 0 & rep2024 == 0 ~ "Stayed Dem",
      rep2020 == 1 & rep2024 == 1 ~ "Stayed Rep",
      rep2020 == 0 & rep2024 == 1 ~ "Dem→Rep",
      rep2020 == 1 & rep2024 == 0 ~ "Rep→Dem",
      TRUE ~ NA_character_
    )
  ) %>%
  pivot_longer(
    cols = c(switch_2016_2020, switch_2020_2024),
    names_to = "period",
    values_to = "vote_switch"
  ) %>%
  mutate(
    period = case_when(
      period == "switch_2016_2020" ~ "2016–2020",
      period == "switch_2020_2024" ~ "2020–2024"
    ),
    vote_switch = factor(
      vote_switch,
      levels = c("Stayed Dem", "Stayed Rep", "Dem→Rep", "Rep→Dem")
    )
  ) %>%
  filter(!is.na(vote_switch))

counts <- anes_switch %>%
  group_by(period, vote_switch) %>%
  summarise(N = n(), .groups = "drop")

ggplot(counts, aes(x = vote_switch, y = N, fill = vote_switch)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = N), vjust = -0.3, size = 4) +
  facet_wrap(~period) +
  scale_fill_manual(values = c(
    "Stayed Dem" = "#1f78b4",
    "Stayed Rep" = "#e31a1c",
    "Dem→Rep"   = "#fb9a99",
    "Rep→Dem"   = "#a6cee3"
  )) +
  labs(
    x = "Vote Switch Category",
    y = "Number of Respondents",
    fill = "Vote Switch"
  ) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

# clean thermometers and discrimination
clean_var <- function(x) { xnum <- as.numeric(x); xnum[xnum < 0] <- NA; xnum }
anes <- anes %>%
  mutate_at(vars(V162314,V162312,V162357,V202482,V202480,V202527,V242518,V242516,V242549),
            clean_var)

# consistent var names
panel <- anes %>%
  transmute(id = V160001_orig,
            rep2016, rep2020, rep2024,
            white2016 = V162314, black2016 = V162312, disc2016 = V162357,
            white2020 = V202482, black2020 = V202480, disc2020 = V202527,
            white2024 = V242518, black2024 = V242516, disc2024 = V242549) %>%
  pivot_longer(cols = starts_with(c("rep","white","black","disc")),
               names_to = c(".value","wave"),
               names_pattern = "([a-z]+)([0-9]+)") %>%
  mutate(year = as.integer(wave),
         wmb = white - black) %>%
  select(id, year, rep = rep, disc, wmb)

# Within-person transformation for FE LPM
panel_fe <- panel %>%
  group_by(id) %>%
  mutate(rep_mean = mean(rep, na.rm=TRUE),
         rep_within = rep - rep_mean,
         disc_mean = mean(disc, na.rm=TRUE),
         disc_within = disc - disc_mean,
         wmb_mean = mean(wmb, na.rm=TRUE),
         wmb_within = wmb - wmb_mean,
         yr2020 = as.numeric(year==2020),
         yr2024 = as.numeric(year==2024),
         yr2020_mean = mean(yr2020, na.rm=TRUE),
         yr2024_mean = mean(yr2024, na.rm=TRUE),
         yr2020_within = yr2020 - yr2020_mean,
         yr2024_within = yr2024 - yr2024_mean) %>%
  ungroup()

fe_data <- panel_fe %>% filter(!is.na(rep_within), !is.na(disc_within), !is.na(wmb_within))

# Table 1: FE LPM
fe_lpm <- lm(rep_within ~ disc_within + wmb_within + yr2020_within + yr2024_within, data = fe_data)
fe_lpm_se <- coeftest(fe_lpm, vcov = vcovCL(fe_lpm, cluster = fe_data$id))
write.csv(tidy(fe_lpm, conf.int=TRUE), file.path(outdir,"fe_lpm_tidy.csv"), row.names=FALSE)
capture.output(summary(fe_lpm), file=file.path(outdir,"fe_lpm_summary.txt"))
capture.output(fe_lpm_se, file=file.path(outdir,"fe_lpm_cluster_se.txt"))

# Table 1: Approx FE-logit
switcher_ids <- panel %>%
  group_by(id) %>%
  summarize(n_nonmiss = sum(!is.na(rep)), n_unique = n_distinct(rep[!is.na(rep)])) %>%
  filter(n_nonmiss >= 2, n_unique > 1) %>%
  pull(id)

panel_sw <- panel_fe %>% filter(id %in% switcher_ids) %>% filter(!is.na(rep), !is.na(disc_within), !is.na(wmb_within))

fe_logit <- glm(rep ~ disc_within + wmb_within + yr2020_within + yr2024_within + factor(id),
                data = panel_sw, family = binomial(link="logit"))
fe_logit_se <- coeftest(fe_logit, vcov = vcovCL(fe_logit, cluster = panel_sw$id))
write.csv(tidy(fe_logit, conf.int=TRUE), file.path(outdir,"fe_logit_tidy.csv"), row.names=FALSE)
capture.output(summary(fe_logit), file=file.path(outdir,"fe_logit_summary.txt"))
capture.output(fe_logit_se, file=file.path(outdir,"fe_logit_cluster_se.txt"))

#Figure 3

panel <- panel %>%
  mutate(
    discrim_01 = scales::rescale(disc, to = c(0,1)),
    wmb_01     = scales::rescale(wmb, to = c(0,1))
  )

panel <- panel %>%
  arrange(id, year) %>%
  group_by(id) %>%
  mutate(
    delta_discrim_01 = discrim_01 - lag(discrim_01),
    delta_wmb_01     = wmb_01 - lag(wmb_01)
  ) %>%
  ungroup() %>%
  filter(!is.na(delta_discrim_01) & !is.na(delta_wmb_01))

delta_seq <- seq(-1, 1, length.out = 200)  

coef_discrim <- 0.01712
se_discrim   <- 0.00489

coef_wmb     <- 0.00031
se_wmb       <- 0.00024

pred_df <- bind_rows(
  data.frame(
    delta = delta_seq,
    predictor = "Perceived Discrimination",
    mean_rep = coef_discrim * delta_seq,
    lower    = (coef_discrim - 1.96*se_discrim) * delta_seq,
    upper    = (coef_discrim + 1.96*se_discrim) * delta_seq
  ),
  data.frame(
    delta = delta_seq,
    predictor = "Relative Warmth",
    mean_rep = coef_wmb * delta_seq,
    lower    = (coef_wmb - 1.96*se_wmb) * delta_seq,
    upper    = (coef_wmb + 1.96*se_wmb) * delta_seq
  )
)

p <- ggplot(pred_df, aes(x = delta, y = mean_rep, group = predictor)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = predictor), alpha = 0.25) +
  geom_line(aes(color = predictor, linetype = predictor), size = 0.9) +
  scale_color_manual(values = c(
    "Perceived Discrimination" = "black",
    "Relative Warmth" = "black"
  )) +
  scale_linetype_manual(values = c(
    "Perceived Discrimination" = "solid",
    "Relative Warmth" = "dotted"
  )) +
  
  labs(
    x = "Racial Attitude Change (Rescaled)",
    y = "Predicted Probability of Voting Republican",
    color = "Predictor",
    fill = "Predictor",
    linetype = "Predictor"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "bottom"
  )

ggsave(file.path(outdir, "combined_rescaled_within_person_merged_legend.png"),
       p, width = 8, height = 5, dpi = 300)

# distributions of racial attitudes across wave
panel_changes <- panel %>%
  arrange(id, year) %>%
  group_by(id) %>%
  mutate(
    delta_disc = disc - lag(disc),
    delta_wmb  = wmb  - lag(wmb),
    prior_year = lag(year),
    period = paste0(prior_year, "-", year)
  ) %>%
  ungroup() %>%
  # keep only the two intervals of interest
  filter(!is.na(period) & period %in% c("2016-2020", "2020-2024")) %>%
  select(id, period, delta_disc, delta_wmb)

df_long <- panel_changes %>%
  pivot_longer(
    cols = c(delta_disc, delta_wmb),
    names_to = "measure",
    values_to = "change_score"
  ) %>%
  filter(!is.na(change_score)) %>%
  mutate(
    measure = as.character(measure),
    measure = dplyr::recode(
      measure,
      "delta_disc" = "Perceived Discrimination",
      "delta_wmb"  = "WMB Warmth"
    )
  )

df_long <- df_long %>%
  group_by(measure) %>%
  mutate(change_rescaled = scales::rescale(change_score, to = c(-1, 1))) %>%
  ungroup()

violin_plot <- ggplot(df_long, aes(x = measure, y = change_rescaled, fill = measure)) +
  geom_violin(trim = FALSE, alpha = 0.45) +
  facet_wrap(~ period, ncol = 2) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(
    x = "",
    y = "Change (rescaled −1 to 1)",
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(outdir, "appendix_change_violin_boxplots_periods.png"),
       plot = violin_plot, width = 8, height = 4.5, dpi = 300)

# Table 2: First-difference and reverse causality models
wide <- panel %>% pivot_wider(id_cols = id, names_from = year, values_from = c(rep,wmb,disc), names_sep="_") %>%
  mutate(d_rep_16_20 = rep_2020 - rep_2016,
         d_rep_20_24 = rep_2024 - rep_2020,
         d_disc_16_20 = disc_2020 - disc_2016,
         d_disc_20_24 = disc_2024 - disc_2020,
         d_wmb_16_20 = wmb_2020 - wmb_2016,
         d_wmb_20_24 = wmb_2024 - wmb_2020)

diffs <- bind_rows(
  wide %>% transmute(id, d_rep = d_rep_16_20, d_disc = d_disc_16_20, d_wmb = d_wmb_16_20),
  wide %>% transmute(id, d_rep = d_rep_20_24, d_disc = d_disc_20_24, d_wmb = d_wmb_20_24)
) %>% filter(!is.na(d_rep) & !is.na(d_disc))

diff_model <- lm(d_rep ~ d_disc , data = diffs)
capture.output(summary(diff_model), file=file.path(outdir,"diff_model_summary.txt"))

rev_model <- lm(d_disc ~ d_rep , data = diffs)
capture.output(summary(rev_model), file=file.path(outdir,"rev_model_summary.txt"))

# Table 3: Interaction with ideo strength
if ("V161126" %in% names(anes)) {
  pidbase <- anes %>% 
    filter(V161126 <= 7) %>%              # Drop respondents with V161126 > 7
    select(id = V160001_orig, pid = V161126) %>% 
    distinct()
  # Merge into panel FE data
  panel_fe2 <- panel_fe %>% 
    left_join(pidbase, by = "id") %>% 
    filter(!is.na(pid)) %>%             
    mutate(pid_strength = abs(pid - median(pid, na.rm = TRUE)))
  interact_data <- panel_fe2 %>% 
    filter(!is.na(disc_within),
           !is.na(pid_strength),
           !is.na(rep_within))
  
  int_model <- lm(rep_within ~ disc_within * pid_strength + 
                    wmb_within + yr2020_within + yr2024_within, 
                  data = interact_data)
  
  int_se <- coeftest(int_model, vcov = vcovCL(int_model, cluster = interact_data$id))
  
  write.csv(tidy(int_model, conf.int = TRUE),
            file.path(outdir, "interaction_tidy.csv"),
            row.names = FALSE)
  
  capture.output(summary(int_model),
                 file = file.path(outdir, "interaction_summary.txt"))
  
  capture.output(int_se,
                 file = file.path(outdir, "interaction_cluster_se.txt"))
}

# marginal effects across 7 point ideological range
marginal_effects <- ggeffect(int_model, terms = c("disc_within", "pid_strength [all]"))

plot_df <- as.data.frame(marginal_effects)

ggplot(plot_df, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high)) +
  geom_line(size = 1) +
  geom_ribbon(alpha = 0.2) +
  labs(
    x = "Ideological Strength (Absolute Distance from Median Partisan ID)",
    y = "Predicted Within-Person Change in Republican Vote",
    title = "Marginal Effect of Changes in Perceived Discrimination by Ideological Strength"
  ) +
  theme_bw(base_size = 14)

#Appendix A: Diagnostics and Var Properties

#Rescaling
rescale_minus1_1 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  2 * (x - rng[1]) / (rng[2] - rng[1]) - 1
}

panel_fe_rescaled <- panel_fe %>%
  mutate(
    disc_within_r = rescale_minus1_1(disc_within),
    wmb_within_r = rescale_minus1_1(wmb_within)
  )

#Descriptive stats
desc_stats <- panel_fe_rescaled %>%
  summarise(
    mean_disc = mean(disc_within_r, na.rm=TRUE),
    sd_disc = sd(disc_within_r, na.rm=TRUE),
    min_disc = min(disc_within_r, na.rm=TRUE),
    max_disc = max(disc_within_r, na.rm=TRUE),
    skew_disc = skewness(disc_within_r, na.rm=TRUE),
    mean_wmb = mean(wmb_within_r, na.rm=TRUE),
    sd_wmb = sd(wmb_within_r, na.rm=TRUE),
    min_wmb = min(wmb_within_r, na.rm=TRUE),
    max_wmb = max(wmb_within_r, na.rm=TRUE),
    skew_wmb = skewness(wmb_within_r, na.rm=TRUE)
  )

write.csv(desc_stats, file.path(outdir, "appendix_measure_desc_stats_minus1_1.csv"), row.names=FALSE)

#Rescaling
panel_long <- panel_fe_rescaled %>%
  select(id, disc_within_r, wmb_within_r) %>%
  pivot_longer(cols = c(disc_within_r, wmb_within_r),
               names_to = "measure",
               values_to = "delta") %>%
  filter(!is.na(delta))

#Density plot
density_plot <- ggplot(panel_long, aes(x = delta, fill = measure)) +
  geom_density(alpha = 0.4) +
  theme_bw() +
  labs(x = "Within-person change (rescaled -1 to 1)", 
       y = "Density",
       title = "Distribution of within-person change in racial attitudes")

ggsave(filename = file.path(outdir, "appendix_density_withinperson_change_minus1_1.png"),
       plot = density_plot,
       width = 7, height = 5)

#Appendix B: Lead models 
panel_lead <- panel_fe %>%
  arrange(id, year) %>%
  group_by(id) %>%
  mutate(
    disc_within_lead = lead(disc_within),
    wmb_within_lead  = lead(wmb_within),
    rep_within_lag   = rep_within  # vote at time t predicts next attitude
  ) %>%
  ungroup() %>%
  filter(!is.na(disc_within_lead) & !is.na(rep_within_lag))  # remove last obs per person

rescale_minus1_1 <- function(x) {
  2 * ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))) - 1
}

panel_lead <- panel_lead %>%
  group_by(id) %>%
  mutate(
    disc_within_lead_rescaled = rescale_minus1_1(disc_within_lead),
    wmb_within_lead_rescaled  = rescale_minus1_1(wmb_within_lead)
  ) %>%
  ungroup()

lead_disc_model <- lm(disc_within_lead_rescaled ~ rep_within_lag, data = panel_lead)
lead_wmb_model  <- lm(wmb_within_lead_rescaled ~ rep_within_lag, data = panel_lead)

summary(lead_disc_model)
summary(lead_wmb_model)
