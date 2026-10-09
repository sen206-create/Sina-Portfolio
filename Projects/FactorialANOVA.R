library(tidyverse)
library(readxl)
library(easystats)

#Setting up data
Heart_rate_data = read_excel("Heart rate data.xlsx")

Heart_rate_long <- Heart_rate_data |>
  pivot_longer(cols = c("DBA", "C57BL/6", "BALB/c", "CBA", "SWR"), names_to = "strain", values_to = "heart_rate")
Heart_rate_long <- Heart_rate_long |>
  mutate(strain = as_factor(strain), sex = as_factor(sex))

Heart_rate_long |>
  group_by(strain, sex) |>
  describe_distribution(ci = 0.95, iterations = 1000) |>
  ungroup() |>
  display()

Heart_rate_long |>
  filter(is.na(heart_rate))

ggplot(Heart_rate_long, aes(x = strain, y = heart_rate, fill = sex, colour = sex)) +
  geom_point(size = 1, position = position_jitterdodge(jitter.width = 0.3)) +
  geom_boxplot(alpha = .3, colour = "black", outliers = FALSE) +
  labs(x = 'Mice Strain', y = 'Heart Rate (BPM)' )+
  theme_bw(base_size = 14)

#Levene test for homogeneity of variance
library(car)

leveneTest(Heart_rate_long$heart_rate, interaction(Heart_rate_long$strain, Heart_rate_long$sex), center = median)

#Fit two-way independent factorial ANOVA
Heart_rate_model <- lm(heart_rate ~ strain*sex, data = Heart_rate_long, na.action = "na.exclude")
anova(Heart_rate_model) |>
  model_parameters(es = "omega", ci = 0.95) |>
  display(digits = 3)

#Standardising residuals
Heart_rate_long <- Heart_rate_long |>
  mutate(st_residual = rstandard(Heart_rate_model))

#Checking normality of standardised residuals
ggplot(Heart_rate_long, aes(x = st_residual)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, colour = "black", fill = "magenta") +
  stat_function(fun = dnorm, args = list(
    mean = mean(Heart_rate_long$st_residual, na.rm = TRUE),
    sd = sd(Heart_rate_long$st_residual, na.rm = TRUE))) +
  labs(x ="Standardised Residuals", y = "Density") +
  theme_bw()

library(qqplotr)
ggplot(Heart_rate_long, aes(sample = st_residual)) +
  stat_qq_band(fill = "orange", alpha = 0.2) +
  stat_qq_line() +
  stat_qq_point() +
  labs(x= "theoretical quantiles", y = "standard residual") +
  theme_bw(base_size = 12)
shapiro.test(Heart_rate_long$st_residual)

#Checkin the distance of standardised residuals 
cumulative_percent <- function(var, cut_off = 1.96){
  ecdf_proportions <- abs(var) |>
    ecdf()
  100*(1 - ecdf_proportions(cut_off))
}

Heart_rate_long |>
  summarize(
    z_1.96 = cumulative_percent(st_residual),
    z_2.58 = cumulative_percent(st_residual, cut_off = 2.58),
    z_3.29 = cumulative_percent(st_residual,cut_off = 3.29)
  )

# Plotting interaction plot
ggplot(Heart_rate_long, aes(x = strain, y = heart_rate, colour = sex, group = sex)) +
  labs(title = 'Mean Resting Heart Rate by Mouse Strain and Sex', x = 'Mice Strain', y = 'Heart Rate (BPM)') +
  stat_summary(fun = mean, geom = "point", size = 2, position = position_dodge(0.1), na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line", size = 1, position = position_dodge(0.1), na.rm = TRUE) +
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", width = 0.1, size = 1, position = position_dodge(0.1), na.rm = TRUE) +
  coord_cartesian(ylim = c(0, 1000)) +
    theme_bw(base_size = 14)

#Simple effects analysis
estimate_contrasts(model = Heart_rate_model, contrast = "sex", by = "strain",p_adjust = "sidak") |>
  display()
estimate_contrasts(model = Heart_rate_model, contrast = "strain", by = "sex",p_adjust = "sidak") |>
  display()

#Checking for influential figures
Heart_rate_influence <- influence.measures(Heart_rate_model)$infmat |>
  as_tibble() |>
  rowid_to_column(var = "case_no")

Heart_rate_influence |>
  filter(
    if_any(starts_with("df"), ~ abs(.) > 1)
  )

#Rank based robust ANOVA
library(Rfit)
Heart_rate_rfit_model <- raov(heart_rate ~ strain:sex, data = Heart_rate_long)
Heart_rate_rfit_model
