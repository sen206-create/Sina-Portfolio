library(tidyverse)
library(readxl)
library(easystats)

#Setting up data
bloodvessel_ACh = read_excel('bloodvessel ACh.xlsx')

bloodvessel_ACh <- bloodvessel_ACh |>
  mutate(rat_ID = factor(row_number()))

bloodvessel_ACh_long = bloodvessel_ACh |>
  pivot_longer(cols = c('Control diameter (microns)','toxin(in)', 'toxin(out)', 'ACh', 'toxin(in) + ACh', 'toxin(out) + ACh'), names_to = 'Treatment', values_to = 'Diameter')
bloodvessel_ACh_long <- bloodvessel_ACh_long |>
  mutate(Treatment = as_factor(Treatment), sex = as_factor(sex))

bloodvessel_ACh_long |>
  group_by(Treatment, sex) |>
  describe_distribution(ci = 0.95, iterations = 1000) |>
  ungroup() |>
  display()

ggplot(bloodvessel_ACh_long, aes(x = Treatment, y = Diameter, colour = sex, group = sex)) +
  stat_summary(fun = "mean", geom ="point", size = 2, position = position_dodge(0.1)) +
  stat_summary(fun = "mean", geom ="line", size = 1, position = position_dodge(0.1)) +
  stat_summary(fun.data = "mean_cl_normal", geom ="errorbar", width = 0.1, position = position_dodge(0.1)) +
  theme_bw(base_size = 14)+
  labs(title = 'Blood Vessel Diameter by Treatment and Sex', x = 'Treatment', y = "Blood Vessel Diamater (microns)")
  
control_diameter = bloodvessel_ACh_long |>
  filter(Treatment == 'Control diameter (microns)')

#Fit a repeated measures ANOVA
library(afex)
control_diameter_afx = aov_4(Diameter ~ Treatment + (Treatment|rat_ID), data = bloodvessel_ACh_long)
model_parameters(control_diameter_afx, es_type = "omega")|>
  display(use_symbols = TRUE)

#Standardising residuals
bloodvessel_ACh_long <- bloodvessel_ACh_long |>
  mutate(st_residual = (residuals(control_diameter_afx) - mean(residuals(control_diameter_afx)))/sd(residuals(control_diameter_afx)))

#Checking normality of standardised residuals
ggplot(bloodvessel_ACh_long, aes(x = st_residual)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, colour = "black", fill = "magenta") +
  stat_function(fun = dnorm, args = list(
    mean = mean(bloodvessel_ACh_long$st_residual, na.rm = TRUE),
    sd = sd(bloodvessel_ACh_long$st_residual, na.rm = TRUE))) +
  labs(x ="Standardised Residuals", y = "Density") +
  theme_bw()

#Q-Q plot
library(qqplotr)
ggplot(bloodvessel_ACh_long, aes(sample = st_residual)) +
  stat_qq_band(fill = "orange", alpha = 0.2) +
  stat_qq_line() +
  stat_qq_point() +
  labs(x= "theoretical quantiles", y = "standard residual") +
  theme_bw(base_size = 12)
shapiro.test(bloodvessel_ACh_long$st_residual)

#Checking the distance of standardised residuals 
cumulative_percent <- function(var, cut_off = 1.96){
  ecdf_proportions <- abs(var) |>
    ecdf()
  100*(1 - ecdf_proportions(cut_off))
}

bloodvessel_ACh_long |>
  summarize(
    z_1.96 = cumulative_percent(st_residual),
    z_2.58 = cumulative_percent(st_residual, cut_off = 2.58),
    z_3.29 = cumulative_percent(st_residual,cut_off = 3.29)
  )

#Non-normality of standardised residuals, will have to transform them to resolve.

#Square root transformation
bloodvessel_ACh_sqrt <- bloodvessel_ACh_long |>
  mutate(sqrt_Diameter = sqrt(Diameter))

#Fit repeated measures ANOVA on transformed data
bloodvessel_ACh_sqrt_afx = aov_4(sqrt_Diameter ~ Treatment + (Treatment|rat_ID), data = bloodvessel_ACh_sqrt)

model_parameters(bloodvessel_ACh_sqrt_afx, es_type = "omega") |>
  display(use_symbols = TRUE)

#Standardising new residuals
bloodvessel_ACh_sqrt <- bloodvessel_ACh_sqrt |>
  mutate(st_residual = (residuals(bloodvessel_ACh_sqrt_afx) - mean(residuals(bloodvessel_ACh_sqrt_afx)))/sd(residuals(bloodvessel_ACh_sqrt_afx)))

#Checking normality of standardised residuals
ggplot(bloodvessel_ACh_sqrt, aes(x = st_residual)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, colour = "black", fill = "magenta") +
  stat_function(fun = dnorm, args = list(
    mean = mean(bloodvessel_ACh_sqrt$st_residual, na.rm = TRUE),
    sd = sd(bloodvessel_ACh_sqrt$st_residual, na.rm = TRUE))) +
  labs(x = "Standardised Residuals", y = "Density") +
  theme_bw()

#Q-Q plot
ggplot(bloodvessel_ACh_sqrt, aes(sample = st_residual)) +
  stat_qq_band(fill = "orange", alpha = 0.2) +
  stat_qq_line() +
  stat_qq_point() +
  labs(x = "Theoretical Quantiles", y = "Standardised Residuals") +
  theme_bw(base_size = 12)

shapiro.test(bloodvessel_ACh_sqrt$st_residual)

#Checking the distance of standardised residuals
bloodvessel_ACh_sqrt |>
  summarize(
    z_1.96 = cumulative_percent(st_residual),
    z_2.58 = cumulative_percent(st_residual, cut_off = 2.58),
    z_3.29 = cumulative_percent(st_residual, cut_off = 3.29)
  )

#Sqrt transformation did not resolve non-normality issue

#Log transformation
bloodvessel_ACh_log <- bloodvessel_ACh_long |>
  mutate(log_Diameter = log(Diameter))

#Fit repeated measures ANOVA on log transformed data
bloodvessel_ACh_log_afx = aov_4(log_Diameter ~ Treatment + (Treatment|rat_ID), data = bloodvessel_ACh_log)

model_parameters(bloodvessel_ACh_log_afx, es_type = "omega") |>
  display(use_symbols = TRUE)

#Standardising residuals
bloodvessel_ACh_log <- bloodvessel_ACh_log |>
  mutate(st_residual = (residuals(bloodvessel_ACh_log_afx) - mean(residuals(bloodvessel_ACh_log_afx)))/sd(residuals(bloodvessel_ACh_log_afx)))

#Q-Q plot
ggplot(bloodvessel_ACh_log, aes(sample = st_residual)) +
  stat_qq_band(fill = "orange", alpha = 0.2) +
  stat_qq_line() +
  stat_qq_point() +
  labs(x = "Theoretical Quantiles", y = "Standardised Residuals") +
  theme_bw(base_size = 12)

shapiro.test(bloodvessel_ACh_log$st_residual)

#Fit a mixed two-way ANOVA
bloodvessel_ACh_mixed_afx = aov_4(log_Diameter ~ Treatment*sex + (Treatment|rat_ID), data = bloodvessel_ACh_log)

model_parameters(bloodvessel_ACh_mixed_afx, es_type = "omega") |>
  display(use_symbols = TRUE)

#Standardising residuals
bloodvessel_ACh_log <- bloodvessel_ACh_log |>
  mutate(st_residual = (residuals(bloodvessel_ACh_mixed_afx) - mean(residuals(bloodvessel_ACh_mixed_afx)))/sd(residuals(bloodvessel_ACh_mixed_afx)))

#Q-Q plot
ggplot(bloodvessel_ACh_log, aes(sample = st_residual)) +
  stat_qq_band(fill = "orange", alpha = 0.2) +
  stat_qq_line() +
  stat_qq_point() +
  labs(x = "Theoretical Quantiles", y = "Standardised Residuals") +
  theme_bw(base_size = 12)

shapiro.test(bloodvessel_ACh_log$st_residual)

#Checking homogeneity of variance
check_homogeneity(bloodvessel_ACh_mixed_afx)

#Checking sphericity and corrected results
summary(bloodvessel_ACh_mixed_afx)
 
#Bonferroni comparisons between treatments
estimate_contrasts(model = bloodvessel_ACh_mixed_afx, contrast = "Treatment", p_adjust = "bonferroni") |>
  filter(p < 0.05) |>
  display()





