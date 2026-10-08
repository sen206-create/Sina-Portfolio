library(tidyverse)
library(easystats)
library(readxl)

Oestrogen_bone_density = read_excel("Oestrogen bone density.xlsx")

#Data is in wide format. Pivot to long so all doses in a single column
oestrogen_data = Oestrogen_bone_density |>
  pivot_longer(cols = c("0", "1", "4", "40", "400", "4000"), names_to = "dose", values_to = "bone_density") |>
  drop_na(bone_density) |>
  mutate(dose=as_factor(dose))

oestrogen_data |>
  group_by(dose) |>
  describe_distribution(ci = 0.95, iterations = 1000) |>
  display()

#Box plot to compare distribution of data
ggplot(oestrogen_data, aes(x = dose, y = bone_density, fill = dose)) +
  geom_boxplot(alpha = .3) +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_x_discrete(labels = c("0", "1", "4", "40", "400", "4000")) +
  labs(x = "Oestrogen dose (mg/kg/day)", y = "Bone Density (mg/square cm)") +
  theme_bw(base_size = 12)

#Fitting anova
oestrogen_model = lm(bone_density ~ dose, data = oestrogen_data)
anova(oestrogen_model) |>
  model_parameters(es = "omega", ci = 0.95) |>
  display()

#Standardising residuals
oestrogen_data = oestrogen_data |>
  mutate(st_residual = rstandard(oestrogen_model))

#Histogram to see distribution of standardised residuals
ggplot(oestrogen_data, aes(x = st_residual)) +
  geom_histogram(aes(y = after_stat(density)), bins = 10, colour = "black", fill = "turquoise") +
  stat_function(fun = dnorm, args = list(
    mean = mean(oestrogen_data$st_residual),
    sd = sd(oestrogen_data$st_residual))) +
  labs(x ="Standardised Residuals", y = "Density") +
  theme_bw()

library(qqplotr)

ggplot(oestrogen_data, aes(sample = st_residual)) +
  stat_qq_band(fill = "orange", alpha = 0.2) +
  stat_qq_line() +
  stat_qq_point() +
  labs(x= "Theoretical Quantiles", y = "Standard Residual") +
  theme_bw(base_size = 12)

shapiro.test(oestrogen_data$st_residual)

cumulative_percent = function(var, cut_off = 1.96){
  ecdf_proportions = abs(var) |>
    ecdf()
  100*(1 - ecdf_proportions(cut_off))
}

oestrogen_data |>
  summarize(
    z_1.96 = cumulative_percent(st_residual),
    z_2.58 = cumulative_percent(st_residual, cut_off = 2.58),
    z_3.29 = cumulative_percent(st_residual, cut_off = 3.29)
  )

#Finding influential data points
oestrogen_inf = influence.measures(oestrogen_model)$infmat |>
  as_tibble() |>
  rowid_to_column(var = "case_no")

oestrogen_inf |>
  filter(
    if_any(starts_with("df"), ~ abs(.) > 1)
  )

#Levene test for equality of variance assumption in anova
library(car)
leveneTest(oestrogen_data$bone_density, oestrogen_data$dose, center = median)

#Rank based linear regression
library(Rfit)
oestrogen_rfit_model = with(oestrogen_data, oneway.rfit(bone_density, dose))
oestrogen_rfit_model

#Planned Contrasts between mid to high doses for question 2.
mid_doses_vs_high_doses = c(0, 0, -1/4, -1/4, 1/4, 1/4)
low_doses_vs_mid_high_doses = c(-2/6, -2/6, 1/6, 1/6, 1/6, 1/6)
dose_0_vs_1 = c(-1/2, 1/2, 0, 0, 0, 0)
dose_4_vs_40 = c(0, 0, -1/2, 1/2, 0, 0)
dose_400_vs_4000 = c(0, 0, 0, 0, -1/2, 1/2)

contrasts(oestrogen_data$dose) = cbind(mid_doses_vs_high_doses, low_doses_vs_mid_high_doses, dose_0_vs_1, dose_4_vs_40, dose_400_vs_4000)
contrasts(oestrogen_data$dose)

oestrogen_model = lm(bone_density ~ dose, data = oestrogen_data)
anova(oestrogen_model) |>
  model_parameters(es = "omega", ci = 0.95) |>
  display()

model_parameters(oestrogen_model) |>
  display()

# Finding trend for bone density to change with oestrogen dose
contrasts(oestrogen_data$dose) = contr.poly(6)
oestrogen_poly = lm(bone_density ~ dose, data = oestrogen_data)
model_parameters(oestrogen_poly) |>
  display()