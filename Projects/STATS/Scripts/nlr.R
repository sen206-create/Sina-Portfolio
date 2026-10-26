#Setting up data for analysis
library(easystats)
library(tidyverse)
library(readxl)

puromycin_inhibition = read_excel("puromycin inhibition.xlsx")
puromycin_inhibition = puromycin_inhibition |>
  rename(Substrate_concentration = 'Substrate concentration (ppm)',
         Reaction_rate = 'Reaction rate (counts/min/min)')

inhibition = puromycin_inhibition |>
  filter(condition == "treated")

control = puromycin_inhibition |>
  filter(condition == "untreated")

ggplot(puromycin_inhibition, aes(x = Substrate_concentration, y = Reaction_rate, colour = condition)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = 'Substrate concentration (ppm)', y = ' Reaction rate (counts/min/min)')
  theme_bw(base_size = 14)

# Non-linear regression with the equation for reaction rate = Vmax x [S]/(KM + [S]
inhibition_model = nls(Reaction_rate ~ Vmax*Substrate_concentration/(Km + Substrate_concentration),
      data = inhibition,
      start = list(Vmax = 200, Km = 0.1))
summary(inhibition_model)

# Residuals and fitted values
inhibition$residual = residuals(inhibition_model)
inhibition = mutate(inhibition, st_residual = (residual - mean(residual))/sd
                     (residual))
inhibition$predicted = fitted(inhibition_model)

#Checking for deviations of residuals from normality and homoscedasticity
ggplot(inhibition, aes(x = st_residual)) +
  geom_histogram(aes(y = after_stat(density)), bins = 8, colour = "black",
                 fill = "magenta") +
  stat_function(fun = dnorm, args = list(
    mean = mean(inhibition$st_residual),
    sd = sd(inhibition$st_residual))) +
  labs(x ="standardised residuals", y = "density") +
  theme_bw()

ggplot(inhibition, aes(sample = st_residual)) +
  geom_qq() +
  geom_qq_line() +
  labs(x = "theoretical quantiles", y = "standardised residuals")

ggplot(inhibition, aes(x = predicted, y = st_residual)) +
  geom_point() +
  labs(x = "predicted values", y = "standardised residuals") +
  geom_abline(slope = 0, intercept = 0, colour = "black") +
  theme_bw()
shapiro.test(inhibition$st_residual)

# Nonlinear regression for the control data 

control_model = nls(Reaction_rate ~ Vmax*Substrate_concentration/(Km + Substrate_concentration),
      data = control,
      start = list(Vmax = 150, Km = 0.1))
summary(control_model)

confint(control_model)

control$residual = residuals(control_model)
control = mutate(control, st_residual = (residual - mean(residual))/sd(residual))
control$predicted = fitted(control_model)

ggplot(control, aes(x = st_residual)) +
  geom_histogram(aes(y = after_stat(density)), bins = 8, colour = "black",
                 fill = "magenta") +
  stat_function(fun = dnorm, args = list(
    mean = mean(control$st_residual),
    sd = sd(control$st_residual))) +
  labs(x ="standardised residuals", y = "density") +
  theme_bw()

ggplot(control, aes(sample = st_residual)) +
  geom_qq() +
  geom_qq_line() +
  labs(x = "theoretical quantiles", y = "standardised residuals")

shapiro.test(control$st_residual)

ggplot(control, aes(x = predicted, y = st_residual)) +
  geom_point() +
  labs(x = "predicted values", y = "standardised residuals") +
  geom_abline(slope = 0, intercept = 0, colour = "black") +
  theme_bw()


#Plotting predicted data
control_pred_data = data.frame(Substrate_concentration = seq(0.02, 1.1, length.out = 1000))
predicted_density = predict(control_model, control_pred_data)
control_pred_data = mutate(control_pred_data, Reaction_rate = predicted_density)

treated_pred_data = data.frame(Substrate_concentration = seq(0.02, 1.1, length.out = 1000))
treated_predicted = predict(inhibition_model, treated_pred_data)
treated_pred_data = mutate(treated_pred_data, Reaction_rate = treated_predicted)

ggplot() +
  geom_point(data = inhibition, aes(x = Substrate_concentration, y = Reaction_rate), colour = "blue") +
  geom_line(data = treated_pred_data, aes(x = Substrate_concentration, y = treated_predicted), colour = "blue") +
  geom_point(data = control, aes(x = Substrate_concentration, y = Reaction_rate), colour = "red") +
  geom_line(data = control_pred_data, aes(x = Substrate_concentration, y = Reaction_rate), colour = "red") +
  labs(x = 'Substrate concentration (ppm)', y = ' Reaction rate (counts/min/min)')
  theme_bw(base_size = 14)
