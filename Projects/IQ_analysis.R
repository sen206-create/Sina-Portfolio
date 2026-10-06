library(easystats)
library(tidyverse)
library(readxl)

iq_data = read_excel("IQ dataset.xlsx")

ggplot(IQ_dataset, aes(x = brain_size, y = PIQ, colour = gender))+
  geom_point()+
  geom_smooth(method=lm, colour = "black", fill = 'orange', alpha = 0.3) +
  labs(x = "brain size(arbitrary units)", y = "Procedural IQ")+
  theme_bw()

#Fit a simple linear regression
PIQ_model = lm(PIQ ~ brain_size, data = IQ_dataset, na.action = na.exclude)
summary(PIQ_model)

#Find the R^2 values, p-values, F-statisitc
model_performance(PIQ_model) |>
  display(digits = 3)
test_wald(PIQ_model) |>
  display(digits = 3)
model_parameters(PIQ_model) |>
  display(digits = 3)


#Add the standard residuals and predicted values to dataset
IQ_dataset = IQ_dataset |>
  mutate(st_residual = rstandard(PIQ_model))
IQ_dataset = IQ_dataset |>
  mutate(predicted = PIQ_model$fitted.values)
head(IQ_dataset)


ggplot(IQ_dataset, aes(x=brain_size, y = st_residual))+
  geom_point()+
  labs(x = "predicted values", y = "standardised residuals")+
  geom_abline(slope = 0, intercept = 0, colour = 'black', linetype = 'dashed')+
  theme_bw()

ggplot(IQ_dataset, aes(x=brain_size))+
         geom_histogram(aes(y = after_stat(density)), bins = 10, colour = 'black', fill = 'turquoise') +
         stat_function(fun= dnorm, args = list(
           mean = mean(IQ_dataset$st_residual),
           sd = sd(IQ_dataset$st_residual)))+
             labs(x='Brain Size', y = 'Density') +
             theme_bw()
library(qqplotr)
ggplot(IQ_dataset, aes(sample = st_residual))+
  stat_qq_band(fill = 'orange', alpha = 0.2)+
  stat_qq_line()+
  stat_qq_point()+
  labs(x = 'theoretical quantiles', y = 'standardised residuals')+
  theme_bw()
shapiro.test(IQ_dataset$st_residual)

IQ_dataset |>
  summarise(
    z_1.96 = cummulative_percent(st_residual),
    z_2.58 = cummulative_percent(st_residual, cut_off = 2.58),
    z_3.29 = cummulative_percent(st_residual, cut_off = 3.29)
  )

PIQ_slr_inf = influence.measures(PIQ_model)$infmat |>
  as_tibble() |>
  rowid_to_column(var = "case no")
head(PIQ_slr_inf)

set.seed(100)
model_parameters(PIQ_model, bootstrap = TRUE, iteration = 2000) |>
  display(digits = 3)
