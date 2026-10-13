library(readxl)
library(tidyverse)
library(easystats)
library(gmodels)

arthritis_activity_data = read_excel("arthritis activity data.xlsx")
male = subset(arthritis_activity_data, sex == "Male")
female = subset(arthritis_activity_data, sex == "Female")

#Expected frequencies in the individual contingency tables
CrossTable(male$treatment, male$outcome, expected = TRUE, sresid = TRUE, format = "SPSS")

fisher.test(table(male$treatment, male$outcome)) #Fishers test as some counts are below 5

CrossTable(female$treatment, female$outcome, expected = TRUE, sresid = TRUE, format = "SPSS")

#3-way contingency
arthritis_sex_contingency = xtabs(~ treatment + outcome + sex, data = arthritis_activity_data)
arthritis_sex_contingency

#Log-linear analysis
library(MASS)
arthritis_sex_saturated = loglm( ~ treatment * outcome * sex, data = arthritis_sex_contingency, fit = TRUE)
arthritis_sex_saturated

#Remove 3-way interaction
arthritis_sex_three_way = loglm( ~ treatment + outcome + sex + treatment:outcome + treatment:sex + outcome:sex, data = arthritis_sex_contingency, fit = TRUE)
arthritis_sex_three_way

#Test significance of change against saturated model
anova(arthritis_sex_saturated, arthritis_sex_three_way)

# p 0.598, and so >0.05 so the three-way interaction treatment:outcome:sex can be removed from the model. 

#Test two-way interactions by removing
arthritis_sex_minus_treatment_outcome = loglm( ~ treatment + outcome + sex + treatment:sex + outcome:sex, data = arthritis_sex_contingency, fit = TRUE)
arthritis_sex_minus_treatment_outcome

anova(arthritis_sex_three_way, arthritis_sex_minus_treatment_outcome)

#Removing the treatment:outcome interaction causes a significant difference in fit, so it needs to stay in the final model
arthritis_sex_minus_treatment_sex = loglm( ~ treatment + outcome + sex + treatment:outcome + outcome:sex, data = arthritis_sex_contingency, fit = TRUE)
arthritis_sex_minus_treatment_sex

anova(arthritis_sex_three_way, arthritis_sex_minus_treatment_sex)

#There is no significant change in fit when treatment:sex is removed from the analysis, so it should be removed from the final model.
arthritis_sex_minus_outcome_sex = loglm( ~ treatment + outcome + sex + treatment:outcome + treatment:sex, data = arthritis_sex_contingency, fit = TRUE)
arthritis_sex_minus_outcome_sex

anova(arthritis_sex_three_way, arthritis_sex_minus_outcome_sex)

mosaicplot(arthritis_sex_contingency, shade = TRUE, main = "Arthritis treatment and outcomes")


