# Statistics with R

> **Work in progress:** A six-section project from my statistics class. Only Section 1 is included so far. Come back for new sections, updated results, and progress as I learn.

[View the R script](IQ_analysis.R)

## Section 1: Simple Linear Regression

I examined the relationship between brain size and PIQ scores in 39 observations., using R to fit a regression and examine how well it describes the data.

```r
PIQ_model = lm(PIQ ~ brain_size, data = IQ_dataset, na.action = na.exclude)
summary(PIQ_model)
```

**Tools:** R, tidyverse, easystats, readxl, and qqplotr.

![](Figures/slrplot.png)


- **Positive association:** Larger brain size was associated with higher PIQ scores (slope = 0.118, 95% CI [0.024, 0.212], p = .015).
- **Modest explanatory power:** The model accounted for approximately 15% of the variation in PIQ (R² = .150; adjusted R² = .127).
- **Bootstrap support:** The bootstrap confidence interval for the slope also excluded zero [0.029, 0.211] supporting positive association.
- **Residual normality:** The Shapiro–Wilk test did not detect a significant departure from normality (W = .962, p = .204).

The fitted equation was:

**Predicted PIQ = 4.061 + 0.118 × brain size**

## Section 2: Multiple Linear Regression

I extend the analysis beyond brain size alone, examining a model that includes brain size and height as predictors of PIQ.

- **Improved fit:** Adding height increased R² from .150 to .285. The improvement was statistically significant (p = .014), while the third model, gender, offered no significant further improvement (p = .554).
- **Preferred model:** The brain size + height model had the lowest AIC and BIC, with an adjusted R² of .245.
- **Predictor associations:** Holding the other predictor constant, brain size was positively associated with PIQ (standardised β = .659, p < .001), while height was negatively associated (β = −.458, p = .013).
- **Diagnostics:** The Shapiro–Wilk test did not detect a significant departure from residual normality (p = .458). Both predictors had VIF = 1.55, suggesting limited multicollinearity.
- **Influence:** Case 14 was flagged by the chosen influence screening rule and warrants further examination.
- **Sensitivity check:** Rank based regression also found a positive brain size association and a negative height association, both statistically significant.

The two-predictor model accounted for approximately **28.5% of the observed variation in PIQ**. The single regression model from section 1 only accounted for 15%.

**Work in progress:** These steps investigate model fit and assumptions; conclusions will be added after reviewing the outputs. Come back for results and the next section.

## Sections 3–6



**Coming as the class progresses.** Topics and code will be added here—check back to follow the project.
