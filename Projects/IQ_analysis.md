# Statistics with R

> **Work in progress:** A six-section project from my statistics class. Only Section 1 is included so far. Come back for new sections, updated results, and progress as I learn.

[View the R script](IQ_analysis.R)

## Section 1: Simple Linear Regression

I’m exploring the relationship between brain size and PIQ score, using R to fit a regression and examine how well it describes the data.

```r
PIQ_model = lm(PIQ ~ brain_size, data = IQ_dataset, na.action = na.exclude)
summary(PIQ_model)
```

### What the script covers

- **Visualisation:** Scatterplot with a fitted regression line.
- **Model interpretation:** Coefficients, confidence intervals, R², and statistical tests.
- **Diagnostics:** Standardised residuals, Q–Q plots, and a Shapiro–Wilk test.
- **Influence:** Identify observations that strongly affect the model.
- **Bootstrapping:** Explore uncertainty in parameter estimates through resampling.

**Tools:** R, tidyverse, easystats, readxl, and qqplotr.

![](Figures/slrplot.png)

## Sections 2–6

**Coming as the class progresses.** Topics and code will be added here—check back to follow the project.
