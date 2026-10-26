# EEG slowing: group comparison

Participants: 100 PD and 49 controls.

Two prespecified features were compared using two-sided Welch t-tests. Bonferroni adjustment was applied across both tests (adjusted p = min(2 × raw p, 1)). Mean-difference confidence intervals are ordinary 95% intervals, not multiplicity-adjusted. Hedges' g uses the pooled standard deviation with the usual approximate bias correction.

## Relative theta power

PD mean: 0.2335; control mean: 0.1438.
PD minus control: 0.0897; 95% CI [0.0608, 0.1186].
Hedges' g: 0.839. Raw p: 7.64319e-09; Bonferroni-adjusted p: 1.52864e-08.

## Dominant frequency 4-13 Hz

PD mean: 7.5925; control mean: 7.2092.
PD minus control: 0.3833; 95% CI [-0.4322, 1.1988].
Hedges' g: 0.177. Raw p: 0.352225; Bonferroni-adjusted p: 0.704449.

## Interpretation limits

Higher theta power and lower dominant frequency are the expected directions for slowing. These are unadjusted group associations, not evidence of causation, cognitive decline, or individual diagnostic accuracy. Age, medication and other possible confounders are not controlled in these comparisons. The 4–13 Hz peak is a whole-scalp spectral maximum, not the paper's wavelet measure or a specifically occipital alpha peak.

Peak-at-search-boundary counts (False/True):

Peak at search boundary  False  True 
Group                                
Control                     35     14
PD                          91      9

## Methods references

https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html

Bonferroni controls family-wise error across the two planned comparisons. These group-comparison p-values are separate from the future classifier permutation test.
