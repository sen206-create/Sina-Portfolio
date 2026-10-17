# Exploring Schizophrenia Classification with Structural MRI and Support Vector Machines

## Why I started this project

After working on EEG and MRI preprocessing, I wanted to learn how to take brain images through to a machine learning analysis. I was interested in [Mourao-Miranda and colleagues' study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3315786/), which used structural MRI and support vector machines to predict the course of psychosis.

I used that paper as a starting point for learning the approach. My question was different: **could a model distinguish participants with schizophrenia from controls in an open MRI dataset?** 

This gave me a chance to learn more about SPM, CAT and the decisions involved in evaluating a classifier.

[View the R scripts and running instructions](MRI_SVM/README.md)

## The dataset and tools

I used T1 structural MRI scans from [UCLA's Consortium for Neuropsychiatric Phenomics dataset, OpenNeuro ds000030](https://openneuro.org/datasets/ds000030). My analysis included 175 participants: 50 with schizophrenia and 125 controls. The scans and diagnostic labels come from the dataset authors.

| Tool | What I used it for |
| --- | --- |
| MATLAB, SPM25 and CAT26 | MRI preprocessing, tissue segmentation and image inspection |
| R and hdf5r | Read the extracted voxel features and calculate the linear kernel |
| R and quadprog | Fit the SVM and check that the optimisation succeeded |
| R | Cross-validation, permutation testing and sensitivity analyses |

## From MRI scans to features

I processed the scans using CAT within SPM, producing modulated grey matter maps in a common space. I checked segmentation reports and alignment, reviewed unusual scans and collected the CAT quality ratings. I retained all 175 participants, keeping the quality and artifact flags for later analysis.

![Original T1 views and tissue segmentation from a CAT report, with its quality ratings](Projects/Figures/MRI_SVM/cat_segmentation.png)

This pilot report shows the original T1 views on the left and the segmentation on the right. Its weighted quality score was 65.43%, grade D. I have included it as an example of why I inspected the outputs rather than treating successful processing as a quality pass. The report header containing local file paths has been cropped out.

![Distribution of CAT quality scores across the 171 scans in the main processing batch](Projects/Figures/MRI_SVM/quality_distribution.png)

The distribution covers the 171 scans in the main batch. The four separately processed pilots are not included in this plot. These scores describe image and preprocessing quality, not classification performance.

I kept the common 1.5 mm grid and smoothed the maps at 8 mm FWHM. A common nonzero mask left 1,524,807 voxel values per participant. These became the model's features. I then calculated a 175 × 175 linear kernel matrix from the voxel values so the SVM could work with comparisons between participants.

## Training the SVM

I used a linear SVM with C = 1 and five folds. Each fold contained 25 controls and 10 participants with schizophrenia. Each model trained on 140 participants and predicted the remaining 35, so everyone received one prediction from a model that had not trained on their diagnosis.

One useful lesson came from the fitting process. The first implementation produced predictions, but numerical checks showed that it had not solved the SVM optimisation properly. I replaced the solver with quadratic programming, keeping the same features, folds and C value. All five corrected fits passed the optimisation checks.

## What the model found

![Confusion matrix showing 108 controls and 16 participants with schizophrenia correctly classified](Projects/Figures/MRI_SVM/confusion_matrix.png)

| Measure | Result |
| --- | ---: |
| Accuracy | 70.9% |
| Sensitivity for schizophrenia | 32.0% |
| Specificity for controls | 86.4% |
| Balanced accuracy | 59.2% |

The model identified controls more successfully than participants with schizophrenia. It missed 34 of the 50 schizophrenia cases. Because controls were the larger group, predicting CONTROL for everyone would already give 71.4% accuracy. Balanced accuracy, the average of sensitivity and specificity, was more informative here.

I also shuffled the labels 1,000 times and repeated the training and evaluation. I ran a second test that restricted shuffling within fold, scanner and ghosting groups.

![Permutation distribution with the observed balanced accuracy marked by a red line](Projects/Figures/MRI_SVM/permutation_test.png)

This is the second test. The red line marks the observed 59.2% balanced accuracy. Both tests gave p ≈ 0.034 after Holm correction for the two comparisons. This supports some discrimination under these testing assumptions, but does not establish clinical usefulness or show that the model learned only disease-related anatomy.

## Checking what might influence the result

![Balanced accuracy of the MRI model, metadata model and retrained sensitivity analyses](Projects/Figures/MRI_SVM/balanced_accuracy.png)

A logistic regression using age, gender, scanner and the ghosting flag reached 68.2% balanced accuracy, higher than the MRI model's 59.2%. When I excluded scans flagged for ghosting and retrained, MRI balanced accuracy was 52.0%. Training on one scanner and testing on the other gave 54.0% and 55.3%.

These comparisons made me more cautious about the original result. The subsets differ in size and composition, so the changes cannot be attributed to artifacts or scanners alone. I have not demonstrated that MRI adds predictive value beyond the metadata.

## What I learnt and the limitations

This project helped me understand how MRI preprocessing connects to machine learning, and why a single accuracy score is not enough. I learnt to check optimisation, compare against simple baselines and investigate how the dataset itself might influence predictions.

The common mask was created using all participants, including those later held out. Feature inclusion was therefore influenced by the test images. The evaluation is exploratory, with no external validation. A descriptive bootstrap gave a balanced accuracy interval of 52.2% to 66.4%, but it resampled the existing predictions and does not capture the full uncertainty from retraining or changing folds.

The reference paper studied illness course rather than schizophrenia versus controls. Its reported accuracy was the average of sensitivity and specificity, so the corresponding metric here is 59.2%, not 70.9%. I also used different preprocessing software, a 1.5 mm grid rather than 2 mm, and five-fold evaluation rather than paired leave-two-out. This project is an adaptation for learning, not a replication or a diagnostic tool.

The [code notes](MRI_SVM/README.md) describe the inputs and execution order. The [aggregate results](MRI_SVM/results/) contain the main metrics and sensitivity comparisons.
