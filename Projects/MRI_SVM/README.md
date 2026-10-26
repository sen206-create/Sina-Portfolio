# SVM analysis in R

These scripts cover the machine learning part of my MRI project, from assigning participants to cross-validation folds to fitting and evaluating the classifier. I prepared the MRI features beforehand using MATLAB, SPM and CAT.

| Script | What I did |
| --- | --- |
| [10_create_cv_folds.R](Scripts/10_create_cv_folds.R) | Divided the 175 participants into five folds, keeping the same proportion of controls and schizophrenia participants in each fold. |
| [11_build_linear_kernel.R](Scripts/11_build_linear_kernel.R) | Calculated the linear kernel from the grey matter voxel features. |
| [16_refit_linear_svm_qp.R](Scripts/16_refit_linear_svm_qp.R) | Trained the SVM with C = 1 and checked that each model had solved the optimisation correctly. |
| [17_complete_evaluation.R](Scripts/17_complete_evaluation.R) | Evaluated the predictions, ran permutation tests and investigated scanner differences, ghosting artifacts and a metadata comparison model. |

The numbering follows my working project. These are the final analysis scripts; the earlier fitting attempt and debugging scripts are not included here.

The scripts require locally prepared MRI features and participant information, which are not included in this repository.

[Return to the project write-up](schizophrenia-mri-svm.md)
