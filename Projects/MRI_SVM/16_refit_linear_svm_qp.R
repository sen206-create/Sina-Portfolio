#Loading the participant folds and linear kernel
project = Sys.getenv("UCLA_MRI_PROJECT", unset = path.expand("~/ucla-mri-r"))
if (!requireNamespace("quadprog", quietly = TRUE)) stop("Run install.packages('quadprog') first.")
base = file.path(project, "results", "evaluation")
participants = read.csv(file.path(base, "cv_folds.csv"), stringsAsFactors = FALSE)
participants = participants[order(participants$feature_row), ]
kernel = readRDS(file.path(project, "feature-prep", "gm_1p5mm_fwhm8", "linear_kernel_raw.rds"))
stopifnot(
  identical(as.integer(participants$feature_row), 1:175),
  identical(dim(kernel), c(175L, 175L)),
  all(is.finite(kernel)),
  isSymmetric(kernel),
  all(participants$diagnosis %in% c("CONTROL", "SCHZ")),
  setequal(participants$fold, 1:5)
)
counts = table(participants$fold, participants$diagnosis)
stopifnot(
  all(counts[, "CONTROL"] == 25),
  all(counts[, "SCHZ"] == 10)
)

#Setting up the linear SVM
C = 1
predictions = rep(NA_character_, 175)
scores = rep(NA_real_, 175)
models = diagnostics = vector("list", 5)

#Training on four folds and predicting the remaining fold
for (f in 1:5) {
  training_rows = which(participants$fold != f)
  test_rows = which(participants$fold == f)
  y = ifelse(participants$diagnosis[training_rows] == "SCHZ", 1, -1)
  n = length(training_rows)
  Q = kernel[training_rows, training_rows] * outer(y, y)
  fit = quadprog::solve.QP(Q, rep(1, n), cbind(y, diag(n), -diag(n)), c(0, rep(0, n), rep(-C, n)),
    meq = 1)

  #Calculating the SVM coefficients and intercept
  alpha = fit$solution
  beta = alpha * y
  g = as.vector(kernel[training_rows, training_rows] %*% beta)
  bound_tol = 1e-08
  low = alpha <= bound_tol
  high = alpha >= C - bound_tol
  free = !low & !high
  if (any(free)) {
    intercept = median(y[free] - g[free])
  } else {
    lower = c((1 - g)[y == 1 & low], (-1 - g)[y == -1 & high])
    upper = c((1 - g)[y == 1 & high], (-1 - g)[y == -1 & low])
    if (!length(lower) || !length(upper))
      stop("Cannot determine finite intercept.")
    intercept = (max(lower) + min(upper))/2
  }

  #Checking that the solver has found a valid solution
  margin = y * (g + intercept)
  violation = ifelse(low, pmax(0, 1 - margin), ifelse(high, pmax(0, margin - 1), abs(margin - 1)))
  norm_sq = as.numeric(crossprod(beta, kernel[training_rows, training_rows] %*% beta))
  primal = 0.5 * norm_sq + C * sum(pmax(0, 1 - margin))
  dual = sum(alpha) - 0.5 * norm_sq
  gap = primal - dual
  box = max(c(0, -alpha, alpha - C))
  equality = abs(sum(beta))
  passed = box < 1e-07 && equality < 1e-07 && max(violation) < 0.01 && abs(gap) < 1e-04 * max(1, abs(primal),
    abs(dual))
  diagnostics[[f]] = data.frame(fold = f, box_violation = box, equality_residual = equality, max_margin_violation = max(violation),
    primal = primal, dual = dual, objective_gap = gap, passed = passed)
  print(diagnostics[[f]], row.names = FALSE)
  if (!passed)
    stop("Optimisation check failed. No performance reported.")

  #Predicting the participants left out of training
  scores[test_rows] = as.vector(kernel[test_rows, training_rows, drop = FALSE] %*% beta + intercept)
  predictions[test_rows] = ifelse(scores[test_rows] > 0, "SCHZ", "CONTROL")
  models[[f]] = list(training_rows = training_rows, test_rows = test_rows, alpha = alpha, beta = beta,
    intercept = intercept, C = C, solver_iterations = fit$iterations)
}
stopifnot(!anyNA(predictions))

#Calculating accuracy, sensitivity, specificity and balanced accuracy
calculate_metrics = function(truth, guess) {
  sensitivity = mean(guess[truth == "SCHZ"] == "SCHZ")
  specificity = mean(guess[truth == "CONTROL"] == "CONTROL")
  data.frame(accuracy = mean(truth == guess), sensitivity = sensitivity, specificity = specificity,
    balanced_accuracy = (sensitivity + specificity)/2)
}

#Comparing the SVM with always predicting control
summary = rbind(cbind(model = "Linear SVM C=1 (QP)", calculate_metrics(participants$diagnosis, predictions)),
  cbind(model = "Always CONTROL", calculate_metrics(participants$diagnosis, rep("CONTROL", 175))))
fold_metrics = do.call(rbind, lapply(1:5, function(f) {
  i = participants$fold == f
  cbind(fold = f, calculate_metrics(participants$diagnosis[i], predictions[i]))
}))

#Saving the predictions, results and fitted models
out = file.path(base, paste0("linear_svm_C1_QP_", format(Sys.time(), "%Y%m%d_%H%M%S")))
if (dir.exists(out)) stop("Output already exists; retry in a second.")
dir.create(out, recursive = TRUE)
participants$predicted_diagnosis = predictions
participants$correct = predictions == participants$diagnosis
participants$decision_score = scores
write.csv(participants, file.path(out, "out_of_fold_predictions.csv"), row.names = FALSE)
write.csv(summary, file.path(out, "summary_metrics.csv"), row.names = FALSE)
write.csv(fold_metrics, file.path(out, "fold_metrics.csv"), row.names = FALSE)
write.csv(do.call(rbind, diagnostics), file.path(out, "optimisation_diagnostics.csv"), row.names = FALSE)
saveRDS(models, file.path(out, "fold_models_qp.rds"))

cat("\nAll five optimisation checks passed.\n")
print(table(Actual = participants$diagnosis, Predicted = factor(predictions, levels = c("CONTROL", "SCHZ"))))
print(summary, row.names = FALSE)
writeLines(normalizePath(out), file.path(base, "latest_qp_run.txt"))
cat("\nSaved results to:", out, "\n")
