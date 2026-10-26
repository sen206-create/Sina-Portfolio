library(quadprog)

#Evaluating the SVM results
main = function() {
  project = Sys.getenv("UCLA_MRI_PROJECT", unset = path.expand("~/ucla-mri-r"))
  run_dir = readLines(file.path(project, "results", "evaluation", "latest_qp_run.txt"))[1]

  #Setting the number of permutations and bootstrap samples
  B = 1000L
  bootstrap_B = 5000L

  #Loading the saved predictions and participant information
  participants = read.csv(file.path(run_dir, "out_of_fold_predictions.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  kernel_file = file.path(project, "feature-prep", "gm_1p5mm_fwhm8", "linear_kernel_raw.rds")
  kernel = readRDS(kernel_file)
  participants = participants[order(participants$feature_row), ]
  y = ifelse(participants$diagnosis == "SCHZ", 1, -1)
  scanner = as.character(participants$ScannerSerialNumber)
  ghost = as.character(participants$ghost_NoGhost)
  out = file.path(run_dir, "extended_evaluation")
  dir.create(out, showWarnings = FALSE)

  #Calculating performance measures
  calculate_metrics = function(truth, guess) {
    sensitivity = if (any(truth == 1))
      mean(guess[truth == 1] == 1)
    else NA_real_
    specificity = if (any(truth == -1))
      mean(guess[truth == -1] == -1)
    else NA_real_
    c(accuracy = mean(truth == guess), sensitivity = sensitivity, specificity = specificity, balanced_accuracy = (sensitivity +
      specificity)/2)
  }

  #Fitting the same C = 1 SVM to a chosen training group
  predict_svm = function(training_rows, test_rows, labels) {
    z = labels[training_rows]
    n = length(training_rows)
    Q = kernel[training_rows, training_rows, drop = FALSE] * outer(z, z)
    fit = quadprog::solve.QP(Q, rep(1, n), cbind(z, diag(n), -diag(n)), c(0, rep(0, n), rep(-1, n)),
      meq = 1)
    a = fit$solution
    beta = a * z
    g = as.vector(kernel[training_rows, training_rows, drop = FALSE] %*% beta)
    low = a <= 1e-08
    high = a >= 1 - 1e-08
    free = !low & !high
    if (any(free))
      b = median(z[free] - g[free])
    else {
      lo = c((1 - g)[z == 1 & low], (-1 - g)[z == -1 & high])
      hi = c((1 - g)[z == 1 & high], (-1 - g)[z == -1 & low])
      if (!length(lo) || !length(hi))
        stop("Intercept interval unavailable.")
      b = (max(lo) + min(hi))/2
    }
    margin = z * (g + b)
    violation = ifelse(low, pmax(0, 1 - margin), ifelse(high, pmax(0, margin - 1), abs(margin - 1)))
    norm_sq = as.numeric(crossprod(beta, kernel[training_rows, training_rows, drop = FALSE] %*% beta))
    primal = 0.5 * norm_sq + sum(pmax(0, 1 - margin))
    dual = sum(a) - 0.5 * norm_sq
    if (max(c(0, -a, a - 1)) >= 1e-07 || abs(sum(beta)) >= 1e-07 || max(violation) >= 0.01 || abs(primal -
      dual) >= 1e-04 * max(1, abs(primal), abs(dual)))
      stop("SVM optimisation failed.")
    ifelse(as.vector(kernel[test_rows, training_rows, drop = FALSE] %*% beta + b) > 0, 1, -1)
  }

  #Running cross-validation for all participants or a subset
  cross_validate = function(labels, indices = seq_len(nrow(participants))) {
    guess = rep(NA_real_, nrow(participants))
    for (f in sort(unique(participants$fold[indices]))) {
      test_rows = indices[participants$fold[indices] == f]
      training_rows = indices[participants$fold[indices] != f]
      if (length(test_rows))
        guess[test_rows] = predict_svm(training_rows, test_rows, labels)
    }
    stopifnot(!anyNA(guess[indices]))
    guess
  }

  #Using the saved cross-validation predictions
  original = ifelse(participants$predicted_diagnosis == "SCHZ", 1, -1)
  observed = calculate_metrics(y, original)

  #Looking at performance by fold, scanner and ghosting flag
  groups = list(fold = as.character(participants$fold), scanner = scanner, ghost = ghost)
  subgroup = do.call(rbind, lapply(names(groups), function(v) {
    do.call(rbind, lapply(sort(unique(groups[[v]])), function(g) {
      i = which(groups[[v]] == g)
      cbind(data.frame(variable = v, group = g, n = length(i), n_CONTROL = sum(y[i] == -1), n_SCHZ = sum(y[i] ==
        1)), as.data.frame(as.list(calculate_metrics(y[i], original[i]))))
    }))
  }))
  write.csv(subgroup, file.path(out, "corrected_subgroup_metrics.csv"), row.names = FALSE)

  #Bootstrapping the fixed predictions, without retraining the models
  set.seed(202610171)
  strata = split(seq_along(y), y)
  boot = replicate(bootstrap_B, {
    i = unlist(lapply(strata, function(ix) ix[sample.int(length(ix), length(ix), replace = TRUE)]),
      use.names = FALSE)
    calculate_metrics(y[i], original[i])
  })
  intervals = data.frame(metric = names(observed), estimate = as.numeric(observed), lower = apply(boot,
    1, quantile, 0.025), upper = apply(boot, 1, quantile, 0.975), interpretation = "Descriptive fixed-prediction bootstrap; not full CV uncertainty")
  write.csv(intervals, file.path(out, "descriptive_bootstrap_intervals.csv"), row.names = FALSE)

  #Collecting the sensitivity analysis results
  sensitivity_results = list()
  summarise_group = function(name, i, guess) {
    cbind(data.frame(analysis = name, n = length(i),
      n_CONTROL = sum(y[i] == -1), n_SCHZ = sum(y[i] == 1)),
      as.data.frame(as.list(calculate_metrics(y[i], guess))))
  }

  #Repeating cross-validation without scans flagged for ghosting
  i = which(ghost == "No_ghost")
  no_ghost_predictions = cross_validate(y, i)
  sensitivity_results[["No_ghost"]] = summarise_group("No_ghost only, retrained five-fold", i, no_ghost_predictions[i])

  #Testing within each scanner and then training on the other scanner
  for (s in sort(unique(scanner))) {
    i = which(scanner == s)
    within = cross_validate(y, i)
    sensitivity_results[[paste("Within scanner", s)]] = summarise_group(paste("Within scanner", s), i, within[i])
    training_rows = which(scanner != s)
    cross = predict_svm(training_rows, i, y)
    sensitivity_results[[paste("Held-out scanner", s)]] = summarise_group(paste("Held-out scanner", s), i, cross)
  }

  #Fitting a logistic regression using participant and scan information
  age = suppressWarnings(as.numeric(participants$age))
  d = data.frame(response = as.integer(y == 1), age = age, gender = factor(participants$gender), scanner = factor(scanner),
    ghost = factor(ghost))
  baseline = rep(NA_real_, length(y))
  for (f in 1:5) {
    training_rows = which(participants$fold != f)
    test_rows = which(participants$fold == f)
    m = glm(response ~ age + gender + scanner + ghost,
      data = d[training_rows, ], family = binomial())
    if (!isTRUE(m$converged)) stop("Metadata model did not converge.")
    probability = predict(m, newdata = d[test_rows, ], type = "response")
    if (anyNA(probability))
      stop("Metadata baseline predictions missing.")
    baseline[test_rows] = ifelse(probability > 0.5, 1, -1)
  }
  sensitivity_results[["Metadata"]] = summarise_group("Metadata-only logistic baseline", seq_along(y), baseline)

  sensitivity_table = do.call(rbind, unname(sensitivity_results))
  write.csv(sensitivity_table, file.path(out, "sensitivity_metrics.csv"), row.names = FALSE)

  #Shuffling diagnosis labels within the two sets of permutation blocks
  schemes = list(within_fold = interaction(participants$fold, drop = TRUE), within_fold_scanner_ghost = interaction(participants$fold,
    scanner, ghost, drop = TRUE))

  permutation_results = list()
  for (scheme in names(schemes)) {
    blocks = split(seq_along(y), schemes[[scheme]])
    null_scores = numeric(B)

    for (j in seq_len(B)) {
      seed_offset = if (scheme == "within_fold") 0L else 100000L
      set.seed(20261017L + j + seed_offset)
      shuffled_labels = y
      for (i in blocks) {
        shuffled_labels[i] = y[i][sample.int(length(i))]
      }
      guess = cross_validate(shuffled_labels)
      null_scores[j] = unname(calculate_metrics(shuffled_labels, guess)["balanced_accuracy"])
    }
    exceed = sum(null_scores >= observed["balanced_accuracy"] - 1e-12)
    monte_carlo_interval = binom.test(exceed, B)$conf.int
    permutation_results[[scheme]] = data.frame(test = scheme, B = B, observed_BA = observed["balanced_accuracy"],
      null_mean = mean(null_scores), exceedances = exceed, p_value = (exceed + 1)/(B + 1), monte_carlo_tail_lower = monte_carlo_interval[1],
      monte_carlo_tail_upper = monte_carlo_interval[2])
    pdf(file.path(out, paste0("null_distribution_", scheme, ".pdf")))
    hist(null_scores, breaks = 25, main = scheme, xlab = "Null balanced accuracy", col = "grey80",
      xlim = range(c(null_scores, observed["balanced_accuracy"])))
    abline(v = observed["balanced_accuracy"], col = "red", lwd = 2)
    dev.off()
  }

  #Applying Holm correction to the two permutation tests
  perm = do.call(rbind, permutation_results)
  perm$p_holm_two_tests = p.adjust(perm$p_value, method = "holm")
  write.csv(perm, file.path(out, "permutation_results.csv"), row.names = FALSE)

  print(intervals, row.names = FALSE)
  print(sensitivity_table, row.names = FALSE)
  print(perm, row.names = FALSE)
}
main()