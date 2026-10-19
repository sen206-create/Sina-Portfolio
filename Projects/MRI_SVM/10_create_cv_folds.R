project = Sys.getenv("UCLA_MRI_PROJECT", unset = path.expand("~/ucla-mri-r"))
feature_dir = file.path(project, "feature-prep", "gm_1p5mm_fwhm8")
participants = read.csv(file.path(feature_dir, "feature_participants.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)

stopifnot(setequal(participants$feature_row, seq_len(nrow(participants))))
participants = participants[order(participants$feature_row), , drop = FALSE]

seed = 20261017L
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(seed)
participants$fold = NA_integer_
for (group in c("CONTROL", "SCHZ")) {
  indices = which(participants$diagnosis == group)
  shuffled = sample(indices, length(indices), replace = FALSE)
  participants$fold[shuffled] = rep(1:5, length.out = length(shuffled))
}

output_dir = file.path(project, "results", "evaluation")
write.csv(participants, file.path(output_dir, "cv_folds.csv"), row.names = FALSE)

print(with(participants, table(fold, diagnosis)))
