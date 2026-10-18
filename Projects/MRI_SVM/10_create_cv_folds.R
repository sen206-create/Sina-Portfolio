# Fixed five-fold splits, stratified by diagnosis.

project = Sys.getenv("UCLA_MRI_PROJECT", unset = path.expand("~/ucla-mri-r"))
participants = read.csv(
  file.path(project, "feature-prep", "gm_1p5mm_fwhm8", "feature_participants.csv")
)

stopifnot(
  all(c("participant_id", "diagnosis", "feature_row") %in% names(participants)),
  nrow(participants) == 175,
  !anyNA(participants$participant_id),
  !anyDuplicated(participants$participant_id),
  !anyNA(participants$diagnosis),
  sum(participants$diagnosis == "CONTROL") == 125,
  sum(participants$diagnosis == "SCHZ") == 50,
  setequal(participants$feature_row, 1:175)
)

participants = participants[order(participants$feature_row), ]

set.seed(20261017)
participants$fold = NA

for (group in c("CONTROL", "SCHZ")) {
  indices = which(participants$diagnosis == group)
  participants$fold[sample(indices)] = rep(1:5, length.out = length(indices))
}

output_dir = file.path(project, "results", "evaluation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file = file.path(output_dir, "cv_folds.csv")

if (file.exists(output_file)) {
  previous = read.csv(output_file)
  keys = c("participant_id", "feature_row", "diagnosis", "fold")

  if (!identical(previous[keys], participants[keys])) {
    stop("Saved folds differ. Existing file was not overwritten.")
  }
} else {
  write.csv(participants, output_file, row.names = FALSE)
}

print(with(participants, table(fold, diagnosis)))
