#Setting up the participant data
project = Sys.getenv("UCLA_MRI_PROJECT", unset = path.expand("~/ucla-mri-r"))
feature_dir = file.path(project, "feature-prep", "gm_1p5mm_fwhm8")
participants = read.csv(file.path(feature_dir, "feature_participants.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)

#Checking participant IDs and feature rows
stopifnot(all(c("participant_id", "diagnosis", "feature_row") %in% names(participants)))
stopifnot(nrow(participants) == 175L, !anyNA(participants$participant_id), !anyDuplicated(participants$participant_id),
    !anyNA(participants$diagnosis), all(participants$diagnosis %in% c("CONTROL", "SCHZ")))
stopifnot(setequal(participants$feature_row, seq_len(nrow(participants))))
participants = participants[order(participants$feature_row), , drop = FALSE]
stopifnot(sum(participants$diagnosis == "CONTROL") == 125L, sum(participants$diagnosis == "SCHZ") ==
    50L)

#Creating five folds with 25 controls and 10 schizophrenia participants in each
seed = 20261017L
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(seed)
participants$fold = NA_integer_
for (group in c("CONTROL", "SCHZ")) {
    indices = which(participants$diagnosis == group)
    shuffled = sample(indices, length(indices), replace = FALSE)
    participants$fold[shuffled] = rep(1:5, length.out = length(shuffled))
}

#Saving the folds for the SVM analysis
output_dir = file.path(project, "results", "evaluation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file = file.path(output_dir, "cv_folds.csv")
if (file.exists(output_file)) {
    previous = read.csv(output_file, stringsAsFactors = FALSE)
    keys = c("participant_id", "feature_row", "diagnosis", "fold")
    if (!identical(previous[keys], participants[keys])) {
        stop("A different split file already exists; do not overwrite established folds.")
    }
    message("Existing splits match; retained cv_folds.csv.")
} else {
    write.csv(participants, output_file, row.names = FALSE, na = "")
}
writeLines(c(paste("Seed:", seed), "RNG: Mersenne-Twister / Inversion / Rejection", "Five folds stratified by diagnosis, at participant level.",
    "Each fold: 25 controls and 10 SCHZ; train on other four folds.", "Planned model: linear SVM, C=1; learned transforms within training folds.",
    "Primary planned metric: balanced accuracy; report sensitivity and specificity.", "Mask uses the full cohort: test images influenced feature inclusion.",
    "These folds are not a scanner-held-out or external-validation design.", "No model training performed by this script."),
    file.path(output_dir, "cv_design.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "split_session_info.txt"))
print(with(participants, table(fold, diagnosis)))
cat("Saved splits:", output_file, "\n")
