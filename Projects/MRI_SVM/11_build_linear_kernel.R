#Setting up the feature file
project = Sys.getenv("UCLA_MRI_PROJECT", unset = path.expand("~/ucla-mri-r"))
feature_dir = file.path(project, "feature-prep", "gm_1p5mm_fwhm8")
output_file = file.path(feature_dir, "linear_kernel_raw.rds")
if (file.exists(output_file)) {
  stop("Kernel already exists. Keep it and proceed to script 16; no rerun needed.")
}
if (!requireNamespace("hdf5r", quietly = TRUE)) {
  stop("Install hdf5r with install.packages('hdf5r'), then rerun.")
}

#Reading the voxel features in blocks to limit memory use
build_kernel = function() {
  feature_data = hdf5r::H5File$new(file.path(feature_dir, "gm_features.mat"), mode = "r")
  on.exit(feature_data$close_all())
  X = feature_data[["X"]]
  stopifnot(all(X$dims == c(175, 1524807)))
  kernel = matrix(0, 175, 175)
  block_size = 10000L
  for (start in seq.int(1L, 1524807L, by = block_size)) {
    end = min(start + block_size - 1L, 1524807L)
    block = X$read(args = list(1:175, start:end))
    if (any(!is.finite(block)))
      stop("Nonfinite voxel values encountered.")
    kernel = kernel + tcrossprod(block)
  }
  kernel
}

#Saving the linear kernel matrix
kernel = build_kernel()
saveRDS(kernel, output_file)
print(dim(kernel))
