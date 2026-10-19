library(hdf5r)

#Setting up the feature file
project = Sys.getenv("UCLA_MRI_PROJECT", unset = path.expand("~/ucla-mri-r"))
feature_dir = file.path(project, "feature-prep", "gm_1p5mm_fwhm8")
output_file = file.path(feature_dir, "linear_kernel_raw.rds")

#Reading the voxel features in blocks to limit memory use
build_kernel = function() {
  feature_data = hdf5r::H5File$new(file.path(feature_dir, "gm_features.mat"), mode = "r")
  on.exit(feature_data$close_all())
  X = feature_data[["X"]]
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
