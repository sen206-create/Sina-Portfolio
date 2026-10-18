# Raw linear kernel, accumulated in voxel blocks.
project = Sys.getenv("UCLA_MRI_PROJECT", unset = path.expand("~/ucla-mri-r"))
feature_dir = file.path(project, "feature-prep", "gm_1p5mm_fwhm8")
output_file = file.path(feature_dir, "linear_kernel_raw.rds")

build_kernel = function() {
  feature_data = hdf5r::H5File$new(
    file.path(feature_dir, "gm_features.mat"),
    mode = "r"
  )
  on.exit(feature_data$close_all())

  X = feature_data[["X"]]
  stopifnot(all(X$dims == c(175, 1524807)))

  kernel = matrix(0, 175, 175)
  block_size = 10000

  for (start in seq(1, 1524807, by = block_size)) {
    end = min(start + block_size - 1, 1524807)
    block = X$read(args = list(1:175, start:end))

    if (any(!is.finite(block))) {
      stop("Nonfinite voxel values encountered.")
    }

    kernel = kernel + tcrossprod(block)
  }

  kernel
}

kernel = build_kernel()
saveRDS(kernel, output_file)
print(dim(kernel))
