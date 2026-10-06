#Independent project to showcase skills learnt in Neurohacking for R.
library(oro.nifti)
library(oro.dicom)
n4_exe = "/opt/anaconda3/envs/ants/bin/N4BiasFieldCorrection"
extraction_exe = "/opt/anaconda3/envs/ants/bin/antsBrainExtraction.sh"
registration_exe = "/opt/anaconda3/envs/ants/bin/antsRegistrationSyN.sh"

Sys.setenv(
  PATH = paste("/opt/anaconda3/envs/ants/bin", Sys.getenv("PATH"), sep = ":")
)
stopifnot(
  file.exists(n4_exe),
  file.exists(extraction_exe),
  file.exists(registration_exe)
)

#FETCHING FILES FOR T1
raw_dir = "~/Desktop/Neurohacking_data-master/BRAINIX/DICOM/T1"
processed_dir = "~/Desktop/Neurohacking_data-master/BRAINIX/Processed"
template_dir = "~/Desktop/Neurohacking_data-master/Template"

all_slices_T1 = readDICOM(raw_dir)

#CONVERSION OF T1 TO NIfTI + BASIC QUALITY CHECK
nii_T1=dicom2nifti(all_slices_T1)
par(mfrow = c(1, 3))
z_range <- range(nii_T1)
# T1 Dimensions: 512x512 with 22 slices
# T1 Pixel dimensions: 0.46875mmx0.46875mmx5mm
# T1 Range: 0-1884
image(
  nii_T1[,,6], 
  col = gray(0:64/64), 
  zlim = z_range,
  asp = 1,
  main="Slice 6")
image(nii_T1[,,12],
      col = gray(0:64/64),
      zlim = z_range,
      asp = 1,
      main="Slice 12")
image(nii_T1[,,18],
      col = gray(0:64/64),
      zlim = z_range,
      asp = 1,
      main="Slice 18")
par(mfrow = c(1, 1))

#SAVING FILE AS NIfTI
writeNIfTI(nim = nii_T1, filename = file.path(processed_dir, "T1_converted"))

#N4 BIAS CORRECTION FOR T1
##N4 altered voxel intensities but there were subtle visual differences on the slices after correction. Dimensions and voxel spacing remain unchanged.

input_file = path.expand(file.path(processed_dir, "T1_converted.nii.gz"))
output_file = path.expand(file.path(processed_dir, "n4_T1_corrected.nii.gz"))

status = system2(
  command = n4_exe,
  args = c(
    "-d", "3",
    "-i", shQuote(input_file),
    "-o", shQuote(output_file),
    "-s", "2",
    "-v", "1"
  )
)
stopifnot(status == 0, file.exists(output_file))
n4img = readNIfTI(output_file, reorient=FALSE)

par(mfrow = c(1,2))
pre_n4img = readNIfTI(input_file, reorient=FALSE)
image(
  pre_n4img[,,18], 
  col=gray(0:64/64), 
  zlim=z_range,
  asp = 1,
  main="Original")
image(
  n4img[,,18], 
  col=gray(0:64/64), 
  zlim=z_range,
  asp = 1,
  main="N4 Corrected")
par(mfrow = c(1,1))

#BRAIN EXTRACTION FOR T1
template_file = path.expand(file.path(template_dir, "JHU_MNI_SS_T1.nii.gz"))
mask_file = path.expand(file.path(template_dir, "JHU_MNI_SS_T1_mask.nii.gz"))
# Template and template mask dimensions: 181 x 217 x 181 voxels

extracted_prefix = path.expand(file.path(processed_dir, "extracted-"))

status = system2(
  command = extraction_exe, 
  args = c(
    "-d", "3",
    "-a", shQuote(input_file),
    "-e", shQuote(template_file),
    "-m", shQuote(mask_file),
    "-o", shQuote(extracted_prefix)
  )
)

stopifnot(status == 0)
brain_file <- paste0(extracted_prefix, "BrainExtractionBrain.nii.gz")
brain_mask_file <- paste0(extracted_prefix, "BrainExtractionMask.nii.gz")
stopifnot(
  file.exists(brain_file),
  file.exists(brain_mask_file)
)
extracted_brain = readNIfTI(brain_file, reorient=FALSE)
brain_mask = readNIfTI(brain_mask_file, reorient=FALSE)

#EXTRACTION QUALITY CHECK

par(mfrow = c(1,3))
for (slice in c(5, 7, 8)) {
  image(
    pre_n4img[, , slice],
    col = gray(0:64 / 64),
    main = paste("Original + mask, slice:", slice),
    zlim = z_range,
    asp = 1
  )
  
  contour(
    brain_mask[, , slice],
    levels = 0.5,
    add = TRUE,
    drawlabels = FALSE,
    col = "red"
  )
}
par(mfrow = c(1, 1))
# Internal mask contours observed near the skull base.
# Possible local over or under-extraction from visual inspection.

#FETCHING FILES FOR T2
raw_dir_T2 = "~/Desktop/Neurohacking_data-master/BRAINIX/DICOM/T2"
all_slices_T2 = readDICOM(raw_dir_T2)

#CONVERSION OF T2 TO NIfTI + BASIC QUALITY CHECK
nii_T2=dicom2nifti(all_slices_T2)
par(mfrow = c(1, 3))
z_range_T2 <- range(nii_T2)

par(mfrow = c(1, 3))

image(
  nii_T2[,,6], 
  col = gray(0:64/64), 
  zlim=z_range_T2,
  asp = 1,
  main="Slice 6")
image(
  nii_T2[,,12], 
  col = gray(0:64/64), 
  zlim=z_range_T2,
  asp = 1,
  main="Slice 12")
image(
  nii_T2[,,18], 
  col = gray(0:64/64), 
  zlim=z_range_T2,
  asp = 1,
  main="Slice 18")

par(mfrow = c(1, 1))

#SAVING FILE AS NIfTI
writeNIfTI(nim = nii_T2, filename = file.path(processed_dir, "T2_converted"))

#N4 BIAS CORRECTION FOR T2
input_file_T2 = path.expand(file.path(processed_dir, "T2_converted.nii.gz"))
output_file_T2 = path.expand(file.path(processed_dir, "n4_T2_corrected.nii.gz"))

status = system2(
  command = n4_exe,
  args = c(
    "-d", "3",
    "-i", shQuote(input_file_T2),
    "-o", shQuote(output_file_T2),
    "-s", "2",
    "-v", "1"
  )
)
stopifnot(status == 0, file.exists(output_file_T2))
n4img_T2 = readNIfTI(output_file_T2, reorient=FALSE)
# N4 corrected T2 Pixel Spacing: 0.4492188x0.4492188x5mm.
# N4 corrected T2 Dimensions: 512x512 with 22 slices

par(mfrow = c(1,2))
pre_n4img_T2 = readNIfTI(input_file_T2, reorient=FALSE)
z_range_T2 <- range(pre_n4img_T2)
image(
  pre_n4img_T2[,,18], 
  col=gray(0:64/64), 
  zlim=z_range_T2,
  asp = 1,
  main="Original")
image(
  n4img_T2[,,18], 
  col=gray(0:64/64), 
  zlim=z_range_T2,
  asp = 1,
  main="N4 Corrected")
par(mfrow = c(1,1))

#FETCHING FILES FOR FLAIR
raw_dir_FLAIR = "~/Desktop/Neurohacking_data-master/BRAINIX/DICOM/FLAIR"
all_slices_FLAIR = readDICOM(raw_dir_FLAIR)

#CONVERSION OF FLAIR DICOM TO NIfTI + BASIC QUALITY CHECK
nii_FLAIR=dicom2nifti(all_slices_FLAIR)
# FLAIR Dimensions: 288x288 with 22 slices
# FLAIR Pixel spacing: 0.7986111x0.7986111x5mm

par(mfrow = c(1, 3))
z_range_FLAIR <- range(nii_FLAIR)

image(
  nii_FLAIR[,,6], 
  col = gray(0:64/64), 
  zlim=z_range_FLAIR,
  asp = 1,
  main="Slice 6")
image(
  nii_FLAIR[,,12], 
  col = gray(0:64/64), 
  zlim=z_range_FLAIR,
  asp = 1,
  main="Slice 12")
image(
  nii_FLAIR[,,18], 
  col = gray(0:64/64), 
  zlim=z_range_FLAIR,
  asp = 1,
  main="Slice 18")
par(mfrow = c(1, 1))

#SAVING FILE AS NIfTI
writeNIfTI(nim = nii_FLAIR, filename = file.path(processed_dir, "FLAIR_converted"))

#N4 BIAS CORRECTION FOR FLAIR

input_file_FLAIR = path.expand(file.path(processed_dir, "FLAIR_converted.nii.gz"))
output_file_FLAIR = path.expand(file.path(processed_dir, "n4_FLAIR_corrected.nii.gz"))

status = system2(
  command = n4_exe,
  args = c(
    "-d", "3",
    "-i", shQuote(input_file_FLAIR),
    "-o", shQuote(output_file_FLAIR),
    "-s", "2",
    "-v", "1"
  )
)
stopifnot(status == 0, file.exists(output_file_FLAIR))
n4img_FLAIR = readNIfTI(output_file_FLAIR, reorient=FALSE)

par(mfrow = c(1,2))
pre_n4img_FLAIR = readNIfTI(input_file_FLAIR, reorient=FALSE)
z_range_FLAIR <- range(pre_n4img_FLAIR)
image(
  pre_n4img_FLAIR[,,18], 
  col=gray(0:64/64), 
  zlim=z_range_FLAIR,
  asp = 1,
  main="Original")
image(
  n4img_FLAIR[,,18], 
  col=gray(0:64/64), 
  zlim=z_range_FLAIR,
  asp = 1,
  main="N4 Corrected")
par(mfrow = c(1,1))

# REGISTRATION

registration_prefix = path.expand(
  file.path(processed_dir, "T2_to_T1_")
)

stopifnot(
  file.exists(registration_exe),
  file.exists(output_file),
  file.exists(output_file_T2)
)

status = system2(
  command = registration_exe,
  args = c(
    "-d", "3",
    "-f", shQuote(output_file),
    "-m", shQuote(output_file_T2),
    "-t", "r",
    "-o", shQuote(registration_prefix)
  )
)

stopifnot(status == 0)

registered_T2_file = paste0(registration_prefix, "Warped.nii.gz")

stopifnot(file.exists(registered_T2_file))

registered_T2 = readNIfTI(
  registered_T2_file,
  reorient = FALSE
)
# REGISTRATION QUALITY CHECK
t1_reference = readNIfTI(output_file, reorient = FALSE)

par(mfrow = c(3, 2))

for (slice in c(6, 12, 18)) {
  image(
    t1_reference[, , slice],
    col = gray.colors(65),
    main = "Reference T1",
    axes = FALSE,
    asp = 1
  )
  
  image(
    registered_T2[, , slice],
    col = gray.colors(65),
    main = "T2 aligned to T1",
    axes = FALSE,
    asp = 1
  )
}

par(mfrow = c(1, 1))

# RIGID REGISTRATION: FLAIR TO T1

flair_registration_prefix = path.expand(
  file.path(processed_dir, "FLAIR_to_T1_")
)

status = system2(
  command = registration_exe,
  args = c(
    "-d", "3",
    "-f", shQuote(output_file),
    "-m", shQuote(output_file_FLAIR),
    "-t", "r",
    "-o", shQuote(flair_registration_prefix)
  )
)

registered_FLAIR_file = paste0(
  flair_registration_prefix, "Warped.nii.gz"
)

stopifnot(
  status == 0,
  file.exists(registered_FLAIR_file)
)

registered_FLAIR = readNIfTI(
  registered_FLAIR_file,
  reorient = FALSE
)


par(mfrow = c(3, 2), mar = c(1, 1, 3, 1))

for (slice in c(6, 12, 18)) {
  image(
    t1_reference[, , slice],
    col = gray.colors(65),
    main = "Reference T1",
    axes = FALSE,
    asp = 1
  )
  
  image(
    registered_FLAIR[, , slice],
    col = gray.colors(65),
    main = "FLAIR aligned to T1",
    axes = FALSE,
    asp = 1
  )
}

par(mfrow = c(1, 1))


#APPLICATION OF T1 BRAIN MASKS TO REGISTERED T2 AND FLAIR

image(
  registered_T2[, , round(dim(registered_T2)[3]*0.5)],
  col = gray(0:64 / 64),
  main = "Registered T2 with T1 brain mask",
  asp = 1
)

contour(
  brain_mask[, , round(dim(registered_T2)[3]*0.5)],
  levels = 0.5,
  add = TRUE,
  drawlabels = FALSE,
  col = "red"
)
stopifnot(
  identical(dim(brain_mask), dim(registered_T2)),
  identical(dim(brain_mask), dim(registered_FLAIR))
)
extracted_T2 = registered_T2 * brain_mask
extracted_FLAIR = registered_FLAIR * brain_mask

image(
  extracted_T2[, , round(dim(extracted_T2)[3]*0.5)],
  col = gray(0:64 / 64),
  main = "Extracted T2",
  asp = 1
)

writeNIfTI(nim = extracted_T2, filename = file.path(processed_dir, "T2_registered_brain"))

image(
  extracted_FLAIR[, , round(dim(extracted_FLAIR)[3]*0.5)],
  col = gray(0:64 / 64),
  main = "Extracted FLAIR",
  asp = 1
)
writeNIfTI(nim = extracted_FLAIR, filename = file.path(processed_dir, "FLAIR_registered_brain"))

stopifnot(
  file.exists(file.path(processed_dir, "T2_registered_brain.nii.gz")),
  file.exists(file.path(processed_dir, "FLAIR_registered_brain.nii.gz"))
)