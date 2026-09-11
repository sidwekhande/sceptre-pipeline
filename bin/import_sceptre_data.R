#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
set.seed(4)

# obtain the command line arguments
trial <- as.logical(args[1])
moi <- args[2]
grna_target_tsv_fp <- args[3]
directories <- args[-(1:3)]

# grna_target_tsv is already the finished data frame (targeting + non-targeting gRNAs together,
# non-targeting rows have grna_target = "non-targeting"); no processing done on it here, that
# prep happens upstream of this pipeline
grna_target_data_frame <- data.table::fread(grna_target_tsv_fp, colClasses = "character") |> as.data.frame()

# import from Cell Ranger output, always odm-backed -- downstream processes in this
# pipeline only support odm-backed sceptre objects
sceptre_object <- sceptre::import_data_from_cellranger(
  directories = directories,
  moi = moi,
  grna_target_data_frame = grna_target_data_frame,
  use_ondisc = TRUE,
  directory_to_write = "."
)

# gene.odm / grna.odm were already written to "." by import_data_from_cellranger() above;
# this writes sceptre_object.rds with a matching @integer_id
sceptre::write_ondisc_backed_sceptre_object(
  sceptre_object = sceptre_object,
  directory_to_write = "."
)
