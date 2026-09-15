#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
set.seed(4)

# obtain the command line arguments
trial <- as.logical(args[1])
moi <- args[2]
per_guide_metadata_tsv_fp <- args[3]
directories <- args[-(1:3)]

# per_guide_metadata_tsv follows the IGVF per-guide metadata submission format (one row per
# guide_id; targeting == "False" for non-targeting/control guides). grna_target is the element's
# own coordinates (intended_target_chr:intended_target_start-intended_target_end) uniformly for
# every targeting guide, including positive controls -- intended_target_name is NOT used here
# since, per spec, it's the regulated gene's ENSEMBL id for promoter-targeting guides rather than
# the element's own locus, which would make positive-control grna_targets inconsistent with the
# element-based grouping every other targeting guide uses.
guide_metadata <- data.table::fread(per_guide_metadata_tsv_fp, colClasses = "character") |> as.data.frame()

grna_target_data_frame <- data.frame(
  grna_id = guide_metadata$guide_id,
  grna_target = ifelse(
    guide_metadata$targeting == "True",
    paste0(guide_metadata$intended_target_chr, ":",
           as.integer(as.numeric(guide_metadata$intended_target_start)), "-",
           as.integer(as.numeric(guide_metadata$intended_target_end))),
    "non-targeting"
  )
)

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
