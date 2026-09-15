#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
set.seed(4)

# obtain the command line arguments
trial <- as.logical(args[1])
moi <- args[2]
per_guide_metadata_tsv_fp <- args[3]
directories <- args[-(1:3)]

# per_guide_metadata_tsv follows the IGVF per-guide metadata submission format (one row per
# guide_id; targeting == "False" for non-targeting/control guides; intended_target_name groups
# guides into elements for element-level analysis -- this becomes sceptre's grna_target)
guide_metadata <- data.table::fread(per_guide_metadata_tsv_fp, colClasses = "character") |> as.data.frame()

grna_target_data_frame <- data.frame(
  grna_id = guide_metadata$guide_id,
  grna_target = ifelse(guide_metadata$targeting == "True", guide_metadata$intended_target_name, "non-targeting")
)

# positive controls: per the spec, if genomic_element == "promoter" the gene's own ENSEMBL id is
# already used as intended_target_name, so grna_target and response_id coincide; if a future
# dataset instead populates putative_target_genes (required for enhancer/insulator/silencer/
# distal-element positive controls), that column should be exploded into one row per gene here too
positive_controls <- guide_metadata[guide_metadata$type == "positive control", ]
positive_control_pairs <- data.frame(
  grna_target = positive_controls$intended_target_name,
  response_id = positive_controls$intended_target_name
) |> dplyr::distinct()
data.table::fwrite(positive_control_pairs, "positive_control_pairs.tsv", sep = "\t")

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
