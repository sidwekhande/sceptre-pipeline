#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
set.seed(4)

# obtain the command line arguments
sceptre_object_fp <- args[1]
response_odm_fp <- args[2]
grna_odm_fp <- args[3]
side <- args[4]
grna_integration_strategy <- args[5]
resampling_approximation <- args[6]
control_group <- args[7]
resampling_mechanism <- args[8]
multiple_testing_method <- args[9]
multiple_testing_alpha <- args[10]
formula_object_fp <- args[11]
discovery_pairs_fp <- args[12]
positive_control_pairs_fp <- args[13]
trial <- as.logical(args[14])
nuclear <- as.logical(args[15])

# load the sceptre object
sceptre_object <- sceptre::read_ondisc_backed_sceptre_object(sceptre_object_fp = sceptre_object_fp,
                                                             response_odm_file_fp = response_odm_fp,
                                                             grna_odm_file_fp = grna_odm_fp)
# process the default arguments
# side
if (identical(side, "default")) {
  side <- c("left", "both", "right")[sceptre_object@side_code + 2L]
}

#  grna_integration_strategy
if (identical(grna_integration_strategy, "default")) {
  grna_integration_strategy <- sceptre_object@grna_integration_strategy
}

# resampling_approximation
if (identical(resampling_approximation, "default")) {
  resampling_approximation <- sceptre_object@resampling_approximation
}

# control_group
if (identical(control_group, "default")) {
  control_group <- if (sceptre_object@control_group_complement) "complement" else "nt_cells"
}

# resampling_mechanism
if (identical(resampling_mechanism, "default")) {
  resampling_mechanism <- if (sceptre_object@run_permutations) "permutations" else "crt"
}

# formula_object
formula_object <- readRDS(formula_object_fp)
if (identical(formula_object, NULL)) {
  formula_object <- sceptre_object@formula_object
}

# discovery_pairs / positive_control_pairs
#
# Subset to just grna_target/response_id -- set_analysis_parameters()/run_qc() only ever use
# these two, and any extra decorative column (e.g. gene_symbol) that differs between
# discovery_pairs and positive_control_pairs makes sceptre's internal rbind() of the two crash
# with an unhelpful "numbers of columns of arguments do not match" deep inside
# compute_pairwise_qc_information, rather than a message naming the actual mismatch.
#
# discovery_pairs/positive_control_pairs are TSVs; an empty (header-only) TSV means "not
# supplied, fall back to what's already on the sceptre object". "nuclear" (trans/genome-wide
# mode) is now passed in explicitly rather than sniffed from the file content, since a TSV can't
# carry the old RDS placeholder's "the whole value is the string 'trans'" sentinel.
required_pair_cols <- c("grna_target", "response_id")
read_pairs_tsv <- function(fp) data.table::fread(fp, colClasses = "character") |> as.data.frame()

discovery_pairs <- read_pairs_tsv(discovery_pairs_fp)
if (nuclear) {
  discovery_pairs <- data.frame(grna_target = character(0), response_id = character(0))
} else {
  if (nrow(discovery_pairs) == 0) {
    discovery_pairs <- sceptre_object@discovery_pairs
  }
  discovery_pairs <- discovery_pairs[, required_pair_cols]
  # discovery_pairs commonly names more candidate elements than ended up with guides in the
  # final library (e.g. dropped during synthesis/QC); restrict to elements sceptre actually knows
  discovery_pairs <- discovery_pairs[which(discovery_pairs$grna_target %in% sceptre_object@grna_target_data_frame$grna_target), ]
  if (trial) {
    n_pairs <- nrow(discovery_pairs)
    discovery_pairs <- discovery_pairs |> dplyr::sample_n(min(100, n_pairs))
  }
}

# positive_control_pairs
if (nuclear) {
  positive_control_pairs <- data.frame(grna_target = character(0), response_id = character(0))
} else {
  positive_control_pairs <- read_pairs_tsv(positive_control_pairs_fp)
  if (nrow(positive_control_pairs) == 0) {
    positive_control_pairs <- sceptre_object@positive_control_pairs
  }
  positive_control_pairs <- positive_control_pairs[, required_pair_cols]
  positive_control_pairs <- positive_control_pairs[which(positive_control_pairs$grna_target %in% sceptre_object@grna_target_data_frame$grna_target), ]
  if (trial) {
    n_pairs <- nrow(positive_control_pairs)
    positive_control_pairs <- positive_control_pairs |> dplyr::sample_n(min(100, n_pairs))
  }
}

# multiple_testing_method
if (nuclear) {
  multiple_testing_method <- "none"
} else if (identical(multiple_testing_method, "default")) {
  multiple_testing_method <- sceptre_object@multiple_testing_method
}

# multiple_testing_alpha
if (identical(multiple_testing_alpha, "default")) {
  multiple_testing_alpha <- sceptre_object@multiple_testing_alpha
} else {
  multiple_testing_alpha <- as.numeric(multiple_testing_alpha)
}

# set the analysis parameters
sceptre_object <- sceptre::set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs,
  positive_control_pairs = positive_control_pairs,
  side = side,
  grna_integration_strategy = grna_integration_strategy,
  formula_object = formula_object,
  resampling_approximation = resampling_approximation,
  control_group = control_group,
  resampling_mechanism = resampling_mechanism,
  multiple_testing_method = multiple_testing_method,
  multiple_testing_alpha = multiple_testing_alpha
)
if (nuclear) sceptre_object@nuclear <- TRUE

# write the sceptre_object
saveRDS(sceptre_object, "sceptre_object.rds")
sink(file = "analysis_summary.txt", append = FALSE)
sceptre::print(sceptre_object)
sink(NULL)
