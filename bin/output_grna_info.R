#!/usr/bin/env Rscript

# obtain the command line arguments
args <- commandArgs(trailingOnly = TRUE)
sceptre_object_fp <- args[1]
response_odm_fp <- args[2]
grna_odm_fp <- args[3]
grna_pod_size <- as.integer(args[4])
trial <- as.logical(args[5])
grna_assignment_method <- args[6]
threshold <- args[7]
umi_fraction_threshold <- args[8]
min_grna_n_umis_threshold <- args[9]
n_em_rep <- args[10]
n_nonzero_cells_cutoff <- args[11]
backup_threshold <- args[12]
probability_threshold <- args[13]
grna_assignment_formula_supplied <- as.logical(args[14])

# load the sceptre object
sceptre_object <- sceptre::read_ondisc_backed_sceptre_object(sceptre_object_fp = sceptre_object_fp,
                                                             response_odm_file_fp = response_odm_fp,
                                                             grna_odm_file_fp = grna_odm_fp)
# get the gRNAs in use and obtain the gRNA-to-pod map
grnas_in_use <- sceptre:::determine_grnas_in_use(sceptre_object, trial)
grna_to_pod_map <- data.frame(grna_id = grnas_in_use,
                              pod_id = sceptre:::get_id_vect(grnas_in_use, grna_pod_size))
grna_pods <- unique(grna_to_pod_map$pod_id)

# resolve the gRNA assignment method and its arguments once, before the
# assignment jobs are scattered across pods
if (identical(grna_assignment_method, "default")) {
  grna_assignment_method <- if (sceptre_object@low_moi) "maximum" else "mixture"
}
grna_assignment_arg_values <- list(
  threshold = threshold,
  umi_fraction_threshold = umi_fraction_threshold,
  min_grna_n_umis_threshold = min_grna_n_umis_threshold,
  n_em_rep = n_em_rep,
  n_nonzero_cells_cutoff = n_nonzero_cells_cutoff,
  backup_threshold = backup_threshold,
  probability_threshold = probability_threshold
)
supplied_arg_names <- names(grna_assignment_arg_values)[
  !vapply(grna_assignment_arg_values, identical, logical(1), "default")
]
if (grna_assignment_formula_supplied) {
  supplied_arg_names <- c(supplied_arg_names, "grna_assignment_formula")
}
applicable_arg_names <- switch(
  grna_assignment_method,
  thresholding = "threshold",
  mixture = c("n_em_rep", "n_nonzero_cells_cutoff", "backup_threshold",
              "probability_threshold", "grna_assignment_formula"),
  maximum = c("umi_fraction_threshold", "min_grna_n_umis_threshold"),
  stop("Unrecognized gRNA assignment method: ", grna_assignment_method)
)
ignored_arg_names <- setdiff(supplied_arg_names, applicable_arg_names)
active_arg_names <- intersect(supplied_arg_names, applicable_arg_names)
grna_assignment_formula_in_use <-
  "grna_assignment_formula" %in% active_arg_names
grna_assignment_args <- list(method = grna_assignment_method)
for (arg_name in setdiff(active_arg_names, "grna_assignment_formula")) {
  grna_assignment_args[[arg_name]] <-
    as.numeric(grna_assignment_arg_values[[arg_name]])
}

# write and save outputs
sceptre:::write_vector(grna_pods, "grna_pods.txt")
saveRDS(grna_to_pod_map, "grna_to_pod_map.rds")
saveRDS(grna_assignment_args, "grna_assignment_args.rds")
sceptre:::write_vector(grna_assignment_method, "grna_assignment_method.txt")
sceptre:::write_vector(
  tolower(grna_assignment_formula_in_use),
  "grna_assignment_formula_in_use.txt"
)
warning_message <- if (length(ignored_arg_names) > 0L) {
  paste0(
    "Ignoring ", paste0("--", ignored_arg_names, collapse = ", "),
    " because ",
    if (length(ignored_arg_names) == 1L) "it does" else "they do",
    " not apply to the ", grna_assignment_method,
    " gRNA assignment method."
  )
} else {
  character()
}
writeLines(warning_message, "grna_assignment_warning.txt")
