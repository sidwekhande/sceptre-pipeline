#!/usr/bin/env Rscript

# obtain the command line arguments
args <- commandArgs(trailingOnly = TRUE)
sceptre_object_fp <- args[1]
response_odm_fp <- args[2]
grna_odm_fp <- args[3]
grna_to_pod_map_fp <- args[4]
grna_pod <- as.integer(args[5])
grna_assignment_args_fp <- args[6]
grna_assignment_formula_fp <- args[7]

# load the sceptre object
sceptre_object <- sceptre::read_ondisc_backed_sceptre_object(sceptre_object_fp = sceptre_object_fp,
                                                             response_odm_file_fp = response_odm_fp,
                                                             grna_odm_file_fp = grna_odm_fp)

# remove disc pairs
sceptre_object@discovery_pairs <- data.frame()
gc() |> invisible()

# load the grna to pod map; determine the grnas in use
grna_to_pod_map <- readRDS(grna_to_pod_map_fp)
grnas_in_use <- subset(grna_to_pod_map, pod_id == grna_pod)$grna_id
sceptre_object@elements_to_analyze <- grnas_in_use
sceptre_object@nf_pipeline <- TRUE

# load the arguments resolved by the one-time preparation step; Nextflow stages
# the user-supplied formula only when that step determines it is applicable
grna_assignment_args <- readRDS(grna_assignment_args_fp)
grna_assignment_formula <- readRDS(grna_assignment_formula_fp)
if (!identical(grna_assignment_formula, NULL)) {
  grna_assignment_args$formula_object <- grna_assignment_formula
}
args_to_pass <- c(
  list(sceptre_object = sceptre_object),
  grna_assignment_args,
  list(parallel = FALSE)
)

# call the gRNA-to-cell assignment function
sceptre_object <- do.call(what = sceptre::assign_grnas, args = args_to_pass)

# save the initial assignment list
grna_assignment_formula <- sceptre_object@grna_assignment_hyperparameters$formula_object 
saveRDS(sceptre_object@initial_grna_assignment_list, "grna_assignments.rds")
saveRDS(grna_assignment_formula, "grna_assignment_formula.rds")
