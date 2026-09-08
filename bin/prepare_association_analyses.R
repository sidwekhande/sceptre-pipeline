#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
sceptre_object_fp <- args[1]
response_odm_fp <- args[2]
grna_odm_fp <- args[3]
n_calibration_pairs <- args[4]
calibration_group_size <- args[5]
pair_pod_size <- as.integer(args[6])

# load the sceptre object
sceptre_object <- sceptre::read_ondisc_backed_sceptre_object(sceptre_object_fp = sceptre_object_fp,
                                                             response_odm_file_fp = response_odm_fp,
                                                             grna_odm_file_fp = grna_odm_fp)

# handle the default arguments
args_to_pass <- list(sceptre_object = sceptre_object)
optional_args_names <- c("n_calibration_pairs", "calibration_group_size")
for (optional_arg_name in optional_args_names) {
  optional_arg_value <- get(x = optional_arg_name)
  if (!identical(optional_arg_value, "default")) {
    args_to_pass[[optional_arg_name]] <- as.numeric(optional_arg_value)
  }
}

# generate the negative control pairs
if (n_calibration_pairs >= 1L) sceptre_object <- do.call(sceptre:::run_calibration_check_pt_1, args = args_to_pass)

# assign pods to negative control, positive control, and discovery pairs
data_table_list <- list(calibration_check = sceptre_object@negative_control_pairs,
                        power_check = sceptre_object@positive_control_pairs_with_info,
                        discovery_analysis = sceptre_object@discovery_pairs_with_info)

# process each of the data tables
#
# Returned as a plain data.frame, not a data.table: sceptre_object's S4 slots for these fields
# validate against data.frame, and while data.table inherits it (is(x, "data.frame") is TRUE),
# the S4 slot-assignment validity check does not accept it -- confirmed by WattEG's
# fit_null_models.R, the first caller to reassign one of these into a fresh sceptre_object, which
# fails with "assignment of an object of class \"data.table\" is not valid for @'...'". Coercing
# here, once, is cheaper than chasing the same failure at every future call site.
process_pair_data_table <- function(data_table) {
  if (nrow(data_table) >= 1L) {
    data_table_pass_qc <- data_table[data_table$pass_qc,]
    data.table::setorderv(data_table_pass_qc, "response_id")
    data_table_pass_qc$pod <- sceptre:::get_id_vect(v = data_table_pass_qc$response_id, pod_size = pair_pod_size)
    data_table_fail_qc <- data_table[!data_table$pass_qc,] |> dplyr::mutate(pod = 1L)
    data_table <- data.table::rbindlist(list(data_table_pass_qc, data_table_fail_qc))
  }
  return(as.data.frame(data_table))
}
processed_data_table_list <- lapply(data_table_list, process_pair_data_table)
sceptre_object@negative_control_pairs <- processed_data_table_list$calibration_check
sceptre_object@positive_control_pairs_with_info <- processed_data_table_list$power_check
sceptre_object@discovery_pairs_with_info <- processed_data_table_list$discovery_analysis

# remove fields no longer needed
sceptre_object@M_matrix <- matrix()

# save the output
for (analysis_name in names(processed_data_table_list)) {
  data_table <- processed_data_table_list[[analysis_name]]
  run_analysis <- tolower(nrow(data_table) >= 1L)
  sceptre:::write_vector(run_analysis, paste0("run_", analysis_name))
  sceptre:::write_vector(unique(data_table$pod), paste0(analysis_name, "_pods"))
}
saveRDS(sceptre_object, "sceptre_object.rds")
