#!/usr/bin/env Rscript

# obtain the command line arguments
args <- commandArgs(trailingOnly = TRUE)
sceptre_object_fp <- args[1]
response_odm_fp <- args[2]
grna_odm_fp <- args[3]
n_nonzero_trt_thresh <- args[4]
n_nonzero_cntrl_thresh <- args[5]
response_n_umis_range_lower <- args[6]
response_n_umis_range_upper <- args[7]
response_n_nonzero_range_lower <- args[8]
response_n_nonzero_range_upper <- args[9]
p_mito_threshold <- args[10]

# load the sceptre object
sceptre_object <- sceptre::read_ondisc_backed_sceptre_object(sceptre_object_fp = sceptre_object_fp,
                                                             response_odm_file_fp = response_odm_fp,
                                                             grna_odm_file_fp = grna_odm_fp)

# process n_nonzero_trt_thresh, n_nonzero_cntrl_thresh, and p_mito_threshold
args_to_pass <- list(sceptre_object)
optional_args_names <- c("n_nonzero_trt_thresh", "n_nonzero_cntrl_thresh", "p_mito_threshold")
for (optional_arg_name in optional_args_names) {
  optional_arg_value <- get(x = optional_arg_name)
  if (!identical(optional_arg_value, "default")) {
    args_to_pass[[optional_arg_name]] <- as.numeric(optional_arg_value)
  }
}

# process response_n_umis_range and response_n_nonzero_range
response_n_umis_range <- c(0.01, 0.99)
response_n_nonzero_range <- c(0.01, 0.99)
if (!identical(response_n_umis_range_lower, "default")) response_n_umis_range[1] <- as.numeric(response_n_umis_range_lower)
if (!identical(response_n_umis_range_upper, "default")) response_n_umis_range[2] <- as.numeric(response_n_umis_range_upper)
if (!identical(response_n_nonzero_range_lower, "default")) response_n_nonzero_range[1] <- as.numeric(response_n_nonzero_range_lower)
if (!identical(response_n_nonzero_range_upper, "default")) response_n_nonzero_range[2] <- as.numeric(response_n_nonzero_range_upper)
args_to_pass[["response_n_umis_range"]] <- response_n_umis_range
args_to_pass[["response_n_nonzero_range"]] <- response_n_nonzero_range

# call the qc function
gc() |> invisible()
sceptre_object <- do.call(what = sceptre::run_qc, args = args_to_pass)

# create plots
p1 <- sceptre::plot_covariates(sceptre_object)
p2 <- sceptre::plot_run_qc(sceptre_object)

# remove fields no longer needed by this pipeline's own downstream steps.
#
# @initial_grna_assignment_list and @covariate_data_frame are no longer pruned here, since
# downstream power-analysis tooling reads them directly -- see the PR description for details.
sceptre_object@discovery_pairs <- data.frame()
sceptre_object@positive_control_pairs <- data.frame()
sceptre_object@grna_assignments_raw <- list()
sceptre_object@import_grna_assignment_info <- list()

# write outputs to disk
saveRDS(sceptre_object, "sceptre_object.rds")
ggplot2::ggsave(filename = "plot_covariates.png", plot = p1, device = "png",
                scale = 1.1, width = 5, height = 4, dpi = 330)
ggplot2::ggsave(filename = "plot_run_qc.png", plot = p2, device = "png",
                scale = 1.1, width = 5, height = 4, dpi = 330)
sink(file = "analysis_summary.txt", append = FALSE)
sceptre::print(sceptre_object)
sink(NULL)
