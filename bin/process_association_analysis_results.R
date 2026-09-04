#!/usr/bin/env Rscript

# 1. obtain the command line arguments
args <- commandArgs(trailingOnly = TRUE)
sceptre_object_fp <- args[1]
response_odm_fp <- args[2]
grna_odm_fp <- args[3]
analysis_type <- args[4]
result_fps <- grep(x = args, pattern = "results*", value = TRUE)
precomp_fps <- grep(x = args, pattern = "precomputations*", value = TRUE)

# 2. load the sceptre object; update completed analysis flags
sceptre_object <- sceptre::read_ondisc_backed_sceptre_object(sceptre_object_fp = sceptre_object_fp,
                                                             response_odm_file_fp = response_odm_fp,
                                                             grna_odm_file_fp = grna_odm_fp)
sceptre_object@functs_called[[analysis_type]] <- TRUE
sceptre_object@last_function_called <- analysis_type

# 3. prune the pairs_with_info fields of the sceptre_object
#
# The discovery-analysis branch no longer prunes @discovery_pairs_with_info or @covariate_matrix
# here, since downstream power-analysis tooling reads them directly -- see the PR description.
if (analysis_type == "run_calibration_check") {
  sceptre_object@negative_control_pairs <- data.frame()
} else if (analysis_type == "run_power_check") {
  sceptre_object@positive_control_pairs_with_info <- data.frame()
}
gc() |> invisible()

# 4. process the results
result_df <- lapply(result_fps, readRDS) |> data.table::rbindlist()
data.table::setorderv(result_df, cols = c("p_value", "response_id"), na.last = TRUE)
process_funct <- switch(analysis_type,
                        run_calibration_check = sceptre:::process_calibration_result,
                        run_power_check = sceptre:::apply_grouping_to_result,
                        run_discovery_analysis = sceptre:::process_discovery_result)
result_df <- process_funct(result_df, sceptre_object)
result_df$pod <- NULL

# 5. add results to the sceptre_object
if (analysis_type == "run_calibration_check") {
  sceptre_object@calibration_result <- result_df
} else if (analysis_type == "run_power_check") {
  sceptre_object@power_result <- result_df
} else { # discovery analysis
  sceptre_object@discovery_result <- result_df
}

# 6. process the response precomputations (if using the complement set as the control group, and not running a discovery analysis)
if (sceptre_object@control_group_complement && analysis_type != "run_discovery_analysis") {
  precomputation_list <- lapply(precomp_fps, readRDS) |> unlist(recursive = FALSE)
  sceptre_object@response_precomputations <- c(sceptre_object@response_precomputations, precomputation_list)
}

# 6. create the plot
p <- sceptre::plot(sceptre_object)

# 7. delete unecessary info from sceptre_object
if (analysis_type == "run_calibration_check") {
  sceptre_object@calibration_result <- result_df[,c("p_value", "log_2_fold_change", "significant")]
} else if (analysis_type == "run_power_check") {
  # pass
} else { # discovery analysis
  # pass
}

# save the outputs
#
# The discovery-analysis branch used to write saveRDS(NULL, ...) here instead -- see the PR
# description for why that's wrong and what it broke.
saveRDS(sceptre_object, "sceptre_object.rds")
saveRDS(result_df, paste0("results_", analysis_type, ".rds"))
ggplot2::ggsave(filename = paste0("plot_", analysis_type, ".png"), plot = p,
                device = "png", scale = 1.1, width = 5, height = 4, dpi = 330)
sink(file = "analysis_summary.txt", append = FALSE)
sceptre::print(sceptre_object)
sink(NULL)
