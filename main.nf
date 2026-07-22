// use DSL 2
nextflow.enable.dsl = 2

/*************************
* DEFAULT PARAMETER VALUES
*************************/
// 1. pipeline meta params
params.pipeline_stop = "run_discovery_analysis"
params.trial = "false"
// 2. set analysis parameters
params.side = "default"
params.grna_integration_strategy = "default"
params.resampling_approximation = "default"
params.control_group = "default"
params.resampling_mechanism = "default"
params.multiple_testing_method = "default"
params.multiple_testing_alpha = "default"
params.formula_object = "${baseDir}/resources/placeholder_file.rds"
params.discovery_pairs = "${baseDir}/resources/placeholder_file.rds"
params.positive_control_pairs = "${baseDir}/resources/placeholder_file.rds"
// 3. gRNA assignment
params.grna_assignment_method = "default"
params.threshold = "default"
params.umi_fraction_threshold = "default"
params.min_grna_n_umis_threshold = "default"
params.n_em_rep = "default"
params.n_nonzero_cells_cutoff = "default"
params.backup_threshold = "default"
params.probability_threshold = "default"
params.grna_assignment_formula = "${baseDir}/resources/placeholder_file.rds"
// 4. QC
params.n_nonzero_trt_thresh = "default"
params.n_nonzero_cntrl_thresh = "default"
params.response_n_umis_range_lower = "default"
params.response_n_umis_range_upper = "default"
params.response_n_nonzero_range_lower = "default"
params.response_n_nonzero_range_upper = "default"
params.p_mito_threshold = "default"
// 5. calibration check
params.n_calibration_pairs = "default"
params.calibration_group_size = "default"
// 6. computation: parallelization
params.grna_pod_size = 150
if ("$params.discovery_pairs" == "trans") params.pair_pod_size = 500000 else params.pair_pod_size = 25000
// 7. computation: time
params.set_analysis_parameters_time = "15m" // set analysis parameters
params.prepare_assign_grnas_time = "15m" // prepare grna assignments
params.assign_grnas_time_per_grna = "2s" // assign grnas
params.combine_assign_grnas_time = "15m" // process grna assignments
params.run_qc_time = "60m" // run qc
params.prepare_association_analysis_time = "15m" // prepare association analyses
params.run_association_analysis_time_per_pair = "0.05s" // run association analysis
params.combine_association_analysis_time = "15m" // process association analysis
// 8. computation: memory
params.set_analysis_parameters_memory = "4GB" // set analysis parameters
params.prepare_assign_grnas_memory = "4GB" // prepare grna assignments
params.assign_grnas_memory = "4GB" // assign grnas
params.combine_assign_grnas_memory = "4GB"  // process grna assignments
if ("$params.discovery_pairs" == "trans") params.run_qc_memory = "4GB" else params.run_qc_memory = "8GB" // run qc
params.prepare_association_analysis_memory = "4GB" // prepare association analyses
params.run_association_analysis_memory = "4GB" // run association analysis
params.combine_association_analysis_memory = "4GB" // process association analysis
// 9. output file type
params.use_parquet = "true" // when running massive-scale "trans" analysis, save results as .parquet (true) or .rds (false)

/********************
* INCLUDE SUBWORKFLOW
*********************/
include { run_analysis_subworkflow as run_analysis_subworkflow_calibration_check } from './run_analysis_subworkflow.nf'
include { run_analysis_subworkflow as run_analysis_subworkflow_power_check } from './run_analysis_subworkflow.nf'
include { run_analysis_subworkflow as run_analysis_subworkflow_discovery_analysis } from './run_analysis_subworkflow.nf'

/*****************************
* GROOVY PROCESSING OF INPUTS
*****************************/
pipeline_steps = ["set_analysis_parameters", "assign_grnas", "run_qc", "run_calibration_check", "run_power_check", "run_discovery_analysis"]
// get the rank of the input params.pipeline_stop; throw an error if not present in the list
def step_rank = pipeline_steps.indexOf(params.pipeline_stop)
if (step_rank == -1) {
    throw new Exception("'$params.pipeline_stop' is not a step of the sceptre pipeline. The parameter 'pipeline_stop' should be set to one of 'assign_grnas', 'run_qc', 'run_calibration_check', 'run_power_check', or 'run_discovery_analysis'.")
}
if ("$params.trial" == "true") {
  println "Running pipeline in trial mode."
}
disc_pairs = params.discovery_pairs
nuclear = false
if ("$params.discovery_pairs" == "trans") {
  disc_pairs = "${baseDir}/resources/trans_placeholder.rds"
  nuclear = true
}

/**********
* PROCESSES
**********/
// PROCESS A: set analysis parameters
process set_analysis_parameters {
  publishDir "${params.output_directory}", mode: 'copy', overwrite: true, pattern: "*.txt"

  time params.set_analysis_parameters_time
  memory params.set_analysis_parameters_memory

  input:
  path "sceptre_object_fp"
  path "response_odm_fp"
  path "grna_odm_fp"
  path "formula_object"
  path "discovery_pairs"
  path "positive_control_pairs"

  output:
  path "sceptre_object.rds", emit: sceptre_object_ch
  path "analysis_summary.txt"

  """
  set_analysis_parameters.R $sceptre_object_fp \
  $response_odm_fp \
  $grna_odm_fp \
  ${params.side} \
  ${params.grna_integration_strategy} \
  ${params.resampling_approximation} \
  ${params.control_group} \
  ${params.resampling_mechanism} \
  ${params.multiple_testing_method} \
  ${params.multiple_testing_alpha} \
  $formula_object \
  $discovery_pairs \
  $positive_control_pairs \
  ${params.trial}
  """
}

// PROCESS B: output gRNA info
process prepare_assign_grnas {
  time params.prepare_assign_grnas_time
  memory params.prepare_assign_grnas_memory

  input:
  path "sceptre_object_fp"
  path "response_odm_fp"
  path "grna_odm_fp"
  val "grna_assignment_formula_supplied"

  output:
  path "grna_to_pod_map.rds", emit: grna_to_pod_map_ch
  path "grna_pods.txt", emit: grna_pods_ch
  path "grna_assignment_args.rds", emit: grna_assignment_args_ch
  path "grna_assignment_method.txt", emit: grna_assignment_method_ch
  path "grna_assignment_formula_in_use.txt", emit: grna_assignment_formula_in_use_ch
  path "grna_assignment_warning.txt", emit: grna_assignment_warning_ch

  """
  output_grna_info.R $sceptre_object_fp \
  $response_odm_fp \
  $grna_odm_fp \
  ${params.grna_pod_size} \
  ${params.trial} \
  ${params.grna_assignment_method} \
  ${params.threshold} \
  ${params.umi_fraction_threshold} \
  ${params.min_grna_n_umis_threshold} \
  ${params.n_em_rep} \
  ${params.n_nonzero_cells_cutoff} \
  ${params.backup_threshold} \
  ${params.probability_threshold} \
  $grna_assignment_formula_supplied
  """
}

// PROCESS C: assign gRNAs
process assign_grnas {
  time {params.assign_grnas_time_per_grna * params.grna_pod_size}
  memory params.assign_grnas_memory

  when:
  grna_assignment_method != "maximum"
  
  input:
  path "sceptre_object_fp"
  path "response_odm_fp"
  path "grna_odm_fp"
  path "grna_to_pod_map"
  val "grna_pod"
  val "grna_assignment_method"
  path "grna_assignment_args"
  path "grna_assignment_formula"

  output:
  path "grna_assignments.rds", emit: grna_assignments_ch
  path "grna_assignment_formula.rds", emit: grna_assignment_formula_ch

  """
  assign_grnas.R $sceptre_object_fp \
  $response_odm_fp \
  $grna_odm_fp \
  $grna_to_pod_map \
  $grna_pod \
  $grna_assignment_args \
  $grna_assignment_formula
  """
}

// PROCESS C.1 assign_grnas_dummy
process dummy_assign_grnas {
  time "1m"
  memory "1GB"
  
  when:
  grna_assignment_method == "maximum"

  input:
  val "grna_assignment_method"

  output:
  path "grna_assignments.rds", emit: grna_assignments_ch
  path "grna_assignment_formula.rds", emit: grna_assignment_formula_ch

  """
  touch grna_assignments.rds
  touch grna_assignment_formula.rds
  """
}


// PROCESS D: process gRNA assignments
process combine_assign_grnas {
  publishDir "${params.output_directory}", mode: 'copy', overwrite: true, pattern: "*.png"
  publishDir "${params.output_directory}", mode: 'copy', overwrite: true, pattern: "*.txt"
  publishDir "${params.output_directory}", mode: 'copy', overwrite: true, pattern: "grna_assignment_matrix.rds"

  time params.combine_assign_grnas_time
  memory params.combine_assign_grnas_memory

  input:
  path "sceptre_object_fp"
  path "response_odm_fp"
  path "grna_odm_fp"
  path "grna_assignment_args"
  path "grna_assignment_formula"
  path "grna_assignments"

  output:
  path "plot_grna_count_distributions.png"
  path "plot_assign_grnas.png"
  path "analysis_summary.txt"
  path "grna_assignment_matrix.rds"
  path "sceptre_object.rds", emit: sceptre_object_ch

  """
  process_grna_assignments.R $sceptre_object_fp \
  $response_odm_fp \
  $grna_odm_fp \
  $grna_assignment_args \
  $grna_assignment_formula \
  grna_assignments*
  """
}

// PROCESS E.1: quality control
process run_qc {
  publishDir "${params.output_directory}", mode: 'copy', overwrite: true, pattern: "*.png"
  publishDir "${params.output_directory}", mode: 'copy', overwrite: true, pattern: "*.txt"

  time params.run_qc_time
  memory params.run_qc_memory

  input:
  path "sceptre_object_fp"
  path "response_odm_fp"
  path "grna_odm_fp"

  output:
  path "plot_covariates.png"
  path "plot_run_qc.png"
  path "analysis_summary.txt"
  path "sceptre_object.rds", emit: sceptre_object_ch

  """
  run_qc.R $sceptre_object_fp \
  $response_odm_fp \
  $grna_odm_fp \
  ${params.n_nonzero_trt_thresh} \
  ${params.n_nonzero_cntrl_thresh} \
  ${params.response_n_umis_range_lower} \
  ${params.response_n_umis_range_upper} \
  ${params.response_n_nonzero_range_lower} \
  ${params.response_n_nonzero_range_upper} \
  ${params.p_mito_threshold}
  """
}

// PROCESS F.1: prepare association analysis
process prepare_association_analysis {
  time params.prepare_association_analysis_time
  memory params.prepare_association_analysis_memory

  input:
  path "sceptre_object_fp"
  path "response_odm_fp"
  path "grna_odm_fp"

  output:
  path "run_calibration_check", emit: run_calibration_check_ch
  path "run_discovery_analysis", emit: run_discovery_analysis_ch
  path "run_power_check", emit: run_power_check_ch
  path "calibration_check_pods", emit: calibration_check_pods_ch
  path "power_check_pods", emit: power_check_pods_ch
  path "discovery_analysis_pods", emit: discovery_analysis_pods_ch
  path "sceptre_object.rds", emit: sceptre_object_ch

  """
  prepare_association_analyses.R $sceptre_object_fp \
  $response_odm_fp \
  $grna_odm_fp \
  ${params.n_calibration_pairs} \
  ${params.calibration_group_size} \
  ${params.pair_pod_size}
  """
}

// PROCESS E.2: quality control trans
process run_qc_trans {
  publishDir "${params.output_directory}", mode: 'copy', overwrite: true, pattern: "*.png"
  publishDir "${params.output_directory}", mode: 'copy', overwrite: true, pattern: "*.txt"

  time params.run_qc_time
  memory params.run_qc_memory

  input:
  path "sceptre_object_fp"
  path "response_odm_fp"
  path "grna_odm_fp"

  output:
  path "plot_covariates.png"
  path "plot_cellwise_qc.png"
  path "analysis_summary.txt"
  path "sceptre_object.rds", emit: sceptre_object_ch

  """
  run_qc_nuclear.R $sceptre_object_fp \
  $response_odm_fp \
  $grna_odm_fp \
  ${params.response_n_umis_range_lower} \
  ${params.response_n_umis_range_upper} \
  ${params.response_n_nonzero_range_lower} \
  ${params.response_n_nonzero_range_upper} \
  ${params.p_mito_threshold}
  """
}

// PROCESS F.2: prepare association analysis trans
process prepare_association_analysis_trans {
  publishDir "${params.output_directory}", mode: 'copy', overwrite: true, pattern: "response_to_pod_map.rds"
  
  time params.prepare_association_analysis_time
  memory params.prepare_association_analysis_memory

  input:
  path "sceptre_object_fp"
  path "response_odm_fp"
  path "grna_odm_fp"

  output:
  path "discovery_analysis_pods", emit: discovery_analysis_pods_ch
  path "response_to_pod_map.rds", emit: response_to_pod_map_ch

  """
  prepare_association_analysis_nuclear.R $sceptre_object_fp \
  $response_odm_fp \
  $grna_odm_fp \
  ${params.pair_pod_size}
  """
}

// PROCESS G: association analysis trans
process run_discovery_analysis_trans {
  publishDir "${params.output_directory}/trans_results", mode: 'copy', overwrite: true, pattern: "result_${pair_pod}*"
  memory params.run_association_analysis_memory
  time {params.run_association_analysis_time_per_pair * params.pair_pod_size}

  input:
  path "sceptre_object_fp"
  path "response_odm_fp"
  path "grna_odm_fp"
  path "response_to_pod_map_fp"
  val "pair_pod"

  output:
  path "result_${pair_pod}*"

  """
  run_association_analysis_nuclear.R $sceptre_object_fp \
  $response_odm_fp \
  $grna_odm_fp \
  $response_to_pod_map_fp \
  $pair_pod \
  ${params.n_nonzero_trt_thresh} \
  ${params.n_nonzero_cntrl_thresh} \
  ${params.trial} \
  ${params.use_parquet}
  """
}

/***************
* MAIN WORKFLOW
***************/
workflow {
  // 0. set analysis parameters
  if (step_rank >= 0) {
    set_analysis_parameters(
      Channel.fromPath(params.sceptre_object_fp, checkIfExists : true),
      Channel.fromPath(params.response_odm_fp, checkIfExists : true),
      Channel.fromPath(params.grna_odm_fp, checkIfExists : true),
      Channel.fromPath(params.formula_object, checkIfExists : true),
      Channel.fromPath(disc_pairs, checkIfExists : true),
      Channel.fromPath(params.positive_control_pairs, checkIfExists : true)
    )
  }

  if (step_rank >= 1) {
  // 2. obtain the gRNA info
  prepare_assign_grnas(
    set_analysis_parameters.out.sceptre_object_ch.first(),
    Channel.fromPath(params.response_odm_fp, checkIfExists : true),
    Channel.fromPath(params.grna_odm_fp, checkIfExists : true),
    params.grna_assignment_formula != "${baseDir}/resources/placeholder_file.rds"
  )

  // 3. process output from above process
  grna_to_pod_map_ch = prepare_assign_grnas.out.grna_to_pod_map_ch.first()
  grna_pods_ch = prepare_assign_grnas.out.grna_pods_ch.splitText().map{it.trim()}
  grna_assignment_args_ch = prepare_assign_grnas.out.grna_assignment_args_ch.first()
  grna_assignment_method_ch = prepare_assign_grnas.out.grna_assignment_method_ch.splitText().map{it.trim()}.first()
  grna_assignment_formula_in_use_ch = prepare_assign_grnas.out.grna_assignment_formula_in_use_ch.splitText().map{it.trim() == "true"}.first()
  grna_assignment_formula_input_ch = grna_assignment_formula_in_use_ch.map{in_use ->
    file(
      in_use ? params.grna_assignment_formula : "${baseDir}/resources/placeholder_file.rds",
      checkIfExists : true
    )
  }
  prepare_assign_grnas.out.grna_assignment_warning_ch
    .splitText()
    .map{it.trim()}
    .filter{it}
    .subscribe{log.warn(it)}

  // 4. assign gRNAs
  assign_grnas(
    set_analysis_parameters.out.sceptre_object_ch.first(),
    Channel.fromPath(params.response_odm_fp).first(),
    Channel.fromPath(params.grna_odm_fp).first(),
    grna_to_pod_map_ch,
    grna_pods_ch,
    grna_assignment_method_ch,
    grna_assignment_args_ch,
    grna_assignment_formula_input_ch
  )
  dummy_assign_grnas(
    grna_assignment_method_ch
  )
  
  // 5. process output from above process
  grna_assignments_ch = assign_grnas.out.grna_assignments_ch.mix(dummy_assign_grnas.out.grna_assignments_ch).collect()
  grna_assignment_formula_ch = assign_grnas.out.grna_assignment_formula_ch.mix(dummy_assign_grnas.out.grna_assignment_formula_ch).first()
  
  // 6. process the gRNA assignments
  combine_assign_grnas(
    set_analysis_parameters.out.sceptre_object_ch.first(),
    Channel.fromPath(params.response_odm_fp).first(),
    Channel.fromPath(params.grna_odm_fp).first(),
    grna_assignment_args_ch,
    grna_assignment_formula_ch,
    grna_assignments_ch
  )
  }

  if (!nuclear) {
    if (step_rank >= 2) {
    // 7. run quality control
    run_qc(
      combine_assign_grnas.out.sceptre_object_ch,
      Channel.fromPath(params.response_odm_fp).first(),
      Channel.fromPath(params.grna_odm_fp).first(),
    )
    }

    if (step_rank >= 3) {
    // 8. prepare association analyses
    prepare_association_analysis(
      run_qc.out.sceptre_object_ch,
      Channel.fromPath(params.response_odm_fp).first(),
      Channel.fromPath(params.grna_odm_fp).first()
    )

    // 9. run calibration check
    calibration_check_pods_ch = prepare_association_analysis.out.calibration_check_pods_ch.splitText().map{it.trim()}
    run_calibration_check_ch = prepare_association_analysis.out.run_calibration_check_ch.splitText().map{it.trim()}.first()
    run_analysis_subworkflow_calibration_check(
      prepare_association_analysis.out.sceptre_object_ch,
      Channel.fromPath(params.response_odm_fp).first(),
      Channel.fromPath(params.grna_odm_fp).first(),
      calibration_check_pods_ch,
      run_calibration_check_ch,
      Channel.from("run_calibration_check").first()
    )
    }

    if (step_rank >= 4) {
    // 10. run power check
    power_check_pods_ch = prepare_association_analysis.out.power_check_pods_ch.splitText().map{it.trim()}
    run_power_check_ch = prepare_association_analysis.out.run_power_check_ch.splitText().map{it.trim()}.first()
    run_analysis_subworkflow_power_check(
      run_analysis_subworkflow_calibration_check.out.first(),
      Channel.fromPath(params.response_odm_fp).first(),
      Channel.fromPath(params.grna_odm_fp).first(),
      power_check_pods_ch,
      run_power_check_ch,
      Channel.from("run_power_check").first()
    )
    }

    if (step_rank >= 5) {
    // 11. run discovery analysis
    discovery_analysis_pods_ch = prepare_association_analysis.out.discovery_analysis_pods_ch.splitText().map{it.trim()}
    run_discovery_analysis_ch = prepare_association_analysis.out.run_discovery_analysis_ch.splitText().map{it.trim()}.first()
    run_analysis_subworkflow_discovery_analysis(
      run_analysis_subworkflow_power_check.out.first(),
      Channel.fromPath(params.response_odm_fp).first(),
      Channel.fromPath(params.grna_odm_fp).first(),
      discovery_analysis_pods_ch,
      run_discovery_analysis_ch,
      Channel.from("run_discovery_analysis").first()
    )
    }
  } else { // going nuclear
    if (step_rank >= 2) {
    // 7. run cellwise qc
    run_qc_trans(
      combine_assign_grnas.out.sceptre_object_ch,
      Channel.fromPath(params.response_odm_fp).first(),
      Channel.fromPath(params.grna_odm_fp).first(),
    )
    }

    if (step_rank >= 5) {
    // 8. prepare association analysis
    prepare_association_analysis_trans(
      run_qc_trans.out.sceptre_object_ch,
      Channel.fromPath(params.response_odm_fp).first(),
      Channel.fromPath(params.grna_odm_fp).first(),
    )
    discovery_analysis_pods_ch = prepare_association_analysis_trans.out.discovery_analysis_pods_ch.splitText().map{it.trim()}
    response_to_pod_map_ch = prepare_association_analysis_trans.out.response_to_pod_map_ch

    // 9. run association analysis
    run_discovery_analysis_trans(
      run_qc_trans.out.sceptre_object_ch,
      Channel.fromPath(params.response_odm_fp).first(),
      Channel.fromPath(params.grna_odm_fp).first(),
      response_to_pod_map_ch,
      discovery_analysis_pods_ch
    )
    }
  }
}
