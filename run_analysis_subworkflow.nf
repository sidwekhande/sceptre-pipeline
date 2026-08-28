// PROCESS: run association analysis
process run_association_analysis {
  time {params.run_association_analysis_time_per_pair * params.pair_pod_size}
  memory params.run_association_analysis_memory
  
  when:
  run_analysis == "true"

  input:
  path "sceptre_object_fp"
  path "response_odm_fp"
  path "grna_odm_fp"
  val "pair_pod"
  val "run_analysis"
  val "analysis_type"

  output:
  path "result.rds", emit: results_ch
  path "precomputations.rds", emit: precomputations_ch

  """
  run_association_analysis.R $sceptre_object_fp \
  $response_odm_fp \
  $grna_odm_fp \
  $pair_pod \
  $analysis_type
  """
}

// PROCESS: process association analysis results
process combine_association_analysis {
  publishDir "${params.output_directory}", mode: 'copy', overwrite: true, pattern: "*.png"
  publishDir "${params.output_directory}", mode: 'copy', overwrite: true, pattern: "*.txt"
  publishDir "${params.output_directory}", mode: 'copy', overwrite: true, pattern: "results_*"
  // This process runs once per subworkflow (calibration_check, power_check, discovery_analysis),
  // chained in that order, and each publishes to the same output_directory with overwrite: true --
  // same pattern already used for analysis_summary.txt and the plots above. So after a full run,
  // whichever object is here is the discovery-analysis one: assigned gRNAs, QC applied, and
  // discovery_result attached -- the only fully-analyzed sceptre_object the pipeline produces.
  publishDir "${params.output_directory}", mode: 'copy', overwrite: true, pattern: "sceptre_object.rds"
  
  time params.combine_association_analysis_time
  memory params.combine_association_analysis_memory

  when:
  run_analysis == "true"

  input:
  path "sceptre_object_fp"
  path "response_odm_fp"
  path "grna_odm_fp"
  path "results"
  path "precomputations"
  val "run_analysis"
  val "analysis_type"

  output:
  path "sceptre_object.rds", emit: sceptre_object_ch
  path "results_*"
  path "*.png"
  path "analysis_summary.txt"

  """
  process_association_analysis_results.R $sceptre_object_fp \
  $response_odm_fp \
  $grna_odm_fp \
  $analysis_type \
  results* \
  precomputations*
  """
}

// PROCESS: dummy channel
process dummy_process {
  time "1m"
  memory "1GB"
  
  when:
  run_analysis == "false"

  input:
  path "sceptre_object_fp"
  val "run_analysis"

  output:
  path "sceptre_object_fp", emit: sceptre_object_ch

  """
  """
}


workflow run_analysis_subworkflow {
  take:
  sceptre_object_ch
  response_odm_fp_ch
  grna_odm_fp_ch
  pods_ch
  run_analysis_ch
  analysis_type

  main:
  run_association_analysis(
    sceptre_object_ch,
    response_odm_fp_ch,
    grna_odm_fp_ch,
    pods_ch,
    run_analysis_ch,
    analysis_type
  )
  combine_association_analysis(
    sceptre_object_ch,
    response_odm_fp_ch,
    grna_odm_fp_ch,
    run_association_analysis.out.results_ch.collect(),
    run_association_analysis.out.precomputations_ch.collect(),
    run_analysis_ch,
    analysis_type
  )
  dummy_process(
    sceptre_object_ch,
    run_analysis_ch
  )
  output_sceptre_ch = combine_association_analysis.out.sceptre_object_ch.mix(dummy_process.out.sceptre_object_ch)

  emit:
  output_sceptre_ch
}
