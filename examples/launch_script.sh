#!/bin/bash
# Limit NF driver to 4 GB memory
nextflow pull timothy-barry/sceptre-pipeline
export NXF_OPTS="-Xms500M -Xmx4G"

##################
# OUTPUT DIRECTORY
##################
output_directory=$HOME"/sceptre_outputs"

###############################################################
# OPTION A: pre-built sceptre_object.rds + gene.odm + grna.odm
###############################################################
data_directory=$HOME"/sceptre_data/"
sceptre_object_fp=$data_directory"sceptre_object.rds"
response_odm_fp=$data_directory"gene.odm"
grna_odm_fp=$data_directory"grna.odm"

nextflow run timothy-barry/sceptre-pipeline -r main \
 --sceptre_object_fp $sceptre_object_fp \
 --response_odm_fp $response_odm_fp \
 --grna_odm_fp $grna_odm_fp \
 --output_directory $output_directory \
 --grna_assignment_method mixture \
 --pair_pod_size 1000 \
 --grna_pod_size 25 \
 --trial

##########################################################################
# OPTION B: build the sceptre object in-pipeline from raw Cell Ranger dirs
##########################################################################
cellranger_dir_1=$data_directory"cellranger_out/sample_1"
cellranger_dir_2=$data_directory"cellranger_out/sample_2"
per_guide_metadata_tsv=$data_directory"per_guide_metadata.tsv"
discovery_pairs_tsv=$data_directory"discovery_pairs.tsv"
positive_control_pairs_tsv=$data_directory"positive_control_pairs.tsv"

nextflow run timothy-barry/sceptre-pipeline -r main \
 --rna_directories "$cellranger_dir_1,$cellranger_dir_2" \
 --per_guide_metadata_tsv $per_guide_metadata_tsv \
 --moi high \
 --discovery_pairs $discovery_pairs_tsv \
 --positive_control_pairs $positive_control_pairs_tsv \
 --output_directory $output_directory \
 --grna_assignment_method mixture \
 --pair_pod_size 1000 \
 --grna_pod_size 25 \
 --trial
