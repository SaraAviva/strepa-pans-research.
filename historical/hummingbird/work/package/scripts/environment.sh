#!/bin/bash
set +u
module load miniconda3 nextflow/25.10.4 singularity-ce/singularity-ce.4.1.4
export PATH="/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f/runtime/bin:$PATH"
export NXF_HOME="/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f/runtime/nextflow_home"
export NXF_TEMP="/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f/runtime/tmp"
export NXF_OFFLINE=true
export NXF_VER=25.10.4
export NXF_OPTS='-Xms256m -Xmx4g -XX:ActiveProcessorCount=8'
export SINGULARITY_CACHEDIR="/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f/runtime/singularity_cache"
export SINGULARITY_TMPDIR="/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f/runtime/tmp"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
set -u
