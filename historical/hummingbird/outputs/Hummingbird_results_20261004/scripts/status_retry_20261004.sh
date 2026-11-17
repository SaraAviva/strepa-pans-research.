#!/bin/bash
set -euo pipefail
run='/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f'
printf 'Run directory: %s\n' "$run"
ids=()
for role in analysis postprocess; do
 if [ -f "$run/provenance/retry_20261004_${role}_job_id.txt" ]; then
  id=$(cat "$run/provenance/retry_20261004_${role}_job_id.txt")
  printf '%s job: %s\n' "$role" "$id"
  ids+=("$id")
 fi
done
if [ "${#ids[@]}" -gt 0 ]; then
 jobs=$(IFS=,; echo "${ids[*]}")
 squeue -j "$jobs" -o '%.18i %.12P %.28j %.10T %.10M %R' || true
 sacct -j "$jobs" --format=JobID,JobName%28,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList
fi
if [ -f "$run/summary/RUN_SUMMARY.txt" ]; then cat "$run/summary/RUN_SUMMARY.txt"; fi
if [ -f "$run/logs/nextflow.log" ]; then tail -18 "$run/logs/nextflow.log"; fi
