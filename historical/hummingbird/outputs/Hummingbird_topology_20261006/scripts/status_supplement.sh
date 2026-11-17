#!/bin/bash
set -eu
squeue -j 796204,796208 -o "%.18i %.12T %.12M %.40R" || true
sacct -j 796170,796163,796204,796208 --format=JobID,State,ExitCode,Elapsed,MaxRSS -P
for f in /hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261006T234134Z_topology_c655cc0c/provenance/verification_submission.json /hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261006T234134Z_topology_c655cc0c/results/tmbed/verification.json /hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261006T234134Z_topology_c655cc0c/summary/run_summary.json; do
 if [ -f "$f" ]; then cat "$f"; fi
done
