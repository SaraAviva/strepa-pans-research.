#!/bin/bash
set -euo pipefail
run='/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f'
sacct -j 792566,792576 --format=JobID,JobName%28,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList
cat "$run/summary/RUN_SUMMARY.txt"
