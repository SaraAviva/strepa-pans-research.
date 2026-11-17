# Reproduction and validation

## Portable annotation-coverage analysis

The README command was run successfully on the included manifest and TSV in a fresh Python 3.12 environment on macOS arm64 on November 17, 2026. All 2,015 manifest proteins were processed. Dependency versions are pinned in `../environments/requirements-annotation.txt`; input fingerprints are recorded in `../provenance/located-inputs.json`.

The historical random seed (20260810) and 3,000 bootstrap resamples were retained. This reruns annotation-count summaries and figures from a saved TSV; it does not repeat InterProScan or establish that every historical annotation task finished. Output hashes and validation details are recorded in `../provenance/local-validation.json`.

The script distinguishes total signature rows from distinct integrated InterPro accessions and retains manifest proteins with no reported matches. These counts have different meanings and should not be labeled interchangeably.

## Original-versus-portable comparison

The PDF-recovered original was also run in the same environment, substituting only its fixed workspace path in a temporary copy. All shared columns for 2,015 proteins and all 68 association metrics matched the portable version exactly after mapping renamed summary keys. Three CSV tables were byte-identical; the isolate-bin table matched after mapping its filename and normalizing dash characters in labels.

The original additionally emits three derived protein columns, a descriptive-summary section, a raw-with-headers table, and a workbook payload. These outputs are omitted by the portable version. Both generated figures, but no pixel-equivalence claim is made. Details are in `../provenance/analysis-version-comparison.json`.

## Historical Hummingbird scripts

`historical/hummingbird/` preserves the original recovered `work/` and `outputs/` relationships. It is an execution archive, not a portable cluster installer. Scripts contain site-specific directories, modules, SLURM settings, and references to run metadata. Personal login details and private spreadsheet IDs have visible placeholders.

The recorded successful October run used InterProScan 6.0.1, InterPro 109.0, Nextflow 25.10.4, and SingularityCE 4.1.4. These are historical versions, not recommendations for current installations. Databases and containers are not distributed here.

The final verifier is `historical/hummingbird/outputs/Hummingbird_results_20261004/scripts/postprocess_final_20261004.py`. Earlier postprocessors remain for traceability and should not replace it. The final supplementary component runner is `historical/hummingbird/outputs/Hummingbird_topology_20261006/scripts/run_components_v4.py`; its preceding versions remain in the development archive.

## Checks requiring no participant data

`tools/check_archive.py` verifies publication-copy hashes, parses Python, and checks Bash syntax without importing or executing the research scripts.

`tools/run_existing_verifier_tests.py` runs the ten existing synthetic final-verifier cases in a temporary directory. It omits the original test script's separate real-data recheck, which requires the full result package. The cases cover successful completion and missing or inconsistent evidence, including failed tasks, changed sequences, and missing structured records.

The original local final result package was separately checked against its 107 delivery checksums; all matched. That result package is not bundled here.

## Limits

The complete research pipeline has not been reproduced end to end. Missing sources, candidate-selection provenance, additional figure code, and unconfirmed database versions are listed in `completeness.md`. The PDF reconstruction has a separate recovery note and must not be mistaken for an original source file.
