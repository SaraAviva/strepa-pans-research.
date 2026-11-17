# Archive coverage and remaining sources

This register distinguishes recovered code from work known only through reports or project notes. The initial archive is not the complete source for every result in the research project.

| Workstream | Preserved here | Remaining gap |
| --- | --- | --- |
| 2,015-protein coverage | Portable script, two inputs, tested dependencies, original PDF and reconstructed code | Original InterProScan batch submission and completion evidence |
| Initial annotation pilots | References in historical source/notes | Original pilot inputs and commands if retained in final methods |
| Outlier/sensitivity analysis | Workbook located in the project upload, not published here | Generating code and settings |
| BLAST | Saved report located, not substituted for source code | Commands, parsers, database versions, and query/response provenance |
| 19-candidate and three-priority reviews | Research scopes documented | Final selection/filtering source and the canonical candidate-table version |
| October InterProScan run | Source, configuration, recovery versions, final verifier, synthetic tests | Portable cluster configuration and full external run metadata |
| Supplementary predictions | Execution and verification source | Full environment and output-to-source mapping |
| Figures and tables | Coverage plotting and workbook generation | Region-track figure source and other final manuscript figure code |
| Upstream tiling and mapping | File index in project notes | Source directory is not readable by the connected Hummingbird account |

## Hummingbird audit

The starting folder was `StrepA_project_revamped`. Thirty-nine source files from the upload's `02_scripts`, both October run script directories, the recovery folder, saved source/configuration provenance, and the runtime wrapper were retrieved read-only and verified against the original local-source hashes. No new source versions were found in those locations. See `../provenance/hummingbird-source-audit.json`.

The hypothetical-protein directory contained 2,015 sequence filenames and no scripts. Third-party installations, containers, Nextflow task/cache trees, and participant datasets were excluded from source collection. This was a targeted source audit, not a claim that every generated file on Hummingbird was downloaded.

## Adjacent-folder check

The adjacent `blast_tasks` folder contained five FASTA files and no scripts. `strepA_tasks` contained older RAxML job scripts and phylogenetic outputs. Two readable job-script versions were retained locally for scope review, but are not included here because their relationship to this antibody-research project has not been confirmed; one additional saved backup was unreadable. The nested `panstripe_tasks` directory contained data/tree files, not source scripts.

## Missing upstream source requested from project maintainers

The project notes identify `SerImmune_project/scripts` as the source location. Access returned **permission denied** on November 17, 2026. The following names are from the notes, not verified files:

- `flare_vs_control_tiling_analysis_4.0_allproteins.R`
- `gene_category_enrichment.R`
- `investigate_flagged_hits.R`
- `map_old_to_new_pangenome.R`
- `map_significant_regions_all_2.0.R`
- `map_significant_regions_all.R`
- `map_significant_region_test.R`
- `run_tiling_all_proteins.sh` (case needs verification)
- `summarize_region_mapping_confidence.R`

To complete the archive, a maintainer can supply copies of these scripts and any BLAST, sensitivity-analysis, or figure-generation sources. Participant files are not needed for source collection. Each additional source should receive a provenance record and be checked before publication.
