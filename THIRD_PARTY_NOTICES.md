# Third-party notices

The MIT license covers project-authored code. It does not relicense external software, database content, or third-party portions of saved configuration files.

## InterProScan

The two saved effective configurations under `historical/hummingbird/outputs/Hummingbird_results_20261004/provenance/` contain upstream InterProScan configuration and contributor metadata together with run-specific settings. InterProScan is developed by the EMBL-EBI InterPro team and is distributed under Apache 2.0. Its license is preserved in `LICENSES/InterProScan-Apache-2.0.txt`.

Source and licensing: https://github.com/ebi-pf-team/interproscan6 and https://github.com/ebi-pf-team/interproscan6/blob/main/LICENSE . Historical source copies were sanitized only as recorded in `provenance/source-manifest.json`.

## Other dependencies and annotation resources

NumPy, pandas, SciPy, Matplotlib, seaborn, Nextflow, SingularityCE, NCBI BLAST, Phobius, TMbed, and the InterPro member databases remain the work of their respective authors. They are external dependencies, not project-authored code. Database, model, and separately licensed application terms apply to their use and distribution.

The repository includes saved annotation output as an analysis input, not the underlying annotation databases or model installations. Cite the tools and database versions actually used in any derived research.
