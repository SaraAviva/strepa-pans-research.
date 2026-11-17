# Data availability and software

The two small annotation-coverage inputs in `../data/annotation_coverage/` contain bacterial protein/genome metadata and saved InterProScan annotations. Their exact bytes and schemas are documented. These inputs enable the portable analysis; they do not reproduce every project stage.

Participant-level datasets, patient classifications, private workbooks, unpublished notes, connection records, credentials, installed databases, binaries, and containers are excluded. The project file index identifies `Cases_controls.csv` as containing human-subject IDs and patient classifications; it is not part of this release.

Historical scripts can reference controlled inputs that are not distributed. Researchers should obtain any needed restricted data through the project maintainers and the applicable access process. No public access to those datasets is implied by publication of the code.

The manifest records privacy substitutions and both source and publication-copy checksums. Public upstream contributor metadata is retained. Private spreadsheet IDs, personal account details, and temporary session paths are replaced with explicit placeholders.

Project-authored code is MIT-licensed. Third-party software and database terms are unchanged; see `../THIRD_PARTY_NOTICES.md`. This repository does not redistribute InterProScan, BLAST, prediction packages, or installed databases.
