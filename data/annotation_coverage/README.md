# Annotation-coverage inputs

These two preserved project files are sufficient for the portable coverage-analysis script. SHA-256 digests are in `../../provenance/located-inputs.json`.

## Manifest

`hypothetical_proteins_manifest_with_counts.csv` contains a header and 2,015 protein records.

| Field | Meaning in the source file |
| --- | --- |
| `gene` | Protein-group identifier used to join annotation records |
| `annotation` | Pre-existing annotation text; may contain multiple descriptions |
| `genome` | Reference genome/assembly identifier |
| `source` | Source feature identifier or reference |
| `seq_length` | Protein sequence length in amino acids |
| `note` | Reserved notes column; empty in all included records |
| `n_isolates` | Number of isolates associated with the protein group |

These are bacterial protein/genome records, not patient or participant records. Original field contents and row order are preserved.

## InterProScan table

`hypothetical_proteins_combined.faa.tsv` contains 2,871 rows and no header. Its 15 fields are: protein identifier, sequence MD5, sequence length, member database, signature accession, signature description, start, end, score, status, run date, integrated InterPro accession, integrated InterPro description, GO terms, and pathways.

One protein may have many signature rows; multiple rows can refer to the same integrated InterPro entry. Proteins with no TSV rows are represented by the manifest. An absent TSV row by itself does not prove that the original annotation job completed successfully.

This is a preserved historical output, not a refreshed database search. Database and sequence-source terms continue to apply; the repository's MIT grant concerns project-authored code.
