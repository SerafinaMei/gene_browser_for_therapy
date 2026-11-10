Comprehensive evidence refresh, 2026-11-10

Public input snapshots:
- IMPC Solr gene: 25,017 records; genotype-phenotype: 67,350 assertions;
  PhenoDigm disease_gene_summary: 17,530; curated disease_model_summary: 108,601.
  The API response counts were checked during pagination. These API record counts
  are not interchangeable with the IMPC release headline phenotype-call total.
- NCBI ClinVar gene summary dated 2026-11-09.
- All-protein-coding UCSC track, downloaded and SHA-256 matched to the local
  lossless BED export (919,143 rows, 19,182 original symbols).
- Dedicated published SPTLC1/2 track: 434 / 541 rows. This is a separate run and is
  shown preferentially for these two genes; the broad track is retained separately.
- Current public WTB, WTC, WTD and KOLF2 bigBeds, verified against GitHub blob IDs.

source-verification: upstream-verification.json contains exact URLs and hashes.
mouse-sources.json contains query URLs, counts, hashes and retrieval timestamps.

The browser maps exact approved HGNC symbols and official previous symbols,
not general aliases. 67 previous names resolve to current entries. Two source
symbols map to DDI2; distinct alleles are counted once and source annotation rows
remain separate in CSV. Twenty source genes omitted by the old catalog are added.
No original research inputs are changed. No targetability inference is made for a
missing gene, and a source PAM flag is not presented as a selected guide count.

To rebuild from the local source folders:
  python3 scripts/refresh_evidence.py --build
Optional --track-bed, --cell-dir and --bigbed-tool set alternative resource paths.
--reuse-tracks reuses only a previously built import whose input hashes still match.
To download fresh IMPC data: python3 scripts/refresh_evidence.py --fetch-mouse
Use --resume-mouse only to resume an interrupted same-release download.
The current release label is 24.0; verify the IMPC release page before future refreshes.

Large raw bigBed/BED files and the intermediate genome-evidence cache are local build
inputs, not required for serving or deployment. Public source URLs above allow them
to be retrieved again. The generated dist/data/evidence-*.json and downloadable
CSVs contain the served records. Do not rerun the older cohort builders alone after
this import; run refresh_evidence.py --build afterwards to refresh the browse index.

## Full RUN_MULTIALLELIC import · 2026-11-10

The source-folder audit found a second cause of undercounting: the public all-protein-coding track is smaller than the full run for many D&D genes. RUN_MULTIALLELIC has larger excision counts for 454 genes when compared with the previously imported headline tracks (for example, RARB 2,427 versus 43 and RYR2 2,201 versus 50). This is a comparison of named datasets with different selection histories, not a rule to take the maximum across runs.

`Grace_pipeline/high_coverage_result/RUN_MULTIALLELIC/summary_files/browser_tracks/per_gene_files/` contains 568 complete BEDs, totaling 135,549 variant records. Every one of the four per-strategy counts reconciles to `summary_files/dnd_summary.csv`, covering 593 genes: 568 with candidates and 25 assessed with zero candidates. The identical summary in `261011_RUN_final` was confirmed by SHA-256. The aggregate BEDs were empty or truncated and are deliberately excluded.

`scripts/import_pipeline.py` validates and imports these complete per-gene files, retains GRCh38 allele identities and source names, and intersects them with the verified public cell-line tracks. The source summary, per-file hashes and validation are retained. Counts and downloads appear as RUN_MULTIALLELIC evidence. These results are prominent for assessed genes; the dedicated SPTLC1/2 track remains prominent for those two genes. Other published tracks and historical selected-guide exports remain separately labelled. No cross-run counts are summed.

The GitHub code checkout contains pipeline code rather than another completed genome-wide result. The later 1,007-gene window-extension and November 4 selected-guide outputs have different transcript, window and expression criteria; their existence does not establish genome-wide completeness. Genes outside the 593-gene run still use available published tracks, with no claim that those enumerate every possible qualifying variant.
