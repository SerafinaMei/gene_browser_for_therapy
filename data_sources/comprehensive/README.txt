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
