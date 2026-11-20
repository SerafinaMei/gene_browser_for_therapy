## Full RUN_MULTIALLELIC import · 2026-11-10

The source-folder audit found a second cause of undercounting: the public all-protein-coding track is smaller than the full run for many D&D genes. RUN_MULTIALLELIC has larger excision counts for 454 genes when compared with the previously imported headline tracks (for example, RARB 2,427 versus 43 and RYR2 2,201 versus 50). This is a comparison of named datasets with different selection histories, not a rule to take the maximum across runs.

`Grace_pipeline/high_coverage_result/RUN_MULTIALLELIC/summary_files/browser_tracks/per_gene_files/` contains 568 complete BEDs, totaling 135,549 variant records. Every one of the four per-strategy counts reconciles to `summary_files/dnd_summary.csv`, covering 593 genes: 568 with candidates and 25 assessed with zero candidates. The identical summary in `261011_RUN_final` was confirmed by SHA-256. The aggregate BEDs were empty or truncated and are deliberately excluded.

`scripts/import_pipeline.py` validates and imports these complete per-gene files, retains GRCh38 allele identities and source names, and intersects them with the verified public cell-line tracks. The source summary, per-file hashes and validation are retained. Counts and downloads appear as RUN_MULTIALLELIC evidence. These results are prominent for assessed genes; the dedicated SPTLC1/2 track remains prominent for those two genes. Other published tracks and historical selected-guide exports remain separately labelled. No cross-run counts are summed.

The GitHub code checkout contains pipeline code rather than another completed genome-wide result. The later 1,007-gene window-extension and November 4 selected-guide outputs have different transcript, window and expression criteria; their existence does not establish genome-wide completeness. Genes outside the 593-gene run still use available published tracks, with no claim that those enumerate every possible qualifying variant.

## Comprehensive evidence update · 2026-11-10

The report now loads a separate, reproducible evidence layer built by `scripts/refresh_evidence.py`. It replaces the headline ClinVar counts with the official November 9, 2026 summary, imports current IMPC phenotype and disease records, and adds the genome-wide target track. The earlier cohort results and selected-guide downloads remain available in expandable sections. See `data_sources/comprehensive/README.txt` for source and rebuild details.

Targeting covers 19,181 current HGNC entries, resolving 67 official previous symbols and adding 20 genes missing from the earlier catalog. The dedicated SPTLC1/2 track is used for those genes, with the distinct all-protein-coding dataset available separately. Counts use unique GRCh38 interval/allele keys. Cell-line intersections require the same allele pair, preserve HIGH/LOW quality, and never substitute locus-only matches. Strategy counts overlap; source PAM flags are separate from selected gRNAs. These are complete imports of the named snapshots, not a claim that the underlying research pipeline exhaustively assessed every possible variant.

The mouse page shows phenotype terms, supporting datasets, allele, zygosity, life stage, sex, manual assertions, physiological systems, and named human disease associations with identifiers. Manual assertions are not treated as statistical p-values. The current PhenoDigm API does not supply phenotype matches for SPTLC1's four associated diseases, so the table says that explicitly rather than inferring matches from the pasted webpage. Historical cohort aggregates remain available separately.

All 919,143 broad-track input rows are accounted for. The current IMPC API snapshot is paginated with response-count checks. SHA-256 hashes, exact source URLs, gene mappings, and cell-line mismatch counts are retained in the source/audit files. `python3 tests/test_comprehensive.py` reconciles every gene's displayed variant/quality/strategy counts with its downloadable CSV, verifies mouse snapshots, and checks current ClinVar/browse consistency.

The footer is “Test version.” The requested generic evidence paragraph and old summary-count footnote have been removed. ClinGen and the historical selected-guide exports retain their explicitly dated original sources; this update does not claim those were rerun.

# Scientific definitions and data coverage

## Scientific definitions

**Collins pHaplo is probability of haploinsufficiency (0–1).** It is never inverted or labelled as a probability of haplosufficiency. pTriplo is separate. ClinGen's `%HI` is the DECIPHER HI percentile rank, where lower values indicate greater predicted HI. It is not a calibrated probability. The browser does not fabricate a combined haplosufficiency probability.

**s_het** uses `post_mean`, `post_lower_95`, and `post_upper_95` from Supplementary Table 1 of the supplied Spence workbook. The interval is a Bayesian **credible** interval. These values are read from the published table, not recalculated. The similarly named `sHet` field in the Collins table is retained only as a source field and is not substituted for GeneBayes.

**ClinGen** disease/MOI/classification records are retained intact from the November 13, 2025 export. HI = 0 means no evidence, not proven haplosufficiency. HI = 30 indicates an autosomal recessive condition; HI = 40 indicates dosage sensitivity unlikely. Gene and regional dosage records remain distinct.

**D&D membership** is displayed as Yes/No under Dosage and constraint. It reports membership in the supplied 593-gene set, not a new biological classification. The definition links to [Ramey et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC13042096/), which describes dominant and dispensable disease alleles as potential targets for allele-specific inactivation while preserving a functional copy. No means absent from this supplied set.

**Variants:** current headline totals come from the official ClinVar gene summary; its date and hash appear in Data sources. `Total_alleles` counts distinct alleles in or overlapping a gene, separately from submissions. The weekly pipeline independently derives the dominant P/LP mapped subset described below. Historical D&D counts remain under Additional counts. Counts from different snapshots or inclusion rules are not used to calculate a dominant fraction.

**Record downloads:** the separate local `haplo_prelim_pipeline/data/ClinVar/variant_summary.txt` supplies GRCh38 variant records matched by NCBI Gene ID. Its snapshot date is unspecified, and downloaded record counts may differ from the dated ClinVar summary above. All source columns are retained. These historical GRCh38 records are separate from the current Dominant assessment download. Available records cover 19,930 genes. Exported text that could be interpreted as a spreadsheet formula is prefixed with an apostrophe.

**Targeting:** exonic indel, NMD-inducing indel, NMD-escaping indel, splice base editing, exon/gene excision, and CRISPRoff counts retain their source field names and definitions. Viable PAM flags are gene-level indicators, not counts of post-PAM variants. Excision pairs are distinct from variants. Strategies can overlap, so counts are not summed. The cell-line dataset is a recovered subset with incomplete source coverage. The 1-23 bp and 1-30 bp settings refer to SNP-to-splice-boundary distance, not editor activity windows; both include exonic SNPs. WTC and WTD have allele-matched genotype calls, while p47 confirms heterozygous loci without allele identities. WTD FILTER is unspecified in the source. Zero means no observed heterozygous candidate in this dataset, not confirmed homozygosity.

The earlier cell-line SNP and guide downloads cover 422 genes and preserve all original association rows, including selection flags, editor, window, sequence background and genotype evidence. A variant or guide can occur in multiple rows. The visible table counts unique variants/guides instead. The selected-guide counts reconcile with every row of the existing gene-level summary. Guide sequences use GRCh38 plus the focal allele only, with neighboring variants unphased; these are candidates, not experimentally validated guides.

**Mouse evidence:** human-to-mouse mapping uses MGI homology identifiers. Successful heterozygous statistical rows from the supplied DND593, Dominant855, HI417 and AR754 exports are deduplicated by `doc_id` across cohorts. `het_rows` includes significant and non-significant rows. `het_significant_rows` preserves the source `significant` flag without applying a new threshold. `het_parameters_tested` counts distinct nonempty `parameter_name` values; `het_procedures_tested` counts procedure names. Rows are statistical results, not animals. `has_hom_record` records a homozygote assertion in the supplied exports. False means not found there, not proven absence of homozygous testing. `het_class` is an explicitly derived all-phenotype summary, not the older mortality/neurobehavior-only enrichment class. The presence of a homozygote record does not imply that heterozygotes are phenotypically normal or documented breeders.

**CMT/PMP22:** the CMT flag matches “Charcot-Marie-Tooth” in ClinGen disease labels and retains all evidence classifications, including limited/disputed assertions. PMP22's gene-level HI evidence remains separate from its regional duplication/CMT1A context. Other CMT-associated genes are not assigned the PMP22-specific disease group.

## Targeting coverage audit (2026-11-01)

The original strategy table used only the 593-gene D&D workbook. The local `Grace_pipeline/high_coverage_result/RUN_MULTIALLELIC/summary_files/dnd_summary.csv` and `Grace_pipeline/261011_RUN_final/summary_files/dnd_summary.csv` also contain 593 unique genes (737 rows), and their SHA-256 hashes are identical. The intermediate `RUN_SMEI_STAGEA_FILTERED_0608` summaries are restricted to the same cohort. Older candidate dictionaries and guide archives are also limited cohorts; no completed genome-wide result was found locally. The quoted Wynton directory is not mounted, so these copies cannot be confirmed as the latest remote revision.

`Grace_pipeline/pipeline_recheck_20261101/` contains a broader completed indel/splice run with 1,007 input genes: the 593 D&D genes, 417 ClinGen HI genes (five overlap), and DSG2/TTN controls. It assessed 1,006 genes, including 413 outside the D&D set. SRY lacks a chromosome-Y allele-frequency input. NMD is unassessed for MIR96, TERC and XIST because a retained transcript lacks coding coordinates. Excision and CRISPRoff were not run in this expanded dataset.

The browser now imports this run without changing D&D membership. Genes with both sources have a dataset selector; the original reference counts remain available. Other assessed genes use the expanded results, with SNP and gRNA downloads when records exist. The expanded cell-line results are also available where the earlier candidate dataset has no records. These datasets use different transcript models and selection criteria and are kept separate.

Expanded strategy counts use distinct common SNPs (AF 0.1-0.9) in shared exonic or 1-23 bp splice regions, before cell-line selection. They do not count ClinVar pathogenic variants. Targetable means at least one candidate in that strategy; viable PAM means at least one candidate has an NGG guide satisfying the source 10 bp SNP-to-PAM rule. NMD counts reconcile to the raw NMD summary. Induction and escape require consistent transcript predictions and need not sum to the pre-NMD count. The cell-line selected-guide counts additionally apply the source expression, genotype, category, coding/NMD and splice-coverage flags. Downloads retain both 1-23 and 1-30 bp settings and all source flags.

A zero is shown only for an assessed strategy with no candidates. Genes outside the available runs show an explanatory message; unassessed NMD and strategies not run are labelled explicitly. `dist/data/targeting-coverage.json` records coverage, exclusions, examined file hashes and source validation. `scripts/build_strategy_screen.py` performs the import as part of the record build.

A genome-wide extension would require another pipeline run. The [upstream pipeline](https://github.com/gramey02/dnd_project) uses configurable gene/exon inputs and SGE wrappers, with Python entry points available for local execution. The configured all-protein-coding exon input and complete phased/biallelic VCF inputs were not found in this workspace, and `qsub` is unavailable locally. GRCh38, GENCODE, GTEx and chromosome 1-22/X allele-frequency files are present, so a broader Stage A run is possible after rebuilding the gene/transcript input; a complete run of all strategies needs its remaining inputs and execution configuration. No genome-wide computation was launched during this browser update.



## ClinGen validity refresh — 2026-11-10

The official complete ClinGen validity export now supplies 3,730 assertions for 3,078 genes, joined by HGNC ID with no unmatched records. It replaces the 2025-11-13 validity snapshot; 847 gene records changed. SPTLC1 includes MONDO:0008086, autosomal dominant, Definitive, evaluated 2026-03-31. Gene-disease validity and dosage sensitivity remain separate sources.

The source CSV is retained at `data_sources/clingen-gene-validity.csv`; its date and SHA-256 are recorded in the manifest. Run `python3 scripts/refresh_clingen.py --fetch` to download and validate a new export, replace gene disease arrays, and update source metadata and the search index. Run without `--fetch` to rebuild from the cached export. The full data builder uses the same cached validity source.


## Population haplotype coverage

The Figure 5A-style plots use 1,093 completed greedy-prioritization results for 562 autosomal genes from the saved GRCh38 analysis. The importer checks completion markers, cumulative gains, matching summary totals, and identical ordered sample headers for every input VCF. Chromosomes are read from the variant files and cross-checked against nonblank master annotations; missing summary annotations are never classified as sex chromosomes. The verified denominator is 2,504 individuals / 5,008 haplotypes, rather than the 5,096 haplotypes in the pasted paper excerpt. These plots describe this saved analysis and do not reproduce the paper's numerical values exactly.

Blue curves count distinct covered haplotypes cumulatively, excluding homozygous genotypes. Gray dashed lines show twice the number of people heterozygous at an eligible site before PAM filtering, divided by 5,008. These are population coverage estimates, not observed patient outcomes or estimates weighted for a particular ancestry.

For single-guide strategies, the red dashed model is min(number of therapies / dominant disease allele count, 1), assuming equally frequent disease alleles. The default comparison uses the current mapped dominant P/LP subset. The Original analysis option uses the matched 2026-08-25 master independently of current ClinVar annotations. For MYH7, that original comparison uses 149 alleles; exon-disruption coverage is 1,206 / 5,008 (24.1%) with the first therapy and 3,344 / 5,008 (66.8%) with all eight. The pre-PAM ceiling is 3,736 / 5,008 (74.6%).

Excision can reuse one guide in a new pair, so its x-axis counts unique gRNAs and has no mutation-specific comparison. Sex-chromosome results for EIF2S3 are withheld pending a ploidy-specific denominator. SPTLC1 has no saved greedy result; candidate variant and cell-line counts are not converted into a population curve.

`data_sources/population-coverage.json.gz` preserves aggregate points, source CSV hashes, the master-file hash, cohort size and an ordered sample-header hash. It contains no sample identifiers or individual genotypes. `scripts/build_population_coverage.py` rebuilds standalone SVGs, gzip CSV downloads and per-gene shards. The paper is linked from each available plot; technical assumptions and data tables are collapsed by default.


## MYH7 coverage interpretation audit (2026-11-18)

The denominator is 2,504 diploid samples multiplied by two phased chromosome copies, yielding 5,008 observed haplotypes. It is not a count of 5,008 unique sequence patterns, nor a count of affected patients. Heterozygosity describes the pair of alleles in a person, not an individual haplotype.

Independent parsing of the saved MYH7 exon-disruption VCFs recovered 5 eligible common biallelic SNP positions before PAM filtering (1,868 distinct heterozygous people) and 4 positions afterward (1,672 distinct heterozygous people). The run's recorded AF_LIMIT is 0.1. Selected alleles at GRCh38 positions 23433544 (G then A), 23429785 (G then A), 23431585 (G then A), and 23415680 (C then T) add 1,206, 1,206, 243, 243, 144, 144, 79 and 79 haplotypes, respectively. These gains were recomputed from phased heterozygous calls after subtracting already covered chromosome copies.

The first two guides provide alternative targeting options for opposite chromosome copies in the same 1,206 heterozygous people. They do not imply 2,412 different people or that both guides should be given to a patient. Final coverage is 3,344 / 5,008 = 66.8%; this corresponds to both chromosome copies represented among the library's alternatives in 1,672 people. A patient-specific choice needs the editing site phased to the disease-causing mutation.

The curve does not test for a ClinVar pathogenic variant in the same 1000 Genomes individual or on the same chromosome copy. The number of haplotypes carrying both an eligible common heterozygous editing site and a confirmed pathogenic MYH7 allele is therefore unknown from this analysis, not zero. The historical 149 dominant disease alleles used by the red equal-frequency comparison are a separate catalogue count, not 149 carriers or haplotypes in this cohort.

## OMIM relationship table availability

The supplied `OMIMgenemap2.txt` was generated on 2026-11-18. Its full parsed relationship records, source gene identifiers and original SHA-256 are cached in `data_sources/omim-gene-phenotypes.json`. This replaces the initial FUS-only table.

The import maps 17,660 browser genes to OMIM entries; 5,127 have phenotype relationships, yielding 7,801 displayed gene-phenotype rows. Gene MIM numbers and phenotype MIM numbers remain separate. Gene MIM IDs provide the primary join. A unique Entrez ID plus an agreeing approved symbol or recorded alias recovers 295 additional OMIM entry assignments, including nine additional genes with phenotype relationships. No association is inferred from genomic proximity or a similar disease name.

The join audit records fallback assignments, any identifier conflicts and unmatched phenotype-bearing MIM entries. Many unmatched entries describe loci, intervals or syndromes without a specific gene assignment. Genes matched to an entry with no phenotype record are distinguished from genes without an imported entry. OMIM qualifiers (`?`, `{}`, `[]`), inheritance strings, cytogenetic locations and mapping keys are preserved exactly. Source dates remain available in Data sources.


## Weekly dominant-variant refresh

`scripts/refresh_weekly.py` refreshes ClinGen validity, ClinVar gene totals and the complete weekly ClinVar VCV XML before publishing. The release is streamed from gzip, checked against NCBI's MD5, and recorded with SHA-256 and release date. The job fails on incomplete releases, unexpected schemas, stale dates, unmapped ClinGen assertions or conflicting gene-summary records. The live site is replaced only after the data checks, interface tests and hosting build pass.

The dominant count is a conservative, reproducible subset:

1. Require a current classified SimpleAllele assigned to exactly one gene with a matching NCBI Gene ID. Included-only records, compound haplotypes/genotypes and multigene alleles do not qualify.
2. Require germline Pathogenic, Likely pathogenic or Pathogenic/Likely pathogenic classification with assertion criteria at both VCV and RCV levels. Conflicting classifications and records without criteria do not qualify.
3. Resolve inheritance for each RCV condition. Prefer explicit inheritance from current, contributing P/LP ClinVar submissions with assertion criteria that match that same condition. Otherwise use exact condition IDs mapped to ClinGen Moderate/Strong/Definitive assertions or OMIM mapping-key-3 relationships. ClinVar's gene-condition mapping bridges MedGen, MONDO and OMIM disease IDs; gene-wide inheritance labels and disease-name similarity are not used.
4. Require dominant-only inheritance (autosomal dominant, X-linked dominant or semidominant). Mixed dominant/recessive mappings, unspecified X-linked modes and missing mappings are unresolved. Every condition in a multi-condition RCV must resolve; a known dominant condition cannot fill an unknown one. Explicit variant-condition inheritance may resolve a disease that has multiple inheritance modes.
5. Count each AlleleID once per gene if at least one qualifying condition resolves as dominant. Unresolved counts exclude alleles already counted as dominant. Known recessive-only records are excluded. Genes without an NCBI Gene ID are unavailable, not zero.

`data_sources/dominant-current.json.gz` retains qualifying and unresolved variant-condition rows with RCV/Variation/Allele identifiers, classification, review status, condition and inheritance basis. Public per-gene downloads reconstruct these same rows; the tests reconcile every count and browse index entry against them. The count is not an exhaustive catalog of disease-causing alleles. Changes from the paper can reflect stricter inclusion rules as well as newer submissions or reclassifications.

The supplied OMIM inheritance snapshot remains dated 2026-11-18 until a new authorized export is imported. ClinGen validity updates do not alter dosage sensitivity records. Editing and cell-line candidates, IMPC data and the original 1000 Genomes results retain their own sources and dates.

`scripts/refresh_population_comparison.py` uses the new mapped count only for the equal-frequency mutation-specific model, min(guides / count, 1). A zero count produces no red model. Every refreshed CSV retains the identical cohort coverage, guide sequence choices and pre-PAM ceiling from the original result. Current model downloads include the clinical release date, count and unresolved count. The original SVG/CSV and its baseline remain available with the Original analysis selector. This refresh does not establish disease status or pathogenic-variant phase in the cohort.
