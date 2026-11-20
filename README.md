# D&D Gene Browser

An interactive gene browser for clinical annotations, dosage sensitivity, evolutionary constraint, variant records, population coverage plots, targeting results and mouse phenotype evidence.

## Edit the website

- **Text, colors, example genes and basic sizing:** [dist/site-config.json](dist/site-config.json)
- **Appearance:** the Shiny-style layout is used throughout. The previous Classic version is preserved in Git.
- **Layout and detailed styling:** [dist/custom.css](dist/custom.css) for Classic, [dist/shiny.css](dist/shiny.css) for Shiny
- **Instructions:** [CUSTOMIZE.md](CUSTOMIZE.md)

The website uses plain HTML, CSS and JavaScript. It has no account requirement or runtime dependency on a hosted editing service.

## Preview locally

```sh
python3 scripts/serve.py
```

Open `http://127.0.0.1:8765/`. On a Mac, you can also double-click `Start Gene Explorer.command`.

## Publish with GitHub Pages

The intended repository is `SerafinaMei/gene_browser_for_therapy`.

1. Push these files to its `main` branch, including the `.github` folder.
2. Open **Settings > Pages** and select **GitHub Actions** as the publishing source.
3. Open **Actions > Deploy gene browser > Run workflow**, or push another change to `main`.
4. When the workflow succeeds, open the address shown under **Settings > Pages**.

For this repository, the expected address is `https://serafinamei.github.io/gene_browser_for_therapy/`. This is an expected address, not confirmation that deployment has completed.

Later commits to `main` automatically rebuild and publish at the same address. Revert a commit to restore an earlier version. You can edit the settings file directly in GitHub's web editor.

GitHub Free requires a public repository for Pages. Private-repository Pages availability depends on the account plan; the website itself is publicly accessible under ordinary Pages hosting. See [GitHub's setup documentation](https://docs.github.com/en/pages/getting-started-with-github-pages/creating-a-github-pages-site).

## Build and check

```sh
node --test tests/*.test.mjs
python3 scripts/prepare_hosting.py
```

The build writes the deployable site to `out/`. GitHub Actions publishes only that folder, not the project documentation or research scripts. All assets use relative paths so the site works under the repository's URL prefix. No theme, provider badge or external script is injected by the application.

Gene annotations are compressed during the build. Download bundles reconstruct the same CSV records in the browser. Current browsers are required. See [SCIENTIFIC_NOTES.md](SCIENTIFIC_NOTES.md) for source definitions and coverage limitations.

The site is labelled **Test version** and includes `noindex` metadata. Anyone with its address can view it and download the data. Search-indexing directives are not access controls.

## Population coverage plots

Gene reports show saved greedy-prioritization curves for 562 autosomal genes, with an editing-strategy selector and SVG/CSV downloads. See MYH7 for an example. Missing curves are identified explicitly, including SPTLC1.

Rebuild these plots with `python3 scripts/build_population_coverage.py` (requires Matplotlib). The committed aggregate cache is sufficient; individual genotypes are not distributed. The script docstring describes importing a new analysis. Run `python3 tests/test_population.py` to reconcile every curve and download.

## OMIM phenotype relationships

Clinical annotations now have an OMIM Gene-Phenotype Relationships table with phenotype MIM IDs, inheritance and mapping keys. The 2026-11-18 genemap2 export supplies 7,801 phenotype relationships across 5,127 browser genes. In total 17,660 genes match an OMIM entry, including entries without a phenotype relationship. Matching uses gene MIM IDs, with a unique Entrez ID plus a matching symbol or alias as fallback. Unmatched loci are recorded in `data_sources/omim-join-audit.json`; ClinGen assertions are not relabelled as OMIM records.

To populate the table from an authorized OMIM export:

```sh
python3 scripts/build_omim_relationships.py --genemap2 /path/to/genemap2.txt --source-date YYYY-MM-DD
python3 tests/test_omim.py
```

The importer joins by gene MIM ID and preserves disease MIM IDs and OMIM qualifiers separately. Source dates remain in Data sources; panel subtitles are omitted throughout the interface.

## Weekly clinical refresh

The Pages workflow refreshes clinical data each **Tuesday at 10:23 UTC** (2:23 am Pacific standard time, 3:23 am Pacific daylight time). GitHub may delay scheduled jobs. The schedule becomes active once the workflow is pushed to the repository's default `main` branch. **Actions > Deploy gene browser > Run workflow** also refreshes by default; uncheck **Refresh ClinGen and ClinVar before publishing** to deploy the committed snapshot only.

The refresh downloads the official ClinGen gene–disease validity export, ClinVar gene summary, and complete weekly ClinVar VCV XML. It verifies the XML checksum, rejects incomplete/stale inputs, rebuilds current dominant P/LP counts, then runs the tests and hosting build. Only successful validation permits a snapshot commit and deployment. A refresh failure leaves the live site on its last successful release; check the failed Actions run for details. GitHub can disable schedules on public repositories after 60 days without activity; re-enable the workflow in Actions if needed.

Current dominant counts are a **mapped subset**, not a complete catalog of all dominant mutations. Each distinct single-gene AlleleID must have P/LP germline classifications with assertion criteria at both variant and variant–condition levels, plus dominant inheritance linked to the same condition. Inheritance comes from qualifying ClinVar submissions or exact condition-ID matches to ClinGen (Moderate/Strong/Definitive) or OMIM (mapping key 3). Mixed and missing inheritance are reported separately. Multigene structural alleles, compound haplotype/genotype records, conflicting classifications and assertions without criteria do not qualify. Counts can decrease after reclassification. The **Dominant assessment** download records the included/unresolved variant–condition rows and their evidence; an allele can have more than one condition row but is counted once per gene.

OMIM stays at the supplied export date until a new authorized export is imported. Dosage sensitivity, editing candidates, cell lines, IMPC and the saved 1000 Genomes analysis are not refreshed by this job. Original D&D counts remain under Additional counts. Population plots default to a mutation-specific comparison using the current mapped dominant count; choose Original analysis to restore the earlier comparison. Only this equal-frequency model is recalculated: the measured cohort coverage, selected guides and original plots stay unchanged. All source dates and count definitions are available under Data sources.

To run the same refresh locally (several GB of downloads and temporary disk space):

```sh
python3 -m pip install matplotlib==3.11.2
python3 scripts/refresh_weekly.py
python3 tests/test_dominant.py
```

The complete VCV XML stays in the temporary work directory and is not committed or served. The repository stores only the relevant per-gene assessments and source hashes. Workflow commits use the repository's `GITHUB_TOKEN`; no personal token is needed. Repository rules must permit these snapshot commits to `main`. Because successful refreshes commit updated snapshots, fetch and pull `main` before pushing later local edits.
