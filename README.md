# D&D Gene Browser

An interactive gene browser for clinical annotations, dosage sensitivity, evolutionary constraint, variant records, targeting results and mouse phenotype evidence.

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
