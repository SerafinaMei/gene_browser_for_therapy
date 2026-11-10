# Customize the website

The website is plain HTML, CSS and JavaScript. No framework or paid editor is required.

## Appearance

The website uses the Shiny-style layout for all visitors, without an appearance selector. Previously saved browser preferences do not override it.

The sidebar, compact tables and muted colors are defined in [dist/shiny.css](dist/shiny.css); [dist/appearance.js](dist/appearance.js) places the search field in the main area on the home page and in the sidebar on reports.

The complete original Classic version is saved under Git tag `style-classic-2026-11-10`. The version with both selectable appearances is also preserved in the commit history.

## Text, colors and size

Edit [dist/site-config.json](dist/site-config.json). On GitHub, open that file and select the pencil icon. On your Mac, open it in your editor or double-click `Customize Website.command`.

| Setting | Controls |
| --- | --- |
| `title` | Header name and browser-tab titles |
| `organization` | Lab or institution in the footer |
| `status_label` | Small status label, such as `Test version`; use `""` to hide it |
| `footer_note` | Footer description |
| `search_heading`, `description` | Home-page heading and description |
| `search_label`, `search_placeholder`, `search_help` | Search field wording |
| `example_genes` | Example gene buttons; use `[]` to hide them |
| `font_family` | `system`, `arial` or `verdana` |
| `body_font_size` | Body text size, from 12 to 20 pixels |
| `content_width` | Maximum page width, from 900 to 1600 pixels |
| `theme` | Text, muted text, borders, surfaces, links, page background and button colors |

Use six-digit hex colors, such as `#285b88`. Keep double quotes around text and commas between entries. The automated checks reject invalid settings before deployment. Text is displayed as text, not interpreted as HTML.

## Layout and detailed styling

Add CSS rules to [dist/custom.css](dist/custom.css). It loads after the main stylesheet. Shiny-specific rules load afterward from `dist/shiny.css`. Examples are included in the file. The page structure is in `dist/index.html`; report rendering is in `dist/app.js`.

## Preview and publish

1. Run `python3 scripts/serve.py` and open `http://127.0.0.1:8765/`.
2. Save your changes and refresh the browser.
3. Commit and push the changes to `main`. The GitHub Actions workflow checks and publishes the website automatically once Pages is enabled.

You can also edit `dist/site-config.json` or `dist/custom.css` directly on GitHub and commit the change there. Watch **Actions > Deploy gene browser** for completion. The address stays the same. Revert a commit to undo an update.

The site remains a public test version. Its HTML asks search engines not to index it; this does not prevent someone from opening or forwarding its address. Keep the status label until you choose to release it.

## Gene data

Appearance changes do not require rebuilding the scientific datasets. The checked-in `dist/data/` and `dist/downloads/` folders contain the browser's data. Rebuilding them requires the research inputs described in `SCIENTIFIC_NOTES.md`; they are not recreated during ordinary GitHub deployments.
