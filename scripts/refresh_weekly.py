"""Refresh official ClinGen validity and ClinVar counts for the weekly Pages build.

A failed refresh exits nonzero; the Actions job must validate before committing
or deploying. OMIM, dosage, targeting and population-cohort inputs retain their
own snapshot dates. No paper-derived counts are relabeled as current.
"""
from __future__ import annotations
import argparse
import csv
import json
import subprocess
import sys
import tempfile
from collections import defaultdict
from datetime import date, datetime, timezone
from pathlib import Path

import refresh_clingen
import refresh_dominant
import refresh_evidence

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'dist/data'
SUMMARY_URL = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/gene_specific_summary.txt'


def check_date(value, maximum_age=21):
    age = (date.today() - date.fromisoformat(value)).days
    if not -1 <= age <= maximum_age:
        raise ValueError(f'Unexpected source date {value}; age {age} days')


def validate_summary(path):
    with path.open() as stream:
        release = datetime.strptime(next(stream).strip().split('dated ')[1], '%B %d, %Y').date().isoformat()
        reader = csv.DictReader(stream, delimiter='\t')
        required = {'GeneID', 'Total_alleles', 'Alleles_reported_Pathogenic_Likely_pathogenic'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('ClinVar gene-summary schema changed')
        rows = list(reader)
    if len(rows) < 10000:
        raise ValueError('Incomplete ClinVar gene summary')
    for row in rows:
        for field in required - {'GeneID'}:
            if row[field] != '-' and not row[field].isdigit():
                raise ValueError(f'Invalid ClinVar count: {field}')
    check_date(release)
    return release


def refresh_summary(work_dir):
    temporary = work_dir / 'gene-summary.txt'
    refresh_dominant.fetch(SUMMARY_URL, temporary)
    release = validate_summary(temporary)
    destination = ROOT / 'data_sources/comprehensive/clinvar_gene_summary.txt'
    destination.write_bytes(temporary.read_bytes())
    genes = {}
    for path in DATA.glob('genes-*.json'):
        genes.update(json.loads(path.read_text()))
    shards = {path: json.loads(path.read_text()) for path in DATA.glob('evidence-?.json')}
    evidence = defaultdict(dict, {symbol: row for shard in shards.values() for symbol, row in shard.items()})
    for row in evidence.values():
        for field in ('clinvar', 'clinvar_records', 'clinvar_date'):
            row.pop(field, None)
    # Clear the old snapshot first so a removed source record cannot survive.
    index_path = DATA / 'index.json'
    index = json.loads(index_path.read_text())
    column = index['fields'].index('clinvar_total')
    for row in index['rows']:
        row[column] = None
    refresh_dominant.save(index_path, index)
    source = refresh_evidence.build_clinvar(genes, evidence)[0]
    source['title'] = 'ClinVar gene summary'
    source['note'] = ('Distinct alleles in or overlapping the gene; P/LP excludes overlapping structural variants. '
                      'Identical duplicate OMIM gene rows are not summed. Variant downloads and the D&D reference counts retain their separate dates.')
    audit = json.loads((DATA / 'clinvar-current-audit.json').read_text())
    if audit['conflicts']:
        raise ValueError(f"Conflicting gene-summary records: {audit['conflicts']}")
    for path in shards:
        letter = path.stem.removeprefix('evidence-')
        refresh_dominant.save(path, {s: row for s, row in evidence.items() if (s[0] if s[0].isalpha() else '_') == letter})
    path = DATA / 'evidence-sources.json'
    sources = json.loads(path.read_text())
    sources['sources'] = [s for s in sources['sources'] if s['id'] != 'clinvar_current'] + [source]
    refresh_dominant.save(path, sources)
    print(f'ClinVar summary {release}: {audit["mapped_genes"]:,} genes', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, default=Path(tempfile.gettempdir()) / 'dnd-weekly')
    parser.add_argument('--xml', type=Path, help='Use an already downloaded, checksum-verified weekly VCV XML')
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    previous = json.loads((DATA / 'manifest.json').read_text())
    previous_dates = {s['id']: s.get('date') for s in previous['sources']}
    subprocess.run([sys.executable, str(ROOT / 'scripts/refresh_clingen.py'), '--fetch'], check=True)
    check_date(refresh_clingen.read_validity(refresh_clingen.CACHE)[0])
    refresh_summary(args.work_dir)
    command = [sys.executable, str(ROOT / 'scripts/refresh_dominant.py')]
    command += ['--xml', str(args.xml)] if args.xml else ['--fetch', '--work-dir', str(args.work_dir)]
    subprocess.run(command, check=True)
    manifest = json.loads((DATA / 'manifest.json').read_text())
    source = next(s for s in manifest['sources'] if s['id'] == 'dominant_current')
    check_date(source['date'])
    for refreshed in manifest['sources']:
        if refreshed['id'] in {'validity', 'dominant_current'} and previous_dates.get(refreshed['id']):
            if refreshed['date'] < previous_dates[refreshed['id']]:
                raise ValueError(f"Source release moved backwards: {refreshed['id']}")
    subprocess.run([sys.executable, str(ROOT / 'scripts/refresh_population_comparison.py')], check=True)
    manifest['weekly_refresh'] = dict(checked_at=datetime.now(timezone.utc).isoformat(), schedule_utc='23 10 * * 2',
                                    sources=['validity', 'clinvar_current', 'dominant_current'])
    refresh_dominant.save(DATA / 'manifest.json', manifest)
    print('Weekly refresh complete. Validate before publishing.', flush=True)


if __name__ == '__main__':
    main()
