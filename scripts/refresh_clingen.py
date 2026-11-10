"""Refresh all gene-disease validity assertions from ClinGen's official CSV.

python3 scripts/refresh_clingen.py --fetch
python3 scripts/refresh_clingen.py
Other evidence sources, including dosage sensitivity, are unchanged.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import re
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
CACHE = PROJECT / 'data_sources/clingen-gene-validity.csv'
OUT = PROJECT / 'dist/data'
URL = 'https://search.clinicalgenome.org/kb/gene-validity/download'


def read_validity(path):
    text = path.read_text(encoding='utf-8-sig')
    match = re.search(r'FILE CREATED: (\d{4}-\d{2}-\d{2})', text)
    if not match:
        raise ValueError('ClinGen export has no snapshot date')
    lines = text.splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('"GENE SYMBOL"'))
    reader = csv.DictReader(io.StringIO('\n'.join(lines[start:])))
    required = {'GENE SYMBOL','GENE ID (HGNC)','DISEASE LABEL','DISEASE ID (MONDO)',
                'MOI','CLASSIFICATION','ONLINE REPORT','CLASSIFICATION DATE','GCEP'}
    if not required.issubset(reader.fieldnames):
        raise ValueError('ClinGen export schema changed')
    rows = [{k: v.strip() or None for k,v in row.items()}
            for row in reader if row['GENE ID (HGNC)'].startswith('HGNC:')]
    if len(rows) < 3000 or any(not row['ONLINE REPORT'] for row in rows):
        raise ValueError('Incomplete ClinGen export')
    if len({row['ONLINE REPORT'] for row in rows}) != len(rows):
        raise ValueError('Duplicate ClinGen assertion URLs')
    return match[1], rows


def source_metadata(path, date, rows):
    return dict(id='validity', title='ClinGen gene–disease validity',
                path=str(path.relative_to(PROJECT)), url=URL, date=date,
                retrieved_at=datetime.now(timezone.utc).isoformat(),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                bytes=path.stat().st_size, records=len(rows),
                note='Disease, inheritance, classification, expert panel and evaluation date from the current complete export.')


def save(path, value):
    path.write_text(json.dumps(value, separators=(',', ':'), allow_nan=False))


def refresh():
    date, rows = read_validity(CACHE)
    shards = {p: json.loads(p.read_text()) for p in OUT.glob('genes-*.json')}
    genes = {s:g for data in shards.values() for s,g in data.items()}
    by_id = {g['hgnc_id']:g for g in genes.values()}
    unmatched = sorted({r['GENE ID (HGNC)'] for r in rows} - by_id.keys())
    if unmatched:
        raise ValueError(f'Unmapped ClinGen HGNC IDs; update HGNC before publishing: {unmatched}')
    grouped = defaultdict(list)
    for row in rows:
        grouped[row['GENE ID (HGNC)']].append(row)
    changed = sum(g['diseases'] != grouped[g['hgnc_id']] for g in genes.values())
    for gene in genes.values():
        gene['diseases'] = grouped[gene['hgnc_id']]
        gene['cmt_assertions'] = [r for r in gene['diseases']
                                 if 'charcot-marie-tooth' in r['DISEASE LABEL'].lower()]
    index = json.loads((OUT/'index.json').read_text())
    symbol_col, cmt_col = (index['fields'].index(k) for k in ('symbol','cmt'))
    for row in index['rows']:
        row[cmt_col] = bool(genes[row[symbol_col]]['cmt_assertions'])
    manifest = json.loads((OUT/'manifest.json').read_text())
    manifest['sources'] = [source_metadata(CACHE,date,rows) if s['id']=='validity' else s
                           for s in manifest['sources']]
    manifest['counts']['clingen'] = sum(bool(g['diseases']) for g in genes.values())
    manifest['audit']['unmatched']['validity'] = []
    manifest['clingen_refreshed_at'] = datetime.now(timezone.utc).isoformat()
    assert sum(len(g['diseases']) for g in genes.values()) == len(rows)
    for path, data in shards.items():
        save(path,data)
    save(OUT/'index.json',index)
    save(OUT/'manifest.json',manifest)
    print(f'ClinGen {date}: {len(rows)} assertions across {manifest['counts']['clingen']} genes; {changed} gene records changed; no unmatched HGNC IDs.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fetch', action='store_true')
    args = parser.parse_args()
    if args.fetch:
        temporary = CACHE.with_suffix('.download')
        try:
            subprocess.run(['curl','--fail','--silent','--show-error','--location',
                            '--connect-timeout','15','--max-time','120',URL,'-o',str(temporary)],check=True)
            read_validity(temporary)
            temporary.replace(CACHE)
        finally:
            temporary.unlink(missing_ok=True)
    refresh()
