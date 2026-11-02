"""Prepare gene-specific CSV downloads and cell-line summaries without changing inputs.

Run after build_data.py. ClinVar is streamed into temporary alphabetic partitions
so the complete variant file does not need to fit in memory.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import tempfile
from collections import defaultdict
from pathlib import Path
from build_strategy_screen import build_screen
from pack_clinvar import pack_clinvar

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT.parent
OUT = PROJECT / 'dist'
LIBRARY = ROOT / 'Grace_pipeline/261011_RUN_final/mxy_local_eda/outputs/20261101_library_update'
CLINVAR = ROOT / 'DnD_project/haplo_prelim_pipeline/data/ClinVar/variant_summary.txt'


def bucket(symbol):
    return symbol[0] if symbol[0].isascii() and symbol[0].isalpha() else '_'


def safe_cell(value):
    """Keep imported text from becoming a spreadsheet formula."""
    text = '' if value is None else str(value)
    return "'" + text if text.startswith(('=', '+', '-', '@', '\t', '\r')) else text


def write_csv(relative, headers, rows):
    path = OUT / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, 'wt', encoding='utf-8-sig', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(headers)
        for row in rows:
            writer.writerow([safe_cell(v) for v in row])


def read_csv(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def source(path, key, title, note, date=None, digest=None):
    if digest is None:
        hasher = hashlib.sha256()
        with path.open('rb') as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                hasher.update(chunk)
        digest = hasher.hexdigest()
    return dict(id=key, path=str(path.relative_to(ROOT)), title=title, note=note,
                date=date, sha256=digest,
                bytes=path.stat().st_size)


def selection(rows, guides):
    selected = [r for r in rows if r['library_candidate_postNGG'] == 'True']
    associations = {r['candidate_association_id'] for r in selected}
    chosen = [r for r in guides if r['candidate_association_id'] in associations]
    return dict(variants=len({r['variant_id'] for r in rows}),
                ngg_variants=len({r['variant_id'] for r in rows if int(float(r['postNGG_unique_guides'])) > 0}),
                selected_variants=len({r['variant_id'] for r in selected}),
                guides=len({r['guide_id'] for r in chosen}),
                editors=sorted({r['editor'] for r in chosen}))


def build_targeting(genes, details):
    # The source uses uppercase keys for C9ORF72, while HGNC uses C9orf72.
    symbols = {s.upper(): s for s in genes}
    variants = defaultdict(list)
    guides = defaultdict(list)
    variant_path = LIBRARY / 'snp_list_prePAM.csv'
    guide_path = LIBRARY / 'grna_detail_postNGG.csv.gz'
    for row in read_csv(variant_path):
        if row['gene'].upper() in symbols:
            variants[symbols[row['gene'].upper()]].append(row)
    for row in read_csv(guide_path):
        if row['gene'].upper() in symbols:
            guides[symbols[row['gene'].upper()]].append(row)
    for symbol, rows in variants.items():
        guide_rows = guides[symbol]
        scopes = []
        for window in ('1-23', '1-30'):
            subset = [r for r in rows if r['window_set'] in (window, 'exonic')]
            if not subset:
                continue
            line_rows = [dict(cell_line='All candidates', **selection(subset, guide_rows))]
            for line, flag in [('WTC', 'het_WTC'), ('WTD', 'het_WTD'), ('p47', 'het_p47_locus_only')]:
                line_rows.append(dict(cell_line=line, **selection([r for r in subset if r[flag] == 'True'], guide_rows)))
            scopes.append(dict(window=window, rows=line_rows))
        exports = {}
        for kind, data in [('variants', rows), ('guides', guide_rows)]:
            if not data:
                continue
            relative = f'downloads/targeting/{symbol}-{kind}.csv.gz'
            headers = list(data[0])
            write_csv(relative, headers, ([r[k] for k in headers] for r in data))
            exports[kind] = dict(path=relative, rows=len(data),
                                 unique=len({r['variant_id' if kind == 'variants' else 'guide_id'] for r in data}))
        details[symbol]['targeting'] = dict(scopes=scopes, exports=exports)
    return [source(variant_path, 'targeting_variants', 'Targeting SNP records',
                   'GRCh38 common SNP candidates, cell-line genotype evidence and selection flags. Each row is a variant/editor/window association. Recovered source coverage is incomplete.', '2026-11-01'),
            source(guide_path, 'targeting_guides', 'gRNA sequences and target associations',
                   'NGG-compatible guides on GRCh38 plus the focal allele; neighboring variants are unphased. Guide/editor/window associations can repeat a guide. Selection and provisional editing-window flags are preserved.', '2026-11-01')]


def build_clinvar(genes, details):
    by_id = {str(g['entrez_id']): s for s, g in genes.items() if g.get('entrez_id')}
    matched = 0
    digest = hashlib.sha256()
    with tempfile.TemporaryDirectory(prefix='dnd-clinvar-') as directory:
        partitions = {}
        with CLINVAR.open('rb') as stream:
            header = stream.readline()
            digest.update(header)
            headers = header.decode('utf-8-sig').rstrip('\r\n').split('\t')
            gene_col, assembly_col = headers.index('GeneID'), headers.index('Assembly')
            for raw in stream:
                digest.update(raw)
                line = raw.decode('utf-8').rstrip('\r\n')
                if not line:
                    continue
                values = line.split('\t')
                if values[assembly_col] != 'GRCh38':
                    continue
                symbols = {by_id[x.strip()] for x in re.split(r'[;,|]', values[gene_col]) if x.strip() in by_id}
                for symbol in symbols:
                    key = bucket(symbol)
                    if key not in partitions:
                        partitions[key] = (Path(directory) / key).open('w')
                    partitions[key].write(symbol + '\t' + line + '\n')
                    matched += 1
        for stream in partitions.values():
            stream.close()
        print(f'ClinVar: {matched:,} gene-record associations; preparing downloads.', flush=True)
        for key in sorted(partitions):
            groups = defaultdict(list)
            with (Path(directory) / key).open() as stream:
                for line in stream:
                    symbol, record = line.rstrip('\n').split('\t', 1)
                    groups[symbol].append(record)
            for symbol, records in groups.items():
                relative = f'downloads/clinvar/{key}/{symbol}.csv.gz'
                write_csv(relative, headers, (r.split('\t') for r in records))
                details[symbol]['clinvar_export'] = dict(path=relative, rows=len(records),
                                                        alleles=len({r.split('\t', 1)[0] for r in records}))
            print(f'ClinVar {key}: {len(groups):,} genes.', flush=True)
    return source(CLINVAR, 'clinvar_variants', 'ClinVar variant records (GRCh38)',
                  'Local variant_summary.txt; snapshot date is not recorded. Downloads retain all source columns for GRCh38 records matched by NCBI Gene ID. These records may differ from the dated gene-summary totals. No variant-level dominant inheritance assignment is supplied.',
                  digest=digest.hexdigest())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--skip-clinvar', action='store_true', help='Reuse existing ClinVar downloads while refreshing targeting data.')
    args = parser.parse_args()
    genes = {}
    for path in (OUT / 'data').glob('genes-*.json'):
        genes.update(json.loads(path.read_text()))
    details = defaultdict(dict)
    sources = []
    if args.skip_clinvar:
        for path in (OUT / 'data').glob('record-details-*.json'):
            for symbol, detail in json.loads(path.read_text()).items():
                if 'clinvar_export' in detail:
                    details[symbol]['clinvar_export'] = detail['clinvar_export']
        previous = json.loads((OUT / 'data/record-sources.json').read_text())
        sources = [s for s in previous['sources'] if s['id'] == 'clinvar_variants']
    else:
        sources.append(build_clinvar(genes, details))
    sources.extend(build_targeting(genes, details))
    sources.extend(build_screen(ROOT, OUT, genes, details, source, write_csv))
    coordinate_path = PROJECT / 'data_sources/gene_coordinates.json'
    cached = json.loads(coordinate_path.read_text())
    for row in cached['records']:
        symbol = row['display_name']
        if symbol not in genes or genes[symbol]['ensembl_id'] != row['id'] or row['assembly_name'] != 'GRCh38':
            raise ValueError(f'Coordinate identity or assembly mismatch: {symbol}')
        if not genes[symbol].get('coordinates'):
            details[symbol]['coordinates'] = f"{row['seq_region_name']}:{row['start']}-{row['end']}"
            details[symbol]['coordinate_source'] = 'ensembl_coordinates'
    coordinate_source = source(coordinate_path, 'ensembl_coordinates', 'Ensembl gene coordinates',
                               'GRCh38 coordinates for SPTLC1, checked against the Ensembl gene ID. Other genes use ClinGen coordinates where supplied, or a UCSC canonical-gene search.', cached['retrieved'])
    coordinate_source['url'] = cached['source']
    sources.append(coordinate_source)
    shards = {bucket(symbol): {} for symbol in genes}
    for symbol, detail in details.items():
        shards[bucket(symbol)][symbol] = detail
    for key, data in shards.items():
        (OUT / f'data/record-details-{key}.json').write_text(json.dumps(data, separators=(',', ':')))
    (OUT / 'data/record-sources.json').write_text(json.dumps({'sources': sources}, indent=2))
    pack_clinvar()
    print(f'Downloads indexed for {len(details):,} genes; targeting summaries for {sum("targeting" in d for d in details.values())} genes.', flush=True)


if __name__ == '__main__':
    main()
