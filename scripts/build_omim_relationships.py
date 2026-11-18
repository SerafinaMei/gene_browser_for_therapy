"""Build OMIM relationship shards from the cached records, or import an authorized genemap2 export.

python3 scripts/build_omim_relationships.py [--genemap2 PATH --source-date YYYY-MM-DD]
Gene MIM numbers join HGNC gene entries; phenotype MIM numbers remain distinct.
"""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CACHE=ROOT/'data_sources/omim-gene-phenotypes.json'
OUT=ROOT/'dist/data'

def parse_phenotype(text):
    # Keep OMIM's ?, {}, and [] qualifiers intact, including descriptions with commas.
    match=re.fullmatch(r'(.*?),\s*(\d{6})\s*\(([1-4])\)(?:,\s*(.*))?',text.strip())
    if match:
        description,mim,key,inheritance=match.groups()
    else:
        match=re.fullmatch(r'(.*?)\s*\(([1-4])\)(?:,\s*(.*))?',text.strip())
        if not match:raise ValueError(f'Unrecognized OMIM phenotype field: {text!r}')
        description,key,inheritance=match.groups();mim=None
    return dict(phenotype=description,phenotype_mim=mim,mapping_key=int(key),inheritance=inheritance or None)

def import_genemap(path,date):
    if not re.fullmatch(r'\d{4}-\d{2}-\d{2}',date):raise ValueError('Supply the export source date as YYYY-MM-DD')
    header=None;records={}
    with path.open(encoding='utf-8-sig') as f:
        for line in f:
            if line.startswith('# Chromosome\t'):
                header=line.lstrip('# ').rstrip('\n').split('\t')
            elif line.startswith('#') or not line.strip():continue
            else:
                if header is None:raise ValueError('Missing genemap2 column header')
                row=dict(zip(header,next(csv.reader([line],delimiter='\t'))))
                gene_mim=row['Mim Number']
                relationships=[]
                for value in row['Phenotypes'].split(';'):
                    if not value.strip():continue
                    record=parse_phenotype(value);record['location']=row['Cyto Location'] or row['Computed Cyto Location'] or None
                    relationships.append(record)
                if gene_mim in records and records[gene_mim]!=relationships:raise ValueError(f'Conflicting records for {gene_mim}')
                records[gene_mim]=relationships
    if header is None:raise ValueError('No genemap2 data found')
    return dict(source=dict(kind='OMIM genemap2 export',date=date,url='https://omim.org/downloads',sha256=hashlib.sha256(path.read_bytes()).hexdigest(),scope='Imported genemap2 records. Coverage follows the supplied export.'),relationships=records)

def build():
    data=json.loads(CACHE.read_text());count=0;rows=0
    for path in OUT.glob('genes-*.json'):
        shard={}
        for symbol,g in json.loads(path.read_text()).items():
            matches=[(mim,data['relationships'][mim]) for mim in g['omim_ids'] if mim in data['relationships']]
            if not matches:continue
            relationships=[dict(r,gene_mim=mim) for mim,records in matches for r in records]
            shard[symbol]={'omim_relationships':dict(records=relationships,source=data['source'])}
            count+=1;rows+=len(relationships)
        (OUT/path.name.replace('genes-','omim-')).write_text(json.dumps(shard,separators=(',',':')))
    manifest=json.loads((OUT/'manifest.json').read_text())
    manifest['sources']=[s for s in manifest['sources'] if s['id']!='omim_relationships']+[dict(id='omim_relationships',title='OMIM gene-phenotype relationships',url=data['source']['url'],date=data['source'].get('date') or data['source'].get('retrieved'),path=str(CACHE.relative_to(ROOT)),sha256=hashlib.sha256(CACHE.read_bytes()).hexdigest(),note=data['source']['kind']+'. '+data['source']['scope'])]
    (OUT/'manifest.json').write_text(json.dumps(manifest,separators=(',',':')))
    print(f'OMIM: {rows} phenotype relationships for {count} genes.')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--genemap2',type=Path);parser.add_argument('--source-date');args=parser.parse_args()
    if args.genemap2:
        if not args.source_date:parser.error('--source-date required with --genemap2')
        CACHE.write_text(json.dumps(import_genemap(args.genemap2,args.source_date),separators=(',',':')))
    build()
