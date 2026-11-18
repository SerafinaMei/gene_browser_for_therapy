"""Build OMIM relationship shards from the cached records, or import an authorized genemap2 export.

python3 scripts/build_omim_relationships.py [--genemap2 PATH --source-date YYYY-MM-DD]
Gene MIM numbers join HGNC gene entries; phenotype MIM numbers remain distinct.
"""
import argparse
import csv
from collections import defaultdict
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
    header=None;records={};identifiers={}
    with path.open(encoding='utf-8-sig') as f:
        for line in f:
            if line.startswith('# Generated: '):
                if line.strip().removeprefix('# Generated: ')!=date:raise ValueError('Source date differs from the export header')
                continue
            if line.startswith('# Chromosome\t'):
                header=line.lstrip('# ').rstrip('\n').split('\t')
            elif line.startswith('#') or not line.strip():continue
            else:
                if header is None:raise ValueError('Missing genemap2 column header')
                row=dict(zip(header,next(csv.reader([line],delimiter='\t'))))
                gene_mim=row.get('MIM Number') or row.get('Mim Number')
                if not gene_mim or not re.fullmatch(r'\d{6}',gene_mim):raise ValueError('Missing or invalid gene MIM number')
                relationships=[]
                for value in row['Phenotypes'].split(';'):
                    if not value.strip():continue
                    record=parse_phenotype(value);record['location']=row['Cyto Location'] or row['Computed Cyto Location'] or None
                    relationships.append(record)
                if gene_mim in records and records[gene_mim]!=relationships:raise ValueError(f'Conflicting records for {gene_mim}')
                records[gene_mim]=relationships
                identifiers[gene_mim]=dict(symbol=row.get('Approved Gene Symbol',''),entrez_id=row.get('Entrez Gene ID',''),ensembl_id=row.get('Ensembl Gene ID',''))
    if header is None:raise ValueError('No genemap2 data found')
    return dict(source=dict(kind='OMIM genemap2 export',date=date,url='https://omim.org/downloads',sha256=hashlib.sha256(path.read_bytes()).hexdigest(),scope='Imported genemap2 records. Coverage follows the supplied export.'),relationships=records,gene_identifiers=identifiers)

def match_genes(data,genes):
    matches={symbol:[mim for mim in g['omim_ids'] if mim in data['relationships']] for symbol,g in genes.items()}
    used={mim for ids in matches.values() for mim in ids}
    by_entrez=defaultdict(list)
    for symbol,g in genes.items():
        if g.get('entrez_id'):by_entrez[g['entrez_id']].append(symbol)
    recovered=[];conflicts=[]
    for mim,ids in data.get('gene_identifiers',{}).items():
        if mim in used or not ids['symbol'] or not ids['entrez_id']:continue
        candidates=by_entrez[ids['entrez_id']]
        if len(candidates)!=1:continue
        symbol=candidates[0]
        if ids['symbol'] not in [symbol,*genes[symbol].get('aliases',[])]:
            conflicts.append(dict(gene_mim=mim,source_symbol=ids['symbol'],browser_symbol=symbol));continue
        matches[symbol].append(mim);used.add(mim)
        recovered.append(dict(symbol=symbol,gene_mim=mim,entrez_id=ids['entrez_id']))
    audit=dict(recovered_by_entrez=recovered,conflicts=conflicts,
               unmatched_phenotype_mims=sorted(mim for mim,rows in data['relationships'].items() if rows and mim not in used))
    return matches,audit

def build():
    data=json.loads(CACHE.read_text());count=0;rows=0;with_phenotypes=0
    genes={}
    for path in OUT.glob('genes-*.json'):genes.update(json.loads(path.read_text()))
    matches,audit=match_genes(data,genes)
    for path in OUT.glob('genes-*.json'):
        shard={}
        for symbol in json.loads(path.read_text()):
            ids=matches[symbol]
            if not ids:continue
            relationships=[dict(r,gene_mim=mim) for mim in ids for r in data['relationships'][mim]]
            shard[symbol]={'omim_relationships':dict(records=relationships,gene_mim_ids=ids,source=data['source'])}
            count+=1;rows+=len(relationships);with_phenotypes+=bool(relationships)
        (OUT/path.name.replace('genes-','omim-')).write_text(json.dumps(shard,separators=(',',':')))
    audit.update(matched_genes=count,genes_with_phenotypes=with_phenotypes,phenotype_relationships=rows)
    (ROOT/'data_sources/omim-join-audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    manifest=json.loads((OUT/'manifest.json').read_text())
    manifest['sources']=[s for s in manifest['sources'] if s['id']!='omim_relationships']+[dict(id='omim_relationships',title='OMIM gene-phenotype relationships',url=data['source']['url'],date=data['source'].get('date') or data['source'].get('retrieved'),path=str(CACHE.relative_to(ROOT)),sha256=hashlib.sha256(CACHE.read_bytes()).hexdigest(),note=data['source']['kind']+'. '+data['source']['scope'])]
    (OUT/'manifest.json').write_text(json.dumps(manifest,separators=(',',':')))
    print(f'OMIM: {rows} phenotype relationships for {with_phenotypes} genes; {count} genes matched overall; {len(audit["recovered_by_entrez"])} recovered by Entrez ID.')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--genemap2',type=Path);parser.add_argument('--source-date');args=parser.parse_args()
    if args.genemap2:
        if not args.source_date:parser.error('--source-date required with --genemap2')
        CACHE.write_text(json.dumps(import_genemap(args.genemap2,args.source_date),separators=(',',':')))
    build()
