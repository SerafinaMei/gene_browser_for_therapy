"""Import genome-wide UCSC targets and current IMPC evidence without altering inputs.

python3 scripts/refresh_evidence.py --fetch-mouse
python3 scripts/refresh_evidence.py --build
Optional inputs are configured with --track-bed, --cell-dir and --bigbed-tool.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, io, json, os, re, subprocess, tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen
from urllib.parse import urlencode

PROJECT = Path(__file__).resolve().parents[1]
CACHE = PROJECT / 'data_sources/comprehensive'
OUT = PROJECT / 'dist'
BASE = 'https://www.ebi.ac.uk/mi/impc/solr/'

def digest(path):
    with open(path, 'rb') as f: return hashlib.file_digest(f, 'sha256').hexdigest()

def save_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, separators=(',', ':')))

def fetch_mouse(resume=False):
    queries = {
        'gene': ('*:*', 'mgi_accession_id,marker_symbol,human_gene_symbol,phenotype_status,significant_top_level_mp_terms,not_significant_top_level_mp_terms', 'mgi_accession_id asc'),
        'genotype-phenotype': ('*:*', 'doc_id,marker_accession_id,marker_symbol,mp_term_id,mp_term_name,top_level_mp_term_name,allele_accession_id,allele_symbol,zygosity,life_stage_name,sex,p_value,assertion_type,statistical_method,procedure_name,parameter_name,phenotyping_center,colony_id', 'doc_id asc'),
        'diseases': ('type:disease_gene_summary', 'disease_id,disease_term,hgnc_gene_id,hgnc_gene_symbol,marker_id,marker_symbol', 'marker_id asc,disease_id asc'),
        'disease-models': ('type:disease_model_summary AND association_curated:true', 'disease_id,marker_id,model_id,model_source,model_description,disease_matched_phenotypes,model_matched_phenotypes,disease_model_max_norm,association_curated', 'marker_id asc,disease_id asc,model_id asc'),
    }
    metadata = {'retrieved_at': datetime.now(timezone.utc).isoformat(), 'release': '24.0', 'release_url': 'https://www.mousephenotype.org/data/release', 'sources': []}
    for name, (query, fields, sort) in queries.items():
        core = 'phenodigm' if name.startswith('disease') else name
        cache = CACHE / (name + '.json.gz')
        docs, expected = (read_cache(name), len(read_cache(name))) if resume and cache.exists() else ([], None)
        while expected is None or len(docs) < expected:
            url = BASE + core + '/select?' + urlencode(dict(q=query, fl=fields, sort=sort, rows=10000, start=len(docs), wt='json'))
            response = subprocess.run(['curl','--compressed','--fail','--silent','--show-error','--retry','3','--retry-all-errors','--max-time','120',url],check=True,capture_output=True)
            data = json.loads(response.stdout)['response']
            if expected is not None and data['numFound'] != expected: raise ValueError('Source changed during pagination; retry snapshot')
            expected = data['numFound']
            if not data['docs'] and len(docs) < expected: raise ValueError('Truncated source')
            docs.extend(data['docs'])
            print(name, len(docs), '/', expected, flush=True)
        path = CACHE / (name + '.json.gz')
        path.write_bytes(gzip.compress(json.dumps(docs, separators=(',', ':')).encode(), mtime=0))
        metadata['sources'].append(dict(id='impc_' + name, title='IMPC ' + name + ' · current API snapshot', url=BASE + core + '/select?' + urlencode(dict(q=query, wt='json')), date=metadata['retrieved_at'][:10], records=len(docs), sha256=digest(path), path=str(path.relative_to(PROJECT)), note='IMPC release 24.0. Complete paginated query; counts checked against numFound.'))
    save_json(CACHE / 'mouse-sources.json', metadata)

def read_cache(name):
    with gzip.open(CACHE / (name + '.json.gz'), 'rt') as f: return json.load(f)

def build_mouse(genes, evidence):
    by_hgnc = {v['hgnc_id']: k for k,v in genes.items()}
    previous={}
    with (PROJECT/'data_sources/hgnc_complete_set.txt').open() as f:
        for g in csv.DictReader(f,delimiter='\t'):
            if g['symbol'] in genes:
                for old in g.get('prev_symbol','').split('|'):
                    if old:previous[old]=g['symbol']
    mouse_to_human = defaultdict(set)
    for row in read_cache('gene'):
        mid = row['mgi_accession_id']
        for symbol in row.get('human_gene_symbol', []):
            if symbol not in genes:symbol=previous.get(symbol,symbol)
            if symbol not in genes: continue
            mouse_to_human[mid].add(symbol)
            m = evidence[symbol].setdefault('mouse_current', {'models': [], 'phenotypes': [], 'diseases': []})
            m['models'].append(row)
    groups = defaultdict(list)
    seen = set()
    for row in read_cache('genotype-phenotype'):
        assert row['doc_id'] not in seen, 'Repeated IMPC document'
        seen.add(row['doc_id'])
        key = tuple(row.get(k, '') for k in ['marker_accession_id','mp_term_id','allele_accession_id','zygosity','life_stage_name'])
        groups[key].append(row)
    for key, rows in groups.items():
        for symbol in mouse_to_human[key[0]]:
            evidence[symbol]['mouse_current']['phenotypes'].append({'term_id': key[1], 'term': rows[0].get('mp_term_name','Phenotype term not supplied: '+rows[0].get('parameter_name','Unspecified parameter')), 'mouse_id':key[0], 'allele':rows[0].get('allele_symbol'), 'zygosity':key[3], 'life_stage':key[4], 'support':rows})
    models = defaultdict(list)
    for row in read_cache('disease-models'): models[(row['marker_id'], row['disease_id'])].append(row)
    for row in read_cache('diseases'):
        symbol = by_hgnc.get(row.get('hgnc_gene_id'))
        if not symbol: continue
        m = evidence[symbol].setdefault('mouse_current', {'models': [], 'phenotypes': [], 'diseases': []})
        matches = models[(row.get('marker_id'),row['disease_id'])]
        row['matched_phenotypes'] = sorted({p for r in matches for p in r.get('disease_matched_phenotypes', [])})
        row['model_sources'] = sorted({r['model_source'] for r in matches})
        m['diseases'].append(row)
    for record in evidence.values():
        m = record.get('mouse_current')
        if not m: continue
        m['phenotypes'].sort(key=lambda p:(p['zygosity'],p['term'],p['life_stage']))
        m['unique_terms'] = len({p['term_id'] for p in m['phenotypes'] if p['term_id']})
        m['supporting_datasets'] = sum(len(p['support']) for p in m['phenotypes'])
    return json.loads((CACHE/'mouse-sources.json').read_text())['sources']

def build_tracks(args, genes, evidence, focused=False):
    # Preserve allele identity. A locus-only intersection is never a qualifying match.
    source_id='sptlc_tracks' if focused else 'genome_tracks'
    label='SPTLC1 / SPTLC2 published track' if focused else 'All-protein-coding UCSC track'
    rows_by_gene, loci = defaultdict(list), defaultdict(list)
    unmatched = defaultdict(int)
    previous = {}
    with (PROJECT/'data_sources/hgnc_complete_set.txt').open() as f:
        for g in csv.DictReader(f,delimiter='\t'):
            if g['symbol'] in genes:
                for old in g.get('prev_symbol','').split('|'):
                    if old: previous[old]=g['symbol']
    renamed={}
    source_rows = 0
    with args.track_bed.open() as f:
        for line in f:
            c = line.rstrip('\n').split('\t')
            if len(c) != 25: raise ValueError('Expected the verified 25-column all-protein-coding track')
            match = re.fullmatch(r'A therapeutically editable common variant for the D&D gene ([A-Za-z0-9_.-]+)',c[9])
            if not match: raise ValueError('Unexpected target gene description')
            source_rows += 1
            symbol = match[1]
            original_symbol=symbol
            if symbol not in genes and symbol in previous:
                symbol=previous[symbol];renamed[original_symbol]=symbol
            if symbol not in genes: unmatched[symbol] += 1; continue
            # chromosome, zero-based start/end, label, alleles, MAF, strategies, PAM flag, excision partners, cell matches
            r = [c[0],int(c[1]),int(c[2]),c[3].strip(),c[10],c[12],c[19],c[21],c[20],{},original_symbol]
            rows_by_gene[symbol].append(r)
            loci[(r[0],r[1],r[2])].append(r)
    source_list = [dict(id=source_id, title=label, path=str(args.track_bed), sha256=digest(args.track_bed), date=None, url=('https://raw.githubusercontent.com/gramey02/DnD_TrackHubs_Public/main/SPTLC_genes/SPTLC_genes.bb' if focused else 'https://genome.ucsc.edu/hubspace/10/gramey02/ProteinCoding_gene_editing/all_protein_coding_genes.bb'), note='Lossless BED export; upstream checked 2026-11-10. Source biological snapshot date unspecified. Distinct GRCh38 chromosome/start/end/allele keys; strategies overlap. PAM flags preserved, no guide selection inferred.')]
    cell_audit = {}
    env = dict(os.environ, DYLD_FALLBACK_LIBRARY_PATH=str(args.bigbed_tool.parent.parent/'lib'))
    for cell in ['WTB','WTC','WTD','KOLF2']:
        source = args.cell_dir / (cell+'_het.bb')
        matched, mismatched, total = set(), set(), 0
        with tempfile.TemporaryDirectory() as temp:
            bed = Path(temp)/'cell.bed'
            subprocess.run([str(args.bigbed_tool),str(source),str(bed)],env=env,check=True)
            with bed.open() as f:
                for line in f:
                    c=line.rstrip('\n').split('\t'); total += 1
                    for r in loci.get((c[0],int(c[1]),int(c[2])),[]):
                        # Cell name contains the observed allele pair, sometimes prefixed by rsID.
                        pair = re.search(r'([ACGT])/([ACGT])', c[3])
                        if not pair or sorted(pair.groups()) != sorted(r[4].split('/')):
                            mismatched.add((c[0],c[1],c[2],c[3])); continue
                        quality=c[-1]
                        if quality not in ['HIGH','LOW']: raise ValueError('Unrecognized cell-line quality')
                        r[9].setdefault(cell,set()).add(quality)
                        matched.add((c[0],c[1],c[2],r[4]))
        cell_audit[cell] = dict(records=total,matched_alleles=len(matched),unmatched_alleles_at_target_loci=len(mismatched))
        source_list.append(dict(id='track_'+cell,title=cell+' heterozygous genome track',path=str(source),sha256=digest(source),date=None,url='https://github.com/gramey02/DnD_TrackHubs_Public',note='Public upstream track retrieved and Git blob hash verified 2026-11-10. Source snapshot date unspecified. Allele-matched intersections; HIGH and LOW retained separately.'))
        print(cell,cell_audit[cell],flush=True)
    headers=['chromosome','start_0based','end','variant_label','alleles','MAF','strategies','PAM_targetable','excision_partners','WTB_quality','WTC_quality','WTD_quality','KOLF2_quality','source_gene']
    for symbol,rows in rows_by_gene.items():
        unique = {}
        for r in rows:
            key=tuple(r[:3]+[r[4]])
            if key not in unique: unique[key]=r[:]
            else:
                old=unique[key]
                old[6]=','.join(sorted(set(old[6].split(','))|set(r[6].split(','))))
                old[7]='Yes' if 'Yes' in [old[7],r[7]] else old[7]
                for cell,qualities in r[9].items():old[9].setdefault(cell,set()).update(qualities)
        rows=list(unique.values())
        strategies=defaultdict(set)
        for r in rows:
            for strategy in r[6].split(','): strategies[strategy.strip()].add(tuple(r[:3]+[r[4]]))
        cells=[]
        for cell in ['WTB','WTC','WTD','KOLF2']:
            matched=[r for r in rows if cell in r[9]]
            high=[r for r in matched if 'HIGH' in r[9][cell]]
            cells.append(dict(cell_line=cell,variants=len(matched),high_quality=len(high),low_only=len(matched)-len(high),pam_high_quality=sum(r[7]=='Yes' for r in high)))
        buf=io.StringIO();writer=csv.writer(buf);writer.writerow(headers)
        for r in rows_by_gene[symbol]: writer.writerow(r[:9]+['|'.join(sorted(r[9].get(cell,[]))) for cell in ['WTB','WTC','WTD','KOLF2']]+[r[10]])
        data=buf.getvalue().encode(); rel=f'downloads/{"sptlc-tracks" if focused else "genome-tracks"}/{symbol}.csv.gz'
        dest=OUT/rel;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(gzip.compress(data,mtime=0))
        if focused and 'genome_targets' in evidence[symbol]:evidence[symbol]['other_genome_targets']=evidence[symbol]['genome_targets']
        evidence[symbol]['genome_targets']=dict(source_id=source_id,dataset_title=label,source_symbols=sorted({r[10] for r in rows_by_gene[symbol]}),variants=len(rows),source_rows=len(rows_by_gene[symbol]),pam_variants=sum(r[7]=='Yes' for r in rows),strategies={k:len(v) for k,v in sorted(strategies.items())},cells=cells,export=dict(path=rel,rows=len(rows_by_gene[symbol]),sha256=hashlib.sha256(data).hexdigest()))
    save_json(OUT/('data/sptlc-track-audit.json' if focused else 'data/genome-track-audit.json'),dict(source_rows=source_rows,mapped_genes=len(rows_by_gene),unmapped_genes=dict(unmatched),renamed_genes=renamed,cells=cell_audit))
    return source_list

def build_clinvar(genes, evidence):
    path = CACHE / 'clinvar_gene_summary.txt'
    by_id = {str(g['entrez_id']):symbol for symbol,g in genes.items()}
    records = defaultdict(list)
    with path.open() as f:
        date = datetime.strptime(next(f).strip().split('dated ')[1], '%B %d, %Y').date().isoformat()
        for row in csv.DictReader(f,delimiter='\t'):
            symbol=by_id.get(row['GeneID'])
            if not symbol: continue
            records[symbol].append({k:(int(v) if v.isdigit() else None) if k not in ['#Symbol','GeneID','Gene_MIM_number'] else v for k,v in row.items()})
    conflicts=[]
    for symbol,rows in records.items():
        counts=[{k:v for k,v in r.items() if k!='Gene_MIM_number'} for r in rows]
        consistent=all(r==counts[0] for r in counts)
        evidence[symbol].update(clinvar=counts[0] if consistent else None,clinvar_records=rows,clinvar_date=date)
        if not consistent:conflicts.append(symbol)
    # Update browse/export counts from the same snapshot as the gene report.
    index_path=OUT/'data/index.json'; index=json.loads(index_path.read_text())
    symbol_col=index['fields'].index('symbol'); count_col=index['fields'].index('clinvar_total')
    for row in index['rows']:
        if row[symbol_col] in records:row[count_col]=(evidence[row[symbol_col]]['clinvar'] or {}).get('Total_alleles')
    save_json(index_path,index)
    save_json(OUT/'data/clinvar-current-audit.json',dict(date=date,mapped_genes=len(records),conflicts=conflicts))
    return [dict(id='clinvar_current',title='ClinVar gene summary · '+date,date=date,url='https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/gene_specific_summary.txt',path=str(path.relative_to(PROJECT)),sha256=digest(path),note='Distinct alleles in/overlapping the gene; not submissions. P/LP excludes overlapping structural variants. Identical duplicate OMIM gene rows are not summed. Historical dominant subset and variant downloads remain separate snapshots.')]

def expand_catalog(genes, track_bed):
    """Include approved source genes omitted by the earlier evidence-only catalog."""
    symbols=set()
    with track_bed.open() as f:
        for line in f:symbols.add(line.split('\t')[9].split()[-1])
    added={}
    with (PROJECT/'data_sources/hgnc_complete_set.txt').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            symbol=r['symbol']
            if symbol in genes or symbol not in symbols or r['status']!='Approved':continue
            added[symbol]=dict(symbol=symbol,name=r['name'],hgnc_id=r['hgnc_id'],ensembl_id=r['ensembl_gene_id'] or None,entrez_id=r['entrez_id'] or None,aliases=sorted(set(filter(None,(r['alias_symbol']+'|'+r['prev_symbol']).split('|')))),cytoband=r['location'],locus_type=r['locus_type'],omim_ids=list(filter(None,r['omim_id'].split('|'))),diseases=[],dosage=[],cmt_assertions=[])
    if not added:return
    genes.update(added)
    for key in {s[0] for s in added}:
        path=OUT/f'data/genes-{key}.json'; shard=json.loads(path.read_text())
        shard.update({s:g for s,g in added.items() if s[0]==key});save_json(path,shard)
    path=OUT/'data/index.json';idx=json.loads(path.read_text())
    for s,g in added.items():idx['rows'].append([s,g['name'],g['hgnc_id'],g['ensembl_id'],g['entrez_id'],g['aliases'],False,False,None,None,None,None,None])
    idx['rows'].sort(key=lambda r:r[0]);save_json(path,idx)
    manifest_path=OUT/'data/manifest.json';manifest=json.loads(manifest_path.read_text());manifest['counts']['genes']=len(genes);save_json(manifest_path,manifest)
    print('Added source genes to catalog:',len(added),flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fetch-mouse',action='store_true');parser.add_argument('--resume-mouse',action='store_true');parser.add_argument('--build',action='store_true');parser.add_argument('--tracks-only',action='store_true');parser.add_argument('--reuse-tracks',action='store_true')
    parser.add_argument('--track-bed',type=Path,default=PROJECT.parent/'github/dnd_browser_maintenance/build/original.bed')
    parser.add_argument('--cell-dir',type=Path,default=CACHE/'cells')
    parser.add_argument('--bigbed-tool',type=Path,default=Path.home()/'miniforge3/envs/dnd-browser-tools/bin/bigBedToBed')
    args=parser.parse_args();CACHE.mkdir(exist_ok=True,parents=True)
    if args.fetch_mouse: fetch_mouse(args.resume_mouse)
    if args.build or args.tracks_only:
        genes={}
        for path in (OUT/'data').glob('genes-*.json'): genes.update(json.loads(path.read_text()))
        expand_catalog(genes,args.track_bed)
        evidence=defaultdict(dict)
        if args.reuse_tracks:
            evidence.update(read_cache('genome-evidence'))
            track_sources=json.loads((CACHE/'tracks-sources.json').read_text())
            for src in track_sources:
                if digest(src['path'])!=src['sha256']:raise ValueError('Track source changed; rebuild tracks')
        else:
            track_sources=build_tracks(args,genes,evidence)
            args.track_bed=CACHE/'SPTLC_genes.bed'
            focused_sources=build_tracks(args,genes,evidence,focused=True)
            track_sources.extend(s for s in focused_sources if s['id']=='sptlc_tracks')
            (CACHE/'genome-evidence.json.gz').write_bytes(gzip.compress(json.dumps(evidence,separators=(',',':')).encode(),mtime=0))
            save_json(CACHE/'tracks-sources.json',track_sources)
        if args.tracks_only:return
        sources=build_mouse(genes,evidence)+track_sources+build_clinvar(genes,evidence)
        shards=defaultdict(dict)
        for symbol,row in evidence.items(): shards[symbol[0] if symbol[0].isalpha() else '_'][symbol]=row
        for key in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ_': save_json(OUT/f'data/evidence-{key}.json',shards[key])
        save_json(OUT/'data/evidence-sources.json',dict(sources=sources,coverage=dict(mouse_genes=sum('mouse_current' in r for r in evidence.values()),target_genes=sum('genome_targets' in r for r in evidence.values()))))
        print('Built evidence for',len(evidence),'genes',flush=True)

if __name__=='__main__':main()
