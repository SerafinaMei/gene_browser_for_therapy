"""Import complete per-gene RUN_MULTIALLELIC records; reject incomplete exports."""
import copy,csv,gzip,hashlib,json,shutil
from pathlib import Path
from collections import defaultdict

ROOT=Path(__file__).resolve().parents[1]
CACHE=ROOT/'data_sources/comprehensive'
RUN=ROOT.parent/'Grace_pipeline/high_coverage_result/RUN_MULTIALLELIC'
FIELDS={'num_indel_vars':'exon disruption','num_crisproff_vars':'epigenetic silencing','num_excision_vars':'excision','num_base_edit_snps':'splice site disruption'}

def digest(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def build_pipeline(args,genes,evidence,build_tracks):
    summary_path=RUN/'summary_files/dnd_summary.csv'
    source_rows=list(csv.DictReader(summary_path.open()))
    summaries={}
    for r in source_rows:
        symbol=r['hgnc_symbol']
        if symbol in summaries:
            assert all(r[k]==summaries[symbol][k] for k in FIELDS),f'Conflicting summary: {symbol}'
        summaries[symbol]=r
    expected={s:{strategy:int(float(r.get(field) or 0)) for field,strategy in FIELDS.items()} for s,r in summaries.items()}
    sources=[];observed={};records=0
    bed=CACHE/'run_multiallelic.bed'
    with bed.open('w') as out:
        for p in sorted((RUN/'summary_files/browser_tracks/per_gene_files').rglob('*.bed')):
            symbol=p.parent.name
            assert symbol in expected,f'Unexpected source gene: {symbol}'
            keys=defaultdict(set)
            with p.open() as f:
                for line in f:
                    c=line.rstrip('\n').split('\t')
                    assert len(c)==25,f'Incomplete BED row in {p}'
                    assert c[9].split()[-1]==symbol
                    key=(c[0],c[1],c[2],c[10])
                    for strategy in c[19].split(','):keys[strategy.strip()].add(key)
                    out.write(line);records+=1
            observed[symbol]={strategy:len(keys[strategy]) for strategy in FIELDS.values()}
            assert observed[symbol]==expected[symbol],f'Incomplete or inconsistent BED: {symbol}'
            sources.append(dict(path=str(p),sha256=digest(p),bytes=p.stat().st_size))
    positive={s for s,v in expected.items() if any(v.values())}
    assert set(observed)==positive,'Candidate genes missing from per-gene BED export'
    shutil.copy2(summary_path,CACHE/'run_multiallelic_summary.csv')
    audit=dict(status='PASS',summary_genes=len(summaries),genes_with_candidates=len(positive),zero_candidate_genes=sorted(set(summaries)-positive),records=records,summary_sha256=digest(summary_path),bed_sha256=digest(bed),files=sources,note='All four strategy counts agree for every gene. Empty/truncated combined BEDs are not used.')
    (ROOT/'dist/data/run-multiallelic-validation.json').write_text(json.dumps(audit,separators=(',',':')))
    a=copy.copy(args);a.track_bed=bed
    imported=build_tracks(a,genes,evidence,dataset_id='multiallelic')
    by_hgnc={g['hgnc_id']:s for s,g in genes.items()}
    for original,row in summaries.items():
        symbol=by_hgnc[row['GENE ID (HGNC)']]
        if original not in positive:
            evidence[symbol]['pipeline_targets']=dict(source_id='run_multiallelic',dataset_title='RUN_MULTIALLELIC · complete per-gene export',source_symbols=[original],variants=0,source_rows=0,pam_variants=0,strategies=expected[original],cells=[dict(cell_line=c,variants=0,high_quality=0,low_only=0,pam_high_quality=0) for c in ['WTB','WTC','WTD','KOLF2']],export=None)
        target=evidence[symbol]['pipeline_targets']
        for strategy,n in expected[original].items():assert target['strategies'].get(strategy,0)==n
        target['strategies']=expected[original]
        target['excision_pairs']=int(float(row.get('num_excision_pairs') or 0))
    return [r for r in imported if r['id']=='run_multiallelic']+[dict(id='run_multiallelic_summary',title='RUN_MULTIALLELIC · gene summary and validation',date=None,path=str(summary_path),sha256=digest(summary_path),url='https://github.com/gramey02/dnd_project',note='593 assessed genes; 568 with candidates, 25 with zero candidates. All per-strategy counts independently reconciled to complete per-gene BEDs. Source run date unspecified; folder name is not treated as a release date.')]
