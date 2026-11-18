"""Import validated saved greedy curves and render gene-specific population plots.

First import: --run-root PATH --master-file PATH
Subsequent builds use data_sources/population-coverage.json.gz.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json,os
from collections import defaultdict
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/dnd-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

PROJECT=Path(__file__).resolve().parents[1]
CACHE=PROJECT/'data_sources/population-coverage.json.gz'
OUT=PROJECT/'dist'
LABELS={'indels':'Exon disruption','CRISPRoff':'Epigenetic silencing','base_edits':'Splice-site disruption','excision':'Excision'}
SUFFIX={'indels':'exon_disr','CRISPRoff':'epi_sil','base_edits':'ss_disr','excision':'excision'}
PAPER='https://pmc.ncbi.nlm.nih.gov/articles/PMC13042096/#F5'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def save_json(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,separators=(',',':'),allow_nan=False))

def validate_curve(points,total,possible,strategy):
    previous=0; guides=0
    for p in points:
        assert p['added']>0 and p['covered']==previous+p['added']
        assert p['guides']-guides in ([1,2] if strategy=='excision' else [1])
        assert 0<=p['covered']<=possible<=total
        previous=p['covered'];guides=p['guides']
    assert points

def vcf_cohort(path):
    with path.open() as f:
        for line in f:
            if line.startswith('#CHROM'):
                names=line.rstrip().split('\t')[9:]
                return len(names),hashlib.sha256('\n'.join(names).encode()).hexdigest()
    raise ValueError(f'No sample header: {path}')

def import_run(root,master_path):
    master={}
    for r in csv.DictReader(master_path.open()):
        symbol=r['hgnc_symbol']
        if symbol in master:
            assert master[symbol]['dominant_mutation_count']==r['dominant_mutation_count']
        master[symbol]=r
    result={'genes':{},'excluded':{},'sources':{'run':str(root),'master':str(master_path),'master_sha256':sha(master_path)},'cohort':'1000 Genomes · GRCh38','baseline_date':'2026-08-25'}
    cohort=None
    for folder in ['non_excision_guides','excision_guides']:
        parent=root/'summary_files/cross_strat_gRNAs'/folder
        for path in sorted((parent/'results').glob('*_gRNAs.csv')):
            symbol,strategy=path.stem.removesuffix('_gRNAs').split('_',1)
            assert strategy in LABELS
            m=master[symbol]
            if m['chrom'] not in ['chr'+str(i) for i in range(1,23)]:
                result['excluded'][symbol]='Sex-chromosome coverage requires a ploidy-specific analysis.'
                continue
            assert (parent/'logs'/f'{symbol}.DONE_{strategy}').exists(),path
            vcfs=[root/s/'excavate/Guide_filtered_vcfs'/f'{symbol}_guide_filtered.vcf' for s in (['donor_base_edits','acceptor_base_edits'] if strategy=='base_edits' else [strategy])]
            existing=[p for p in vcfs if p.exists()]
            assert existing,path
            for vcf in existing:
                current=vcf_cohort(vcf)
                if cohort is None:cohort=current
                assert cohort==current,(vcf,cohort,current)
            points=[]
            for r in csv.DictReader(path.open()):
                points.append(dict(guides=int(r['num_guides']),added=int(r['haplotypes_added']),covered=int(r['cumulative_haplotypes_targeted']),allele1=r['snp1_allele'],allele2=r['snp2_allele']))
            possible=2*int(m['num_hets_preCas_'+SUFFIX[strategy]])
            validate_curve(points,cohort[0]*2,possible,strategy)
            assert points[-1]['covered']==float(m['num_haps_targeted_all_guides_'+strategy.lower()])
            g=result['genes'].setdefault(symbol,dict(hgnc_id=m['GENE ID (HGNC)'],chromosome=m['chrom'],strategies={},dominant_mutations=None))
            n=m['dominant_mutation_count']
            if n and n.lower() not in ('nan','none') and float(n)>0:g['dominant_mutations']=int(float(n))
            g['strategies'][strategy]=dict(points=points,possible=possible,source=str(path.relative_to(root)),sha256=sha(path))
    result['samples'],result['sample_order_sha256']=cohort
    result['haplotypes']=cohort[0]*2
    CACHE.write_bytes(gzip.compress(json.dumps(result,separators=(',',':')).encode(),mtime=0))
    print(f"Imported {len(result['genes'])} genes; {sum(len(g['strategies']) for g in result['genes'].values())} curves; {result['samples']} matched samples; {len(result['excluded'])} sex-chromosome genes withheld.",flush=True)

def draw(symbol,strategy,g,series,total,path):
    points=series['points'];x=[0]+[p['guides'] for p in points];y=[0]+[p['covered']/total*100 for p in points]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none','svg.hashsalt':'dnd-population','axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(9,4.6),layout='constrained')
    ax.plot(x,y,color='#326b8b',lw=2.2,marker='o',markersize=3,label='Mutation-agnostic targeting')
    ax.axhline(series['possible']/total*100,color='#737b82',ls='--',lw=1.4,label='Heterozygous coverage before PAM filtering')
    if g['dominant_mutations'] and strategy!='excision':
        ax.plot(x,[min(k/g['dominant_mutations'],1)*100 for k in x],color='#b34e4e',ls='--',lw=1.7,label='Mutation-specific model (equal-frequency alleles)')
    ax.set(xlim=(0,max(x)),ylim=(0,102),xlabel='Number of prioritized gRNAs' if strategy=='excision' else 'Number of prioritized therapies (one gRNA each)',ylabel='Cohort haplotypes covered (%)',title=f'{symbol} · {LABELS[strategy]}')
    ax.xaxis.set_major_locator(MaxNLocator(nbins=10,integer=True));ax.set_yticks([0,20,40,60,80,100]);ax.grid(axis='y',color='#e4e8eb',lw=.7);ax.set_axisbelow(True)
    ax.legend(loc='upper left',bbox_to_anchor=(0,1.23),frameon=False,fontsize=8,ncols=1)
    fig.savefig(path,format='svg',metadata={'Date':None,'Description':f'{total} haplotypes; saved greedy prioritization results.'})
    path.write_text(path.read_text().replace("font-family: 'DejaVu Sans'","font-family: 'DejaVu Sans', Arial, sans-serif"))
    path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines())+"\n")
    if symbol=='MYH7' and strategy=='indels':fig.savefig('/tmp/MYH7-population-coverage.png',dpi=160)
    plt.close(fig)

def build():
    with gzip.open(CACHE,'rt') as f:data=json.load(f)
    shards={p.stem.removeprefix('genes-'):{} for p in (OUT/'data').glob('genes-*.json')}
    genes={}
    for p in (OUT/'data').glob('genes-*.json'):genes.update(json.loads(p.read_text()))
    by_id={g['hgnc_id']:s for s,g in genes.items()}
    plotdir=OUT/'plots/population';plotdir.mkdir(parents=True,exist_ok=True)
    downloaddir=OUT/'downloads/population';downloaddir.mkdir(parents=True,exist_ok=True)
    for i,(original,g) in enumerate(sorted(data['genes'].items())):
        symbol=by_id[g['hgnc_id']]
        record=dict(samples=data['samples'],haplotypes=data['haplotypes'],cohort=data['cohort'],dominant_mutations=g['dominant_mutations'],baseline_date=data['baseline_date'],strategies={})
        for strategy,series in g['strategies'].items():
            validate_curve(series['points'],data['haplotypes'],series['possible'],strategy)
            name=f'{symbol}-{strategy}'
            plot=plotdir/(name+'.svg')
            draw(symbol,strategy,g,series,data['haplotypes'],plot)
            csvpath=downloaddir/(name+'.csv.gz')
            buf=io.StringIO();writer=csv.writer(buf)
            writer.writerow(['prioritized_gRNAs','new_haplotypes','cumulative_haplotypes','cohort_haplotypes','coverage_percent','prePAM_possible_haplotypes','selected_allele_1','selected_allele_2','mutation_specific_model_percent'])
            for p in [dict(guides=0,added=0,covered=0,allele1='',allele2='')]+series['points']:
                baseline=min(p['guides']/g['dominant_mutations'],1)*100 if g['dominant_mutations'] and strategy!='excision' else ''
                writer.writerow([p['guides'],p['added'],p['covered'],data['haplotypes'],p['covered']/data['haplotypes']*100,series['possible'],p['allele1'],p['allele2'],baseline])
            csvpath.write_bytes(gzip.compress(buf.getvalue().encode(),mtime=0))
            record['strategies'][strategy]=dict(label=LABELS[strategy],plot=str(plot.relative_to(OUT)),export=dict(path=str(csvpath.relative_to(OUT))),possible=series['possible'],points=series['points'])
        shards[symbol[0]][symbol]={'population_coverage':record}
        if i%100==0:print('Rendered',i+1,flush=True)
    for symbol,reason in data['excluded'].items():
        if symbol in genes:shards[symbol[0]][symbol]={'population_coverage':{'unavailable':reason}}
    for key,rows in shards.items():save_json(OUT/'data'/f'population-{key}.json',rows)
    manifest=json.loads((OUT/'data/manifest.json').read_text())
    manifest['sources']=[s for s in manifest['sources'] if s['id']!='population_coverage']+[dict(id='population_coverage',title='Population haplotype coverage',url=PAPER,date=data['baseline_date'],path=str(CACHE.relative_to(PROJECT)),sha256=sha(CACHE),note=f"Saved GRCh38 greedy results for {len(data['genes'])} autosomal genes, {data['samples']:,} samples and {data['haplotypes']:,} haplotypes. Completion markers, cohort sample order, cumulative counts and summary totals checked. The mutation-specific model assumes equally frequent dominant disease alleles; its historical counts are separate from current ClinVar totals. Excision is shown per unique gRNA, including reuse, not per independent therapy. Allele suffix 0 denotes reference and 1 alternate. No population curve is inferred from cell-line counts.")]
    save_json(OUT/'data/manifest.json',manifest)
    print('Population plots complete.',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run-root',type=Path);parser.add_argument('--master-file',type=Path);args=parser.parse_args()
    if args.run_root:
        if not args.master_file:parser.error('--master-file is required with --run-root')
        import_run(args.run_root,args.master_file)
    build()
