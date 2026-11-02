"""Import the completed D&D/HI indel and splice screen as a separate dataset."""
import ast
import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

RUN_RELATIVE = 'Grace_pipeline/pipeline_recheck_20261101'


def read(path):
    with (gzip.open if path.suffix == '.gz' else open)(path, 'rt', encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def positions(text):
    # The upstream CSV serializes numbers as np.float64(...). Parse as data,
    # without importing or executing the serialized expression.
    values = ast.literal_eval(re.sub(r'np\.float64\(([-+0-9.eE]+)\)', r'\1', text))
    if not isinstance(values, list) or not all(isinstance(v, (int, float)) and int(v) == v for v in values):
        raise ValueError('Unexpected NMD position list')
    return {int(v) for v in values}


def build_screen(root, out, genes, details, source, write_csv):
    run = root / RUN_RELATIVE
    validation = json.loads((run / 'results/validation.json').read_text())
    if validation.get('status') != 'PASS':
        raise ValueError('Expanded targeting screen has not passed source validation')
    annotations = read(run / 'results/gene_annotations.csv')
    by_symbol = {s.upper(): s for s in genes}
    by_id = {g['hgnc_id']: s for s, g in genes.items()}
    mappings = {}
    for row in annotations:
        symbol = next((by_id[row[k]] for k in ('dosage_HGNC_ID', 'clingen_HGNC_ID') if row.get(k) in by_id), None) or by_symbol.get(row['gene'])
        if not symbol:
            raise ValueError(f'Unmapped screen gene: {row["gene"]}')
        mappings[row['gene']] = symbol
    if len(set(mappings.values())) != len(mappings):
        raise ValueError('Duplicate source genes map to the same HGNC gene')
    models = {r['hgnc_symbol'] for r in read(run / 'inputs/pipeline_exons.csv')}
    excluded_nmd = {r['gene'] for r in read(run / 'inputs/nmd_model_exclusions.csv')}
    nmd = {r['gene']: r for r in read(run / 'stage_a/indels/NMD/NMD_induction_summary.csv')}
    nmd_positions = {r['gene']: {kind: positions(r[f'vars_consistently_{kind}_NMD']) for kind in ('inducing', 'escaping')}
                     for r in read(run / 'stage_a/indels/NMD/NMD_induction_var_info.csv')}
    variants, guides = defaultdict(list), defaultdict(list)
    variant_path = run / 'results/snps_annotated_prePAM.csv'
    guide_path = run / 'results/all_guide_associations.csv.gz'
    for row in read(variant_path):
        variants[row['gene']].append(row)
    for row in read(guide_path):
        guides[row['gene']].append(row)

    for original, symbol in mappings.items():
        rows, grnas = variants[original], guides[original]
        assessed = original in models
        exon = [r for r in rows if r['editor'] == 'Indel']
        if original in nmd:
            if len({r['variant_id'] for r in exon}) != int(nmd[original]['num_vars_considered']):
                raise ValueError(f'Exonic candidate count mismatch: {original}')
            for kind in ('inducing', 'escaping'):
                if len(nmd_positions[original][kind]) != int(nmd[original][f'num_vars_{kind}_NMD']):
                    raise ValueError(f'NMD count mismatch: {original}')
        elif exon and original not in excluded_nmd:
            raise ValueError(f'Missing NMD results for exonic candidates: {original}')

        def measurement(key, label, candidates, available=True, applicable=True):
            available = available and assessed
            ids = {r['variant_id'] for r in candidates}
            pam = any(r['variant_id'] in ids and r['passes_10bp_SNP_PAM'] == 'True' for r in grnas)
            return dict(key=key, label=label, count=len(ids) if available else None,
                        targetable=bool(ids) if available and applicable else None,
                        pam=pam if available and applicable else None,
                        status='assessed' if available else 'not_assessed')

        inducing = nmd_positions.get(original, {}).get('inducing', set())
        escaping = nmd_positions.get(original, {}).get('escaping', set())
        measures = [measurement('indel_pre', 'Exonic indel, before NMD assessment', exon),
                    measurement('indel_nmd', 'Exonic indel, predicted NMD', [r for r in exon if int(r['pos']) in inducing], original not in excluded_nmd),
                    measurement('indel_escape', 'Exonic indel, predicted NMD escape', [r for r in exon if int(r['pos']) in escaping], original not in excluded_nmd, False),
                    measurement('splice', 'Splice base editing (1-23 bp)', [r for r in rows if r['window_set'] == '1-23'])]
        for key, label in [('excision','Exon / gene excision'),('crisproff','CRISPRoff'),('excision_pairs','Exon / gene excision pairs')]:
            measures.append(dict(key=key,label=label,count=None,targetable=None,pam=None,status='not_run'))
        exports = {}
        for kind, data in [('variants',rows),('guides',grnas)]:
            if not data:
                continue
            relative = f'downloads/expanded/{symbol}-{kind}.csv.gz'
            headers = list(data[0])
            write_csv(relative, headers, ([r[k] for k in headers] for r in data))
            exports[kind] = dict(path=relative,rows=len(data),unique=len({r['variant_id' if kind=='variants' else 'guide_id'] for r in data}))
        details[symbol]['strategy_screen'] = dict(assessed=assessed,rows=measures,exports=exports,
            source_gene=original,source_id='expanded_screen',
            exclusion=None if assessed else 'No chromosome Y allele-frequency input in this run.')

        # Cell-line counts stay separate from population-level strategy counts.
        # Only genes without the earlier cell-line dataset use these as fallback.
        if assessed:
            scopes=[]
            for window in ('1-23','1-30'):
                subset=[r for r in rows if r['window_set'] in (window,'exonic')]
                matching=[r for r in grnas if r['window_set'] in (window,'exonic')]
                line_rows=[]
                for line,flag in [('All candidates',None),('WTC','het_WTC'),('WTD','het_WTD'),('p47','het_p47_locus_only')]:
                    eligible=subset if flag is None else [r for r in subset if r[flag]=='True']
                    ids={r['variant_id'] for r in eligible}
                    guide_subset=[r for r in matching if r['variant_id'] in ids]
                    selected=[r for r in guide_subset if r['primary_qualified']=='True']
                    line_rows.append(dict(cell_line=line,variants=len(ids),
                        ngg_variants=len({r['variant_id'] for r in guide_subset if r['passes_10bp_SNP_PAM']=='True'}),
                        selected_variants=len({r['variant_id'] for r in selected}),guides=len({r['guide_id'] for r in selected}),
                        editors=sorted({r['editor'] for r in selected})))
                scopes.append(dict(window=window,rows=line_rows))
            details[symbol]['expanded_targeting']=dict(scopes=scopes,exports=exports,dataset='expanded')

    audited=[]
    for relative in ['Grace_pipeline/high_coverage_result/RUN_MULTIALLELIC/summary_files/dnd_summary.csv',
                     'Grace_pipeline/261011_RUN_final/summary_files/dnd_summary.csv',
                     'Grace_pipeline/RUN_SMEI_STAGEA_FILTERED_0608/indels/ubiq_region_CommonVars/CommonVars_ALL_summary.txt']:
        p=root/relative
        with p.open() as f:
            reader=csv.DictReader(f,delimiter='\t' if p.suffix=='.txt' else ',')
            records=list(reader)
        field='hgnc_symbol' if 'hgnc_symbol' in records[0] else 'gene'
        audited.append(dict(path=relative,genes=len({r[field] for r in records}),rows=len(records),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    audit=dict(previous_strategy_genes=sum('dnd' in g for g in genes.values()),
        examined_run_summaries=audited,expanded_input_genes=len(mappings),expanded_assessed_genes=len(models),
        extra_input_genes=sum('dnd' not in genes[s] for s in mappings.values()),
        extra_assessed_genes=sum('dnd' not in genes[mappings[s]] for s in models),
        missing_chromosome_inputs=sorted(set(mappings)-models),nmd_model_exclusions=sorted(excluded_nmd),
        expanded_strategies=['exonic indel','NMD','splice base editing'],
        expanded_strategies_not_run=['excision','CRISPRoff'],
        all_gene_run_found=False,remote_wynton_path_mounted=(Path('/wynton/group/capra/projects/dnd_project_results/RUN_MULTIALLELIC/RUN_MULTIALLELIC')).exists(),
        source_validation=validation)
    (out/'data/targeting-coverage.json').write_text(json.dumps(audit,indent=2))
    inputs=[('expanded_screen',run/'results/gene_annotations.csv','Expanded indel and splice screen','D&D, ClinGen HI and two additional control genes. Results are separate from the D&D reference workbook; CRISPRoff and excision were not run.'),
            ('expanded_variants',variant_path,'Expanded-screen SNP records','Common SNPs in shared exonic and splice regions. Splice distances are 1-23 bp and 1-30 bp; AF is 0.1-0.9. These are candidate counts before cell-line expression and genotype selection.'),
            ('expanded_guides',guide_path,'Expanded-screen gRNA records','NGG guides with source sequence, geometry, NMD and selection flags. The primary screen uses the 10 bp SNP-to-PAM rule; broader guide records retain their filter flags.'),
            ('expanded_models',run/'inputs/pipeline_exons.csv','Expanded-screen transcript models','Defines assessed genes and chromosomes; absence from variant records alone is never used to infer a tested zero.'),
            ('expanded_nmd',run/'stage_a/indels/NMD/NMD_induction_summary.csv','Expanded-screen NMD results','Predictions consistent across retained transcripts. Mixed predictions need not be classified as inducing or escaping NMD.'),
            ('expanded_nmd_positions',run/'stage_a/indels/NMD/NMD_induction_var_info.csv','Expanded-screen NMD variant positions','Positions with consistent NMD induction or escape across retained transcripts.'),
            ('expanded_nmd_exclusions',run/'inputs/nmd_model_exclusions.csv','Expanded-screen NMD exclusions','NMD cannot be evaluated when a retained transcript lacks coding coordinates.')]
    return [source(path,key,title,note,date='2026-11-01') for key,path,title,note in inputs]
