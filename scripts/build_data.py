"""Build a static, auditable gene browser from source snapshots; never edit inputs."""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import openpyxl

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT.parent
OUT = PROJECT / 'dist/data'
SOURCES = []
AUDIT = {'unmatched': defaultdict(list), 'conflicts': []}
GENES = {}
BY_HGNC, BY_ENSEMBL, BY_ENTREZ = {}, {}, {}
ALIASES = defaultdict(set)


def clean(v):
    if isinstance(v, float) and not math.isfinite(v):
        return None
    if isinstance(v, str):
        v = v.strip()
        if v in ('', '-', '‐', '—', "'-", 'NoneReported', 'NA', 'NaN', 'nan'):
            return None
    return v


def number(v):
    v = clean(v)
    if v is None:
        return None
    try:
        n = float(v)
        return int(n) if n.is_integer() else n
    except (ValueError, TypeError):
        return None


def source(key, path, title, note, url=None, date=None):
    p = ROOT / path
    SOURCES.append(dict(id=key, path=path, title=title, note=note, url=url, date=date,
                        sha256=hashlib.sha256(p.read_bytes()).hexdigest(), bytes=p.stat().st_size))
    return p


def table(p, sep=',', skip=0):
    op = gzip.open if str(p).endswith('.gz') else open
    with op(p, 'rt', encoding='utf-8-sig', newline='') as f:
        for _ in range(skip):
            next(f)
        for r in csv.DictReader(f, delimiter=sep):
            yield {k: clean(v) for k, v in r.items() if k is not None}


def sheet(p, name):
    w = openpyxl.load_workbook(p, read_only=True, data_only=True)
    it = w[name].values
    headers = next(it)
    rows = [{k: clean(v) for k, v in zip(headers, row) if k} for row in it]
    w.close()
    return rows


def gene(symbol=None, hgnc=None, ensembl=None, entrez=None):
    s = BY_HGNC.get(hgnc) or BY_ENSEMBL.get(ensembl) or BY_ENTREZ.get(str(entrez))
    if not s and symbol in GENES:
        s = symbol
    if not s and len(ALIASES.get(symbol, [])) == 1:
        s = next(iter(ALIASES[symbol]))
    return GENES.get(s)


def unique_record(g, key, row, src):
    if key in g and g[key] != row:
        AUDIT['conflicts'].append({'gene': g['symbol'], 'source': src, 'field': key})
        raise ValueError(f'Conflicting gene-level {key}: {g["symbol"]}')
    g[key] = row


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    p = source('hgnc', 'DnD_gene_browser/data_sources/hgnc_complete_set.txt', 'HGNC identifiers',
               'Approved symbols, aliases, stable IDs, names and OMIM gene identifiers. Retrieved 2026-11-01.',
               'https://www.genenames.org/download/archive/', '2026-11-01')
    for r in table(p, '\t'):
        if r['status'] != 'Approved':
            continue
        s = r['symbol']
        aliases = sorted(set((r['alias_symbol'] or '').split('|') + (r['prev_symbol'] or '').split('|')) - {''})
        GENES[s] = dict(symbol=s, name=r['name'], hgnc_id=r['hgnc_id'], ensembl_id=r['ensembl_gene_id'],
                        entrez_id=r['entrez_id'], aliases=aliases, cytoband=r['location'], locus_type=r['locus_type'],
                        omim_ids=(r['omim_id'] or '').split('|') if r['omim_id'] else [], diseases=[], dosage=[])
        BY_HGNC[r['hgnc_id']] = s
        if r['ensembl_gene_id']:
            BY_ENSEMBL[r['ensembl_gene_id']] = s
        if r['entrez_id']:
            BY_ENTREZ[r['entrez_id']] = s
        for a in aliases:
            ALIASES[a].add(s)

    p = source('validity', 'DnD_project/data/Clingen-Gene-Disease-Summary-2025-11-13.csv', 'ClinGen gene–disease validity',
               'All assertion rows retained with disease, MOI, classification, expert panel and evaluation date.',
               'https://search.clinicalgenome.org/kb/gene-validity', '2025-11-13')
    for r in table(p, skip=4):
        if not (r.get('GENE ID (HGNC)') or '').startswith('HGNC:'):
            continue
        g = gene(hgnc=r['GENE ID (HGNC)'])
        if not g:
            AUDIT['unmatched']['validity'].append(r['GENE SYMBOL']); continue
        g['diseases'].append(r)

    p = source('dosage', 'DnD_project/data/Clingen-Dosage-Sensitive-HI-DETAIL.csv', 'ClinGen dosage sensitivity',
               'Gene records only. Regions remain separate. %HI is the DECIPHER HI percentile rank (lower indicates greater predicted HI), not a probability.',
               'https://search.clinicalgenome.org/kb/gene-dosage')
    score_map = {'NoEvidence': 0, 'LittleEvidence': 1, 'SomeEvidence': 2, 'SufficientEvidence': 3, 'AutosomalRecessive': 30, 'SensitivityUnlikely': 40}
    regions = []
    for r in table(p):
        if r['Type'] != 'Gene':
            if 'PMP22' in (r['Gene Symbol /Region Name'] or ''):
                regions.append(r)
            continue
        g = gene(hgnc=r['HGNC/Dosage ID'], symbol=r['Gene Symbol /Region Name'])
        if not g:
            AUDIT['unmatched']['dosage'].append(r['Gene Symbol /Region Name']); continue
        g['dosage'].append(r)
        g['hi_score'] = score_map.get(r['HI Score'])
        g['hi_evidence'] = r['HI Score']
        g['hi_percentile'] = number(r['%HI'])
        g['coordinates'] = r['GRCh38']

    p = source('shet', 'DnD_project/data/Spense_S-het_supp_table.xlsx', 'GeneBayes s_het · Spence et al.',
               'Supplementary Table 1: post_mean and 95% posterior credible interval. Evolutionary selection against heterozygous loss of function; not an HI probability.',
               'https://pmc.ncbi.nlm.nih.gov/articles/PMC10312940/')
    for r in sheet(p, 'Supplementary Table 1'):
        g = gene(hgnc=r['hgnc'], ensembl=r['ensg'])
        if g:
            unique_record(g, 'shet', r, 'shet')
        else:
            AUDIT['unmatched']['shet'].append(r['ensg'])

    p = source('collins', 'DnD_project/data/Cell_paper/Cell_Haplo_pred.xlsx', 'Collins et al. · Cell 2022',
               'Table S7: pHaplo is the predicted probability of haploinsufficiency (0–1); pTriplo is probability of triplosensitivity. Higher pHaplo means greater predicted HI.',
               'https://pmc.ncbi.nlm.nih.gov/articles/PMC9742861/', '2022')
    for r in sheet(p, 'Table S7'):
        g = gene(symbol=r['Gene'])
        if g:
            g.setdefault('collins_records', []).append(r)
        else:
            AUDIT['unmatched']['collins'].append(r['Gene'])
    for g in GENES.values():
        if len(g.get('collins_records', [])) == 1:
            g['collins'] = g['collins_records'][0]
        elif len(g.get('collins_records', [])) > 1:
            AUDIT['conflicts'].append({'gene':g['symbol'], 'source':'collins', 'field':'pHaplo',
                                      'resolution':'All historical records retained; headline prediction withheld.'})

    p = source('dnd', 'DnD_project/data/dnd_593_2025_12_23_w_dplp_muts.xlsx', 'D&D targeting and variant snapshot',
               'Gene-level targeting counts and local pathogenic/likely pathogenic (P/LP) counts. Duplicate disease rows are not summed. Dominant P/LP counts retain source annotation; variant-level dominant assignments are not supplied.',
               'https://github.com/gramey02/dnd_project', '2025-12-23 (filename)')
    ignored = {'DISEASE LABEL', 'DISEASE ID (MONDO)', 'MOI', 'CLASSIFICATION', 'Gene Group'}
    for r in sheet(p, 'dnd_593_2025_12_23_w_dplp_muts'):
        g = gene(hgnc=r['GENE ID (HGNC)'], symbol=r['hgnc_symbol'])
        if g:
            unique_record(g, 'dnd', {k:v for k,v in r.items() if k not in ignored}, 'dnd')
    dictionary = {r['Column Name']: r['Description'] for r in sheet(p, 'Column_Descriptions')}

    p = source('dominant_legacy', 'DnD_project/data/ClinVar_Omim_dom_annot_reduced_cols.csv', 'Legacy ClinVar–OMIM dominant counts',
               'Original dominant_mutation_count retained separately; source snapshot date and per-variant inclusion rules are unavailable. Do not combine with the D&D workbook dominant P/LP count.')
    for r in table(p):
        g = gene(symbol=r['gene'])
        if g:
            g['legacy_dominant_count'] = number(r['dominant_mutation_count'])

    p = source('clinvar', 'DnD_gene_browser/data_sources/clinvar_gene_summary.txt', 'ClinVar gene summary',
               'Total_alleles counts distinct alleles in/overlapping a gene, not submissions. P/LP summary excludes overlapping structural variants. File statistics can lag the report date. Not directly comparable to the older D&D dominant subset.',
               'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/gene_specific_summary.txt', '2026-10-24')
    for r in table(p, '\t', skip=1):
        g = gene(entrez=r['GeneID'])
        if g:
            g.setdefault('clinvar_records', []).append({k:number(v) if k not in ('#Symbol', 'GeneID', 'Gene_MIM_number') else v for k,v in r.items()})
    for g in GENES.values():
        records = g.get('clinvar_records', [])
        if records:
            counts = [{k:v for k,v in r.items() if k != 'Gene_MIM_number'} for r in records]
            if all(r == counts[0] for r in counts):
                g['clinvar'] = counts[0]
            else:
                AUDIT['conflicts'].append({'gene':g['symbol'],'source':'clinvar','resolution':'All records retained; conflicting total withheld.'})

    # Homology is joined by MGI identifiers, never by capitalizing mouse symbols.
    p = source('homology', 'DnD_project/MousePhenotype/HOM_MouseHumanSequence.rpt', 'MGI human–mouse homology',
               'DB Class Key links human HGNC IDs to mouse MGI identifiers. One-to-many homologs are retained.',
               'https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt')
    classes = defaultdict(list)
    for r in table(p, '\t'):
        classes[r['DB Class Key']].append(r)
    mouse_to_genes = defaultdict(set)
    for rows in classes.values():
        humans = [gene(hgnc=r['HGNC ID']) for r in rows if r['NCBI Taxon ID'] == '9606']
        for r in rows:
            if r['NCBI Taxon ID'] != '10090':
                continue
            for g in filter(None, humans):
                mouse_to_genes[r['Mouse MGI ID']].add(g['symbol'])

    mouse = defaultdict(lambda: {'records': {}, 'hom': set(), 'sources': set(), 'mouse_symbols': set(), 'covered': False})
    for cohort in ('DND593', 'HI417', 'Dominant855', 'AR754'):
        name = f'{cohort}_heterozygous_statistical_results_all_genes.txt.gz'
        p = source('mouse_' + cohort, 'DnD_project/MousePhenotype/' + name, f'IMPC heterozygous results · {cohort}',
                   'Successful statistical results, significant and non-significant. Rows deduplicated by doc_id across cohorts; parameters by parameter_name. Source significant flag is preserved.')
        mapname = {'DND593':'dnd593', 'HI417':'HI417', 'Dominant855':'dominant855', 'AR754':'AR754'}[cohort]
        mapping = ROOT / f'DnD_project/MousePhenotype/human_mouse_homologs_{mapname}.txt'
        source('mouse_map_' + cohort, str(mapping.relative_to(ROOT)), f'Mouse cohort coverage · {cohort}', 'Defines source-cohort coverage; an absent gene is not assigned zero mouse results.')
        for r in table(mapping, '\t'):
            g = gene(hgnc=r.get('GENE ID (HGNC)'), symbol=r.get('hgnc_symbol'))
            if g:
                mouse[g['symbol']]['covered'] = True
        for r in table(p, '\t'):
            if r['zygosity'] != 'heterozygote' or r['status'] != 'Successful':
                continue
            for s in mouse_to_genes[r['marker_accession_id']]:
                m = mouse[s]
                m['covered'] = True
                m['records'][r['doc_id'] or json.dumps(r, sort_keys=True)] = r
                m['sources'].add('mouse_' + cohort)
                m['mouse_symbols'].add(r['marker_symbol'])
        homfile = f'{cohort}_MGI_IMPC_phenotypes_merged.txt' if cohort != 'AR754' else 'MGI_IMPC_phenotypes_AR754_het_hom.txt'
        p = source('assertions_' + cohort, 'DnD_project/MousePhenotype/' + homfile, f'MGI/IMPC assertions · {cohort}',
                   'Homozygote assertion presence only. Absence from these exports does not prove a lack of homozygous testing.')
        for r in table(p, '\t'):
            if r['zygosity'] != 'homozygote':
                continue
            for s in mouse_to_genes[r['marker_accession_id']]:
                mouse[s]['hom'].add(r['marker_accession_id'])
                mouse[s]['sources'].add('assertions_' + cohort)

    for s,m in mouse.items():
        if not m['covered'] and not m['hom']:
            continue
        rows = list(m['records'].values())
        sig = [r for r in rows if r['significant'] == 'True']
        result = dict(hgnc_symbol=s, het_rows=len(rows), het_significant_rows=len(sig),
                      het_nonsignificant_rows=sum(r['significant']=='False' for r in rows),
                      het_significance_unknown_rows=sum(r['significant'] not in ('True','False') for r in rows),
                      has_hom_record=bool(m['hom']),
                      het_class='Significant result present' if sig else 'Non-significant results only' if rows else 'No successful het results in snapshot',
                      het_parameters_tested=len({r['parameter_name'] for r in rows if r['parameter_name']}),
                      het_parameters_significant=len({r['parameter_name'] for r in sig if r['parameter_name']}),
                      het_procedures_tested=len({r['procedure_name'] for r in rows if r['procedure_name']}),
                      mouse_symbols=sorted(m['mouse_symbols']), sources=sorted(m['sources']))
        if not m['covered']:
            for k in ('het_rows','het_significant_rows','het_nonsignificant_rows','het_significance_unknown_rows','het_parameters_tested','het_parameters_significant','het_procedures_tested'):
                result[k] = None
            result['het_class'] = 'Heterozygous cohort coverage unavailable'
        GENES[s]['mouse'] = result

    p = source('library', 'Grace_pipeline/261011_RUN_final/mxy_local_eda/outputs/20261101_library_update/gene_summary.csv', 'November library design recovery',
               'Scenario-specific counts from the recovered October export. Pre-PAM SNPs, post-NGG SNPs and selected guides are separate populations. Not a complete source rerun.', date='2026-11-01')
    keep = ['scenario','category','prePAM_SNPs','postNGG_SNPs','selected_SNPs','selected_unique_guides','selected_unique_spacers','editors','selected_provisional_4_8_guides']
    for r in table(p):
        g = gene(symbol=r['gene'])
        if g:
            g.setdefault('library', []).append({k:number(r[k]) if k not in ('scenario','category','editors') else r[k] for k in keep})

    for g in GENES.values():
        g['cmt_assertions'] = [r for r in g['diseases'] if 'charcot-marie-tooth' in (r['DISEASE LABEL'] or '').lower()]
        if g['symbol'] == 'PMP22':
            g['pmp22_regions'] = regions

    # Compact search index plus bounded shards avoids one large request per gene.
    fields = ['symbol','name','hgnc_id','ensembl_id','entrez_id','aliases','dnd_member','cmt','hi_score','shet','phaplo','clinvar_total','dominant_plp']
    index, shards = [], defaultdict(dict)
    selected = {s:g for s,g in GENES.items() if g['diseases'] or g['dosage'] or any(k in g for k in ('shet','collins_records','dnd','mouse','clinvar'))}
    for s,g in sorted(selected.items()):
        shards[s[0] if s[0].isascii() and s[0].isalpha() else '_'][s] = g
        index.append([s,g['name'],g['hgnc_id'],g['ensembl_id'],g['entrez_id'],g['aliases'],'dnd' in g,bool(g['cmt_assertions']),g.get('hi_score'),
                      g.get('shet',{}).get('post_mean'),g.get('collins',{}).get('pHaplo'),g.get('clinvar',{}).get('Total_alleles'),g.get('dnd',{}).get('dominant_plp_mutations')])
    for key, data in shards.items():
        (OUT / f'genes-{key}.json').write_text(json.dumps(data, separators=(',',':'), allow_nan=False))
    (OUT / 'index.json').write_text(json.dumps({'fields':fields,'rows':index}, separators=(',',':')))
    counts = dict(genes=len(selected), dnd=sum('dnd' in g for g in selected.values()),
                  clingen=sum(bool(g['diseases']) for g in selected.values()), dosage=sum(bool(g['dosage']) for g in selected.values()),
                  shet=sum('shet' in g for g in selected.values()), collins=sum('collins' in g for g in selected.values()),
                  mouse=sum('mouse' in g for g in selected.values()), clinvar=sum('clinvar' in g for g in selected.values()))
    manifest = dict(built_at=datetime.now(ZoneInfo('America/Los_Angeles')).isoformat(), counts=counts,
                    sources=SOURCES, dnd_dictionary=dictionary, audit=AUDIT,
                    notes=['No single probability of haplosufficiency is calculated.',
                           'ClinGen %HI is a percentile rank. Collins pHaplo is a probability of haploinsufficiency.',
                           'No evidence (HI=0) is not evidence of haplosufficiency.',
                           'Dominant variants are the supplied D&D dominant P/LP subset; all-class dominant counts are unavailable.',
                           'Targeting counts describe research eligibility, not validated therapeutic efficacy.',
                           'Mouse counts are unique successful statistical rows, not individual animals or human clinical observations.'])
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(json.dumps({'counts':counts,'unmatched':{k:len(v) for k,v in AUDIT['unmatched'].items()},'conflicts':AUDIT['conflicts']}, indent=2))


if __name__ == '__main__':
    main()
