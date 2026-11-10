"""Reconcile browser download summaries with the original record exports."""
import csv
import gzip
import json
import hashlib
import io
import unittest
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT.parent
LIBRARY = ROOT / 'Grace_pipeline/261011_RUN_final/mxy_local_eda/outputs/20261101_library_update'
DETAILS = {}
for path in (PROJECT / 'dist/data').glob('record-details-*.json'):
    DETAILS.update(json.loads(path.read_text()))


def read(path):
    with (gzip.open if path.suffix == '.gz' else open)(path, 'rt', encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def export_rows(record):
    if record.get('format') != 'columnar-gzip':
        return read(PROJECT/'dist'/record['path'])
    with gzip.open(PROJECT/'dist'/record['path'],'rt') as stream:
        packet=json.load(stream)
    columns=packet['genes'][record['record_key']]
    text=io.StringIO(newline='')
    writer=csv.writer(text)
    writer.writerow(packet['headers'])
    writer.writerows(zip(*columns))
    assert hashlib.sha256(text.getvalue().encode('utf-8-sig')).hexdigest()==record['sha256']
    return [dict(zip(packet['headers'],row)) for row in zip(*columns)]


class RecordTests(unittest.TestCase):
    def test_every_scope_reconciles_to_independent_gene_summary(self):
        symbols = {s.upper(): s for s in DETAILS}
        tested = 0
        for row in read(LIBRARY / 'gene_summary.csv'):
            symbol = symbols[row['gene'].upper()]
            scope = next(s for s in DETAILS[symbol]['targeting']['scopes'] if s['window'] == row['scenario'].split(' ')[0])
            combined = scope['rows'][0]
            for key, source in [('variants','prePAM_SNPs'),('ngg_variants','postNGG_SNPs'),('selected_variants','selected_SNPs'),('guides','selected_unique_guides')]:
                self.assertEqual(combined[key], int(row[source]), (symbol, row['scenario'], key))
            tested += 1
        self.assertEqual(tested, 839)

    def test_cell_line_counts_use_matching_variants_and_deduplicate_guides(self):
        variants = read(LIBRARY / 'snp_list_prePAM.csv')
        guides = read(LIBRARY / 'grna_detail_postNGG.csv.gz')
        for gene in ('AARS1','ABCC8','ACD'):
            for scope in DETAILS[gene]['targeting']['scopes']:
                for summary in scope['rows'][1:]:
                    field = {'WTC':'het_WTC','WTD':'het_WTD','p47':'het_p47_locus_only'}[summary['cell_line']]
                    rows = [r for r in variants if r['gene']==gene and r['window_set'] in (scope['window'],'exonic') and r[field]=='True']
                    self.assertEqual(summary['variants'],len({r['variant_id'] for r in rows}))
                    selected = {r['candidate_association_id'] for r in rows if r['library_candidate_postNGG']=='True'}
                    self.assertEqual(summary['guides'],len({r['guide_id'] for r in guides if r['candidate_association_id'] in selected}))
        lines = {r['cell_line']:r for r in DETAILS['AARS1']['targeting']['scopes'][0]['rows']}
        self.assertEqual(lines['WTC']['variants'],0)
        self.assertEqual(lines['WTD']['selected_variants'],1)
        self.assertEqual(lines['p47']['guides'],10)

    def test_downloads_preserve_records_sequences_and_flags(self):
        for gene in ('AARS1','ABCC8','ACD'):
            for kind, source in [('variants','snp_list_prePAM.csv'),('guides','grna_detail_postNGG.csv.gz')]:
                original = [r for r in read(LIBRARY / source) if r['gene']==gene]
                exported = read(PROJECT/'dist'/DETAILS[gene]['targeting']['exports'][kind]['path'])
                sanitized = [{k:("'"+v if v.startswith(('=','+','-','@','\t','\r')) else v) for k,v in r.items()} for r in original]
                self.assertEqual(exported,sanitized)

    def test_clinvar_exports_are_gene_matched_grch38_records(self):
        for gene,gene_id in [('AARS1','16'),('SPTLC1','10558'),('PMP22','5376')]:
            meta=DETAILS[gene]['clinvar_export']
            rows=export_rows(meta)
            self.assertEqual(len(rows),meta['rows'])
            self.assertEqual(len({r['#AlleleID'] for r in rows}),meta['alleles'])
            self.assertTrue(all(r['Assembly']=='GRCh38' for r in rows))
            self.assertTrue(all(gene_id in r['GeneID'].replace(',', ';').split(';') for r in rows))
            self.assertTrue(all('VariationID' in r and 'ClinicalSignificance' in r for r in rows))
        self.assertNotIn('targeting',DETAILS['SPTLC1'])

    def test_sptlc1_position_matches_cached_ensembl_identity(self):
        cached=json.loads((PROJECT/'data_sources/gene_coordinates.json').read_text())['records'][0]
        gene=json.loads((PROJECT/'dist/data/genes-S.json').read_text())['SPTLC1']
        self.assertEqual(cached['id'],gene['ensembl_id'])
        self.assertEqual(cached['assembly_name'],'GRCh38')
        self.assertEqual(DETAILS['SPTLC1']['coordinates'],'9:92000087-92115491')

    def test_expanded_strategy_counts_match_source_nmd_without_dnd_membership(self):
        run=ROOT/'Grace_pipeline/pipeline_recheck_20261101'
        nmd=read(run/'stage_a/indels/NMD/NMD_induction_summary.csv')
        screens={v['strategy_screen']['source_gene']:v['strategy_screen'] for v in DETAILS.values() if 'strategy_screen' in v}
        self.assertEqual(len(screens),1007)
        for source in nmd:
            rows={r['key']:r for r in screens[source['gene']]['rows']}
            for key,column in [('indel_pre','num_vars_considered'),('indel_nmd','num_vars_inducing_NMD'),('indel_escape','num_vars_escaping_NMD')]:
                self.assertEqual(rows[key]['count'],int(source[column]),(source['gene'],key))
        pmp=json.loads((PROJECT/'dist/data/genes-P.json').read_text())['PMP22']
        self.assertNotIn('dnd',pmp)
        rows={r['key']:r for r in screens['PMP22']['rows']}
        self.assertEqual(rows['indel_pre']['count'],2)
        self.assertEqual(rows['indel_nmd']['count'],0)
        self.assertEqual(rows['indel_escape']['count'],2)

    def test_not_in_run_unassessed_and_tested_zero_remain_distinct(self):
        self.assertNotIn('strategy_screen',DETAILS['SPTLC1'])
        self.assertFalse(DETAILS['SRY']['strategy_screen']['assessed'])
        self.assertIsNone(DETAILS['SRY']['strategy_screen']['rows'][0]['count'])
        mecp=DETAILS['MECP2']['strategy_screen']
        self.assertTrue(mecp['assessed'])
        self.assertEqual(mecp['rows'][0]['count'],0)
        terc={r['key']:r for r in DETAILS['TERC']['strategy_screen']['rows']}
        self.assertEqual(terc['indel_nmd']['status'],'not_assessed')
        self.assertIsNone(terc['indel_nmd']['count'])
        for detail in DETAILS.values():
            for row in detail.get('strategy_screen',{}).get('rows',[]):
                if row['key'] in ['excision','excision_pairs','crisproff']:
                    self.assertEqual(row['status'],'not_run')
                    self.assertIsNone(row['count'])

    def test_expanded_splice_counts_and_pam_flags_reconcile_to_raw_records(self):
        run=ROOT/'Grace_pipeline/pipeline_recheck_20261101/results'
        snps=read(run/'snps_pipeline_prePAM.csv')
        guides=read(run/'all_guide_associations.csv.gz')
        for gene in ['PMP22','BRCA1','SCN1A','MECP2']:
            ids={r['variant_id'] for r in snps if r['gene']==gene and r['window_set']=='1-23'}
            row=next(r for r in DETAILS[gene]['strategy_screen']['rows'] if r['key']=='splice')
            self.assertEqual(row['count'],len(ids))
            self.assertEqual(row['pam'],any(r['gene']==gene and r['variant_id'] in ids and r['passes_10bp_SNP_PAM']=='True' for r in guides))
        exported=read(PROJECT/'dist'/DETAILS['PMP22']['strategy_screen']['exports']['guides']['path'])
        self.assertEqual({r['guide_id'] for r in exported},{r['guide_id'] for r in guides if r['gene']=='PMP22'})


if __name__ == '__main__':
    unittest.main(verbosity=2)
