"""Independent comparisons against original inputs and semantic invariants."""
import csv
import gzip
import json
import unittest
from pathlib import Path
import openpyxl

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT.parent
GENES = {}
for p in (PROJECT / 'dist/data').glob('genes-*.json'):
    GENES.update(json.loads(p.read_text()))


class SourceTests(unittest.TestCase):
    def test_shet_direct_from_hgnc_id(self):
        w = openpyxl.load_workbook(ROOT/'DnD_project/data/Spense_S-het_supp_table.xlsx',read_only=True,data_only=True)
        rows = w['Supplementary Table 1'].values
        keys = next(rows)
        source = {r[1]:dict(zip(keys,r)) for r in rows}
        for s in ('AARS1','PMP22','SPTLC1','BRCA1','ABCB6'):
            self.assertEqual(GENES[s]['shet'],source[GENES[s]['hgnc_id']])
        w.close()

    def test_collins_is_original_probability_not_complement(self):
        w=openpyxl.load_workbook(ROOT/'DnD_project/data/Cell_paper/Cell_Haplo_pred.xlsx',read_only=True,data_only=True)
        raw={r[0]:r[1] for r in list(w['Table S7'].values)[1:]}
        for s in ('AARS1','PMP22','SPTLC1','BRCA1'):
            self.assertEqual(GENES[s]['collins']['pHaplo'],raw[GENES[s]['collins']['Gene']])
        self.assertGreater(GENES['PMP22']['collins']['pHaplo'],.9)
        w.close()

    def test_all_clingen_assertions_preserved(self):
        p=ROOT/'DnD_project/data/Clingen-Gene-Disease-Summary-2025-11-13.csv'
        with p.open() as f:
            for _ in range(4):next(f)
            rows=[r for r in csv.DictReader(f) if r['GENE ID (HGNC)'].startswith('HGNC:')]
        self.assertEqual(len(rows),sum(len(g['diseases']) for g in GENES.values()))
        a=[r for r in rows if r['GENE ID (HGNC)']=='HGNC:20']
        self.assertEqual(GENES['AARS1']['diseases'],a)

    def test_duplicate_disease_rows_do_not_multiply_variant_counts(self):
        w=openpyxl.load_workbook(ROOT/'DnD_project/data/dnd_593_2025_12_23_w_dplp_muts.xlsx',read_only=True,data_only=True)
        it=w.worksheets[0].values;keys=next(it)
        raw=[dict(zip(keys,r)) for r in it]
        rows=[r for r in raw if r['hgnc_symbol']=='AARS1']
        self.assertGreater(len(rows),1)
        self.assertEqual(GENES['AARS1']['dnd']['dominant_plp_mutations'],4)
        self.assertEqual(GENES['AARS1']['dnd']['total_plp'],87)
        self.assertEqual(sum('dnd' in g for g in GENES.values()),593)
        w.close()

    def test_clinvar_total_is_alleles_not_submissions_or_plp(self):
        with (PROJECT/'data_sources/clinvar_gene_summary.txt').open() as f:
            next(f);raw=list(csv.DictReader(f,delimiter='\t'))
        for s in ('AARS1','PMP22','SPTLC1','AKAP17A'):
            r=next(r for r in raw if r['GeneID']==GENES[s]['entrez_id'])
            self.assertEqual(GENES[s]['clinvar']['Total_alleles'],int(r['Total_alleles']))
        self.assertNotEqual(GENES['AARS1']['clinvar']['Total_alleles'],GENES['AARS1']['dnd']['total_plp'])
        self.assertEqual(len(GENES['AKAP17A']['clinvar_records']),2)
        self.assertEqual(GENES['AKAP17A']['clinvar']['Total_alleles'],175)

    def test_mouse_deduplicated_successful_rows_direct(self):
        rows={}
        for cohort in ('DND593','Dominant855','HI417','AR754'):
            p=ROOT/f'DnD_project/MousePhenotype/{cohort}_heterozygous_statistical_results_all_genes.txt.gz'
            with gzip.open(p,'rt') as f:
                for r in csv.DictReader(f,delimiter='\t'):
                    if r['marker_accession_id']=='MGI:2384560' and r['status']=='Successful' and r['zygosity']=='heterozygote':rows[r['doc_id']]=r
        m=GENES['AARS1']['mouse']
        self.assertEqual(m['het_rows'],len(rows))
        self.assertEqual(m['het_significant_rows'],sum(r['significant']=='True' for r in rows.values()))
        self.assertEqual(m['het_parameters_tested'],len({r['parameter_name'] for r in rows.values()}))

    def test_gene_wide_semantics_and_missingness(self):
        for g in GENES.values():
            if 'collins' in g:
                self.assertTrue(0<=g['collins']['pHaplo']<=1)
            if g.get('hi_percentile') is not None:self.assertTrue(0<=g['hi_percentile']<=100)
            if 'shet' in g:
                s=g['shet'];self.assertLessEqual(s['post_lower_95'],s['post_mean']);self.assertLessEqual(s['post_mean'],s['post_upper_95'])
            if 'mouse' in g and g['mouse']['het_rows'] is not None:
                m=g['mouse'];self.assertEqual(m['het_rows'],m['het_significant_rows']+m['het_nonsignificant_rows']+m['het_significance_unknown_rows'])
                self.assertLessEqual(m['het_parameters_tested'],m['het_rows'])
        self.assertNotIn('dnd',GENES['PMP22'])
        self.assertEqual(GENES['PMP22']['hi_score'],3)
        self.assertEqual(GENES['AARS1']['hi_score'],0)
        self.assertNotIn('mouse',GENES['SPTLC1'])
        self.assertGreater(len(GENES['LCOR']['collins_records']),1)
        self.assertNotIn('collins',GENES['LCOR'])


if __name__=='__main__':unittest.main(verbosity=2)
