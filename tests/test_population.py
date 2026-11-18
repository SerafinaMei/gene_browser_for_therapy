"""Reconcile every published curve and download against the imported aggregate cache."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import unittest

PROJECT = Path(__file__).resolve().parents[1]
DIST = PROJECT/'dist'
CACHE = PROJECT/'data_sources/population-coverage.json.gz'
with gzip.open(CACHE,'rt') as f:
    SOURCE = json.load(f)
RECORDS = {}
for path in (DIST/'data').glob('population-*.json'):
    RECORDS.update(json.loads(path.read_text()))
GENES = {}
for path in (DIST/'data').glob('genes-*.json'):
    GENES.update(json.loads(path.read_text()))
BY_ID = {g['hgnc_id']:s for s,g in GENES.items()}


class PopulationTests(unittest.TestCase):
    def test_cohort_and_missingness(self):
        self.assertEqual((SOURCE['samples'],SOURCE['haplotypes']),(2504,5008))
        self.assertEqual(len(SOURCE['genes']),560)
        self.assertEqual(sum(len(g['strategies']) for g in SOURCE['genes'].values()),1089)
        self.assertEqual(len(SOURCE['excluded']),3)
        self.assertNotIn('SPTLC1',RECORDS)
        for gene in SOURCE['excluded']:
            self.assertNotIn('strategies',RECORDS[gene]['population_coverage'])
        self.assertEqual({p.name.replace('genes-','population-') for p in (DIST/'data').glob('genes-*.json')},
                         {p.name for p in (DIST/'data').glob('population-*.json')})
        manifest=json.loads((DIST/'data/manifest.json').read_text())
        provenance=next(s for s in manifest['sources'] if s['id']=='population_coverage')
        self.assertEqual(provenance['sha256'],hashlib.sha256(CACHE.read_bytes()).hexdigest())

    def test_all_curves_and_downloads_preserve_source_values(self):
        for source_gene in SOURCE['genes'].values():
            symbol=BY_ID[source_gene['hgnc_id']]
            record=RECORDS[symbol]['population_coverage']
            self.assertEqual(record['haplotypes'],SOURCE['haplotypes'])
            for strategy,source in source_gene['strategies'].items():
                with self.subTest(gene=symbol,strategy=strategy):
                    published=record['strategies'][strategy]
                    self.assertEqual(published['points'],source['points'])
                    self.assertEqual(published['possible'],source['possible'])
                    svg=(DIST/published['plot']).read_text()
                    self.assertIn('<svg',svg)
                    self.assertIn(symbol,svg)
                    with gzip.open(DIST/published['export']['path'],'rt') as f:
                        rows=list(csv.DictReader(f))
                    self.assertEqual(len(rows),len(source['points'])+1)
                    self.assertEqual(int(rows[0]['cumulative_haplotypes']),0)
                    covered=0; guides=0
                    for p,row in zip(source['points'],rows[1:]):
                        self.assertGreater(p['added'],0)
                        covered+=p['added']
                        self.assertEqual(p['covered'],covered)
                        self.assertIn(p['guides']-guides,[1,2] if strategy=='excision' else [1])
                        guides=p['guides']
                        self.assertLessEqual(covered,source['possible'])
                        self.assertLessEqual(source['possible'],5008)
                        self.assertEqual(int(row['cumulative_haplotypes']),covered)
                        self.assertEqual(int(row['new_haplotypes']),p['added'])
                        self.assertEqual(int(row['prioritized_gRNAs']),guides)
                        self.assertEqual(int(row['cohort_haplotypes']),5008)
                        self.assertEqual(int(row['prePAM_possible_haplotypes']),source['possible'])
                        self.assertEqual(row['selected_allele_1'],p['allele1'])
                        self.assertEqual(row['selected_allele_2'],p['allele2'])
                        self.assertAlmostEqual(float(row['coverage_percent']),covered/5008*100)
                        if source_gene['dominant_mutations'] and strategy!='excision':
                            self.assertAlmostEqual(float(row['mutation_specific_model_percent']),min(guides/source_gene['dominant_mutations'],1)*100)
                        else:
                            self.assertEqual(row['mutation_specific_model_percent'],'')

    def test_myh7_reference_values_and_excision_cost(self):
        gene=RECORDS['MYH7']['population_coverage']
        self.assertEqual(gene['dominant_mutations'],149)
        self.assertEqual(gene['strategies']['indels']['possible'],3736)
        self.assertEqual([p['added'] for p in gene['strategies']['indels']['points']],[1206,1206,243,243,144,144,79,79])
        guides=[p['guides'] for p in gene['strategies']['excision']['points']]
        self.assertEqual(guides[:6],[2,4,6,8,9,10])


if __name__=='__main__':
    unittest.main()
