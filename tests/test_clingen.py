"""Reconcile the shipped disease assertions with the official cached export."""
import csv
import hashlib
import json
import unittest
from collections import defaultdict
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
CACHE = PROJECT/'data_sources/clingen-gene-validity.csv'
GENES = {}
for path in (PROJECT/'dist/data').glob('genes-*.json'):
    GENES.update(json.loads(path.read_text()))


class ClinGenTests(unittest.TestCase):
    def test_every_source_assertion_is_preserved_under_hgnc_id(self):
        with CACHE.open() as f:
            for _ in range(4):
                next(f)
            rows = [{k:v.strip() or None for k,v in r.items()} for r in csv.DictReader(f)
                    if r['GENE ID (HGNC)'].startswith('HGNC:')]
        expected = defaultdict(list)
        for row in rows:
            expected[row['GENE ID (HGNC)']].append(row)
        self.assertEqual(len(rows),sum(len(g['diseases']) for g in GENES.values()))
        self.assertEqual(set(expected),{g['hgnc_id'] for g in GENES.values() if g['diseases']})
        for gene in GENES.values():
            self.assertEqual(gene['diseases'],expected[gene['hgnc_id']])

    def test_snapshot_provenance_and_index_match_records(self):
        manifest=json.loads((PROJECT/'dist/data/manifest.json').read_text())
        source=next(s for s in manifest['sources'] if s['id']=='validity')
        self.assertEqual(source['sha256'],hashlib.sha256(CACHE.read_bytes()).hexdigest())
        self.assertIn('FILE CREATED: '+source['date'],CACHE.read_text())
        self.assertEqual(manifest['counts']['clingen'],sum(bool(g['diseases']) for g in GENES.values()))
        index=json.loads((PROJECT/'dist/data/index.json').read_text())
        for values in index['rows']:
            row=dict(zip(index['fields'],values))
            self.assertEqual(row['cmt'],bool(GENES[row['symbol']]['cmt_assertions']))

    def test_validity_refresh_preserves_dosage_records(self):
        # A new disease assertion must not manufacture a dosage-sensitivity score.
        import importlib.util
        import tempfile
        from unittest.mock import patch
        spec=importlib.util.spec_from_file_location('refresh_clingen',PROJECT/'scripts/refresh_clingen.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as folder:
            out=Path(folder)
            original={'symbol':'GENE','hgnc_id':'HGNC:1','dosage':[{'score':3}], 'diseases':[], 'cmt_assertions':[]}
            (out/'genes-G.json').write_text(json.dumps({'GENE':original}))
            (out/'index.json').write_text(json.dumps({'fields':['symbol','cmt'],'rows':[['GENE',False]]}))
            (out/'manifest.json').write_text(json.dumps({'sources':[{'id':'validity'}],'counts':{},'audit':{'unmatched':{}}}))
            assertion={'GENE ID (HGNC)':'HGNC:1','DISEASE LABEL':'Example','DISEASE ID (MONDO)':'MONDO:1','MOI':'AD'}
            with patch.object(module,'OUT',out), patch.object(module,'read_validity',return_value=('2026-11-20',[assertion])), patch.object(module,'source_metadata',return_value={'id':'validity'}):
                module.refresh()
            updated=json.loads((out/'genes-G.json').read_text())['GENE']
            self.assertEqual(updated['diseases'],[assertion])
            self.assertEqual(updated['dosage'],original['dosage'])


if __name__=='__main__':
    unittest.main()
