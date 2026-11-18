import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('omim_builder',ROOT/'scripts/build_omim_relationships.py')
builder=importlib.util.module_from_spec(spec);spec.loader.exec_module(builder)

class OmimTests(unittest.TestCase):
    def test_parser_preserves_qualifiers_and_optional_identifiers(self):
        r=builder.parse_phenotype('?A phenotype, with complications, 608030 (3), Autosomal dominant, Autosomal recessive')
        self.assertEqual(r,dict(phenotype='?A phenotype, with complications',phenotype_mim='608030',mapping_key=3,inheritance='Autosomal dominant, Autosomal recessive'))
        self.assertEqual(builder.parse_phenotype('{Susceptibility to a phenotype} (2)')['phenotype_mim'],None)
        with self.assertRaises(ValueError):builder.parse_phenotype('Unparseable field')

    def test_genemap_import_uses_gene_mim_not_symbol_or_phenotype_mim(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'genemap2.txt'
            p.write_text('# Generated: example\n# Chromosome\tCyto Location\tComputed Cyto Location\tMim Number\tPhenotypes\nchr16\t16p11.2\t\t137070\tAmyotrophic lateral sclerosis 6, 608030 (3), Autosomal dominant, Autosomal recessive; Essential tremor, hereditary, 4, 614782 (3), Autosomal dominant\n')
            data=builder.import_genemap(p,'2026-11-18')
            self.assertEqual(set(data['relationships']),{'137070'})
            self.assertEqual([r['phenotype_mim'] for r in data['relationships']['137070']],['608030','614782'])

    def test_shipped_fus_matches_supplied_records_and_distinguishes_gene_ids(self):
        cache=json.loads(builder.CACHE.read_text())
        gene=json.loads((ROOT/'dist/data/omim-F.json').read_text())['FUS']['omim_relationships']
        self.assertEqual(gene['source'],cache['source'])
        self.assertEqual(gene['records'],[dict(r,gene_mim='137070') for r in cache['relationships']['137070']])
        self.assertEqual({r['phenotype_mim'] for r in gene['records']},{'608030','614782'})
        self.assertEqual(gene['records'][0]['inheritance'],'AD, AR')
        self.assertEqual({p.name.replace('genes-','omim-') for p in (ROOT/'dist/data').glob('genes-*.json')},{p.name for p in (ROOT/'dist/data').glob('omim-*.json')})

if __name__=='__main__':unittest.main()
