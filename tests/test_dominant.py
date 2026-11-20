"""Check condition-specific inheritance and reconcile every published current count."""
import csv
import gzip
import hashlib
import json
import sys
import unittest
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import refresh_dominant as pipeline


def record(conditions=('C1',), classification='Pathogenic', review='criteria provided, single submitter', inheritance=None):
    root = ET.Element('VariationArchive', RecordType='classified', VariationID='10', VariationName='example')
    ET.SubElement(root, 'RecordStatus').text = 'current'
    classified = ET.SubElement(root, 'ClassifiedRecord')
    allele = ET.SubElement(classified, 'SimpleAllele', AlleleID='20')
    ET.SubElement(ET.SubElement(allele, 'GeneList'), 'Gene', GeneID='1')
    rcv = ET.SubElement(ET.SubElement(classified, 'RCVList'), 'RCVAccession', Accession='RCV1')
    traits = ET.SubElement(rcv, 'ClassifiedConditionList')
    for condition in conditions:
        ET.SubElement(traits, 'ClassifiedCondition', DB='MedGen', ID=condition).text = condition
    germline = ET.SubElement(ET.SubElement(rcv, 'RCVClassifications'), 'GermlineClassification')
    ET.SubElement(germline, 'Description').text = classification
    ET.SubElement(germline, 'ReviewStatus').text = review
    aggregate=ET.SubElement(ET.SubElement(classified,'Classifications'),'GermlineClassification')
    ET.SubElement(aggregate,'Description').text=classification
    ET.SubElement(aggregate,'ReviewStatus').text=review
    if inheritance:
        assertion = ET.SubElement(ET.SubElement(classified, 'ClinicalAssertionList'), 'ClinicalAssertion', ID='30', ContributesToAggregateClassification='true')
        scv_class = ET.SubElement(assertion, 'Classification')
        ET.SubElement(scv_class, 'GermlineClassification').text = 'Pathogenic'
        ET.SubElement(scv_class, 'ReviewStatus').text = 'criteria provided, single submitter'
        ET.SubElement(ET.SubElement(assertion, 'AttributeSet'), 'Attribute', Type='ModeOfInheritance').text = inheritance
        ET.SubElement(ET.SubElement(assertion, 'TraitSet'), 'XRef', DB='MedGen', ID=conditions[0])
    return root


class InheritanceTests(unittest.TestCase):
    def evaluate(self, r, mapping=None, aliases=None):
        return pipeline.evaluate_record(r, {'1': 'GENE'}, mapping or {}, aliases or {})

    def test_exact_condition_join_not_gene_wide_inheritance(self):
        mapping = {('GENE', 'MedGen:C1'): {'AD'}, ('GENE', 'MedGen:C2'): {'AR'}}
        self.assertEqual(self.evaluate(record(), mapping)[1][0][3], 'dominant')
        self.assertIsNone(self.evaluate(record(('C2',)), mapping))
        self.assertEqual(self.evaluate(record(('C3',)), mapping)[1][0][3], 'inheritance_unresolved')

    def test_mixed_inheritance_requires_explicit_variant_evidence(self):
        mapping = {('GENE', 'MedGen:C1'): {'AD', 'AR'}}
        self.assertEqual(self.evaluate(record(), mapping)[1][0][3], 'inheritance_unresolved')
        self.assertEqual(self.evaluate(record(inheritance='Autosomal dominant inheritance'), mapping)[1][0][3], 'dominant')
        self.assertIsNone(self.evaluate(record(inheritance='Autosomal recessive inheritance'), mapping))

    def test_unknown_second_condition_cannot_inherit_first_mapping(self):
        mapping = {('GENE', 'MedGen:C1'): {'AD'}}
        self.assertEqual(self.evaluate(record(('C1', 'C2')), mapping)[1][0][3], 'inheritance_unresolved')

    def test_x_linked_unspecified_is_unresolved(self):
        self.assertEqual(self.evaluate(record(), {('GENE','MedGen:C1'):{'XL'}})[1][0][3], 'inheritance_unresolved')

    def test_conflicts_and_missing_assertion_criteria_are_excluded(self):
        self.assertIsNone(self.evaluate(record(classification='Conflicting classifications of pathogenicity')))
        self.assertIsNone(self.evaluate(record(review='no assertion criteria provided')))
        self.assertIsNone(self.evaluate(record(classification='Uncertain significance')))
        r=record();r.find('ClassifiedRecord/Classifications/GermlineClassification/Description').text='Conflicting classifications of pathogenicity'
        self.assertIsNone(self.evaluate(r))

    def test_multigene_and_obsolete_records_excluded(self):
        r = record()
        ET.SubElement(r.find('ClassifiedRecord/SimpleAllele/GeneList'), 'Gene', GeneID='2')
        self.assertIsNone(self.evaluate(r))
        r = record(); r.find('RecordStatus').text = 'replaced'
        self.assertIsNone(self.evaluate(r))

    def test_explicit_inheritance_must_match_condition(self):
        r = record(inheritance='Autosomal dominant inheritance')
        r.find('ClassifiedRecord/ClinicalAssertionList/ClinicalAssertion/TraitSet/XRef').set('ID', 'C2')
        self.assertEqual(self.evaluate(r)[1][0][3], 'inheritance_unresolved')

    def test_condition_alias_join_is_gene_specific(self):
        aliases = {('GENE','MedGen:C1'):{'OMIM:123456'}}
        self.assertEqual(self.evaluate(record(), {('GENE','OMIM:123456'):{'AD'}}, aliases)[1][0][3], 'dominant')
        self.assertEqual(self.evaluate(record(), {('OTHER','OMIM:123456'):{'AD'}}, aliases)[1][0][3], 'inheritance_unresolved')

    def test_incomplete_xml_never_produces_a_snapshot(self):
        import tempfile
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'partial.xml'
            path.write_text('<ClinVarVariationRelease ReleaseDate="2026-11-15">'+ET.tostring(record(),encoding='unicode'))
            with self.assertRaises(ET.ParseError):
                pipeline.parse_release(path, ({},{'1':'GENE'},{},{}), minimum_records=1)
            path.write_text(path.read_text()+'</ClinVarVariationRelease>')
            with self.assertRaisesRegex(ValueError,'Incomplete ClinVar release'):
                pipeline.parse_release(path, ({},{'1':'GENE'},{},{}), minimum_records=2)


class PublishedSnapshotTests(unittest.TestCase):
    def test_all_counts_exports_and_browse_index_reconcile(self):
        with gzip.open(pipeline.CACHE, 'rt') as f:
            cache = json.load(f)
        self.assertGreater(cache['source']['records'], 1000000)
        self.assertEqual(cache['headers'], pipeline.HEADERS)
        index = json.loads((ROOT / 'dist/data/index.json').read_text())
        index_counts = {r[index['fields'].index('symbol')]: r[index['fields'].index('dominant_plp')] for r in index['rows']}
        checked = 0
        for path in (ROOT / 'dist/data').glob('dominant-*.json'):
            letter = path.stem.removeprefix('dominant-')
            with gzip.open(ROOT / f'dist/downloads/dominant/{letter}.json.gz', 'rt') as f:
                export = json.load(f)
            genes = json.loads((ROOT / f'dist/data/genes-{letter}.json').read_text())
            shard = json.loads(path.read_text())
            self.assertEqual(set(shard), set(genes))
            for symbol, wrapped in shard.items():
                count = wrapped['dominant_current']
                rows = cache['genes'].get(symbol, [])
                dominant = {r[0] for r in rows if r[3] == 'dominant'}
                unresolved = {r[0] for r in rows if r[3] == 'inheritance_unresolved'} - dominant
                self.assertEqual(count['count'], len(dominant) if genes[symbol].get('entrez_id') else None, symbol)
                self.assertEqual(count['unresolved'], len(unresolved) if genes[symbol].get('entrez_id') else None, symbol)
                self.assertEqual(count['date'], cache['source']['date'])
                self.assertEqual(index_counts[symbol], count['count'], symbol)
                if rows:
                    self.assertTrue(all(isinstance(value,str) for row in rows for value in row),symbol)
                    self.assertEqual(rows, [list(r) for r in zip(*export['genes'][symbol])])
                checked += 1
        self.assertEqual(checked, len(index_counts))
        source = next(s for s in json.loads((ROOT / 'dist/data/manifest.json').read_text())['sources'] if s['id'] == 'dominant_current')
        self.assertEqual(source['sha256'], hashlib.sha256(pipeline.CACHE.read_bytes()).hexdigest())
        self.assertEqual(source['date'], cache['source']['date'])

    def test_current_comparison_preserves_coverage_and_updates_only_model(self):
        checked=0
        for path in (ROOT/'dist/data').glob('population-*.json'):
            letter=path.stem.removeprefix('population-')
            dominant=json.loads((ROOT/f'dist/data/dominant-{letter}.json').read_text())
            for symbol,wrapper in json.loads(path.read_text()).items():
                population=wrapper['population_coverage']
                for strategy,series in population.get('strategies',{}).items():
                    count=dominant[symbol]['dominant_current']
                    if strategy=='excision' or not count['assessed']:
                        self.assertNotIn('current_comparison',series)
                        continue
                    current=series['current_comparison']
                    self.assertEqual(current['count'],count['count'])
                    self.assertEqual(current['date'],count['date'])
                    with gzip.open(ROOT/'dist'/current['export']['path'],'rt') as f:rows=list(csv.DictReader(f))
                    with gzip.open(ROOT/'dist'/series['export']['path'],'rt') as f:original=list(csv.DictReader(f))
                    self.assertEqual(len(rows),len(original))
                    for row,old in zip(rows,original):
                        for field in old.keys()-{'mutation_specific_model_percent'}:
                            self.assertEqual(row[field],old[field],(symbol,field))
                        if count['count']:
                            self.assertAlmostEqual(float(row['mutation_specific_model_percent']),min(int(row['prioritized_gRNAs'])/count['count'],1)*100)
                        else:
                            self.assertEqual(row['mutation_specific_model_percent'],'')
                        self.assertEqual(row['ClinVar_release'],count['date'])
                    self.assertTrue((ROOT/'dist'/current['plot']).is_file())
                    checked+=1
        self.assertGreater(checked,500)


if __name__ == '__main__':
    unittest.main()
