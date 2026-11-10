"""Verify the broad import, deduplication, and source/UI reconciliation."""
import csv,gzip,hashlib,json,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'dist'
class ComprehensiveEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.genes={}
        for p in (OUT/'data').glob('evidence-?.json'):cls.genes.update(json.loads(p.read_text()))
    def test_mouse_snapshot_complete_and_sptlc1(self):
        manifest=json.loads((ROOT/'data_sources/comprehensive/mouse-sources.json').read_text())
        for source in manifest['sources']:
            path=ROOT/source['path']
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),source['sha256'])
            with gzip.open(path,'rt') as f:self.assertEqual(len(json.load(f)),source['records'])
        m=self.genes['SPTLC1']['mouse_current']
        self.assertEqual(m['unique_terms'],3)
        self.assertEqual(m['supporting_datasets'],5)
        self.assertEqual(len(m['diseases']),4)
        self.assertEqual({p['term'] for p in m['phenotypes'] if p['zygosity']=='heterozygote'}, {'abnormal embryo size','embryonic growth retardation'})
        self.assertTrue(all(x['assertion_type']=='manual' for p in m['phenotypes'] for x in p['support']))
    def test_every_gene_count_matches_downloaded_alleles(self):
        checked=0
        datasets=[(s,g[k]) for s,g in self.genes.items() for k in ['genome_targets','other_genome_targets','pipeline_targets'] if k in g]
        for symbol,t in datasets:
            if t['export'] is None:
                self.assertEqual(t['variants'],0);self.assertTrue(all(n==0 for n in t['strategies'].values()));continue
            raw=gzip.decompress((OUT/t['export']['path']).read_bytes())
            self.assertEqual(hashlib.sha256(raw).hexdigest(),t['export']['sha256'],symbol)
            rows=list(csv.DictReader(raw.decode().splitlines()))
            keys={(r['chromosome'],r['start_0based'],r['end'],r['alleles']) for r in rows}
            self.assertEqual(t['variants'],len(keys),symbol)
            self.assertEqual(t['export']['rows'],len(rows),symbol)
            grouped={}
            for r in rows:
                key=tuple(r[k] for k in ['chromosome','start_0based','end','alleles'])
                if key not in grouped:grouped[key]=dict(r)
                else:
                    old=grouped[key];old['strategies']=','.join(sorted(set(old['strategies'].split(','))|set(r['strategies'].split(','))))
                    if r['PAM_targetable']=='Yes':old['PAM_targetable']='Yes'
            rows=list(grouped.values())
            self.assertEqual(t['pam_variants'],sum(r['PAM_targetable']=='Yes' for r in rows),symbol)
            for strategy,n in t['strategies'].items():self.assertEqual(n,sum(strategy in [x.strip() for x in r['strategies'].split(',')] for r in rows),symbol)
            for c in t['cells']:
                matched=[r for r in rows if r[c['cell_line']+'_quality']]
                high=[r for r in matched if 'HIGH' in r[c['cell_line']+'_quality'].split('|')]
                self.assertEqual(c['variants'],len(matched),symbol)
                self.assertEqual(c['high_quality'],len(high),symbol)
                self.assertEqual(c['low_only'],len(matched)-len(high),symbol)
                self.assertEqual(c['pam_high_quality'],sum(r['PAM_targetable']=='Yes' for r in high),symbol)
            checked+=1
        self.assertGreater(checked,18000)
    def test_multiallelic_summary_and_all_gene_coverage(self):
        audit=json.loads((OUT/'data/run-multiallelic-validation.json').read_text())
        self.assertEqual(audit['status'],'PASS')
        self.assertEqual((audit['summary_genes'],audit['genes_with_candidates'],audit['records']),(593,568,135549))
        self.assertEqual(sum('pipeline_targets' in g for g in self.genes.values()),593)
        self.assertEqual(self.genes['RARB']['pipeline_targets']['strategies']['excision'],2427)
        self.assertEqual(self.genes['RYR2']['pipeline_targets']['strategies']['excision'],2201)
        self.assertNotIn('pipeline_targets',self.genes['SPTLC1'])
        by_source={original:g['pipeline_targets'] for g in self.genes.values() if 'pipeline_targets' in g for original in g['pipeline_targets']['source_symbols']}
        fields={'num_indel_vars':'exon disruption','num_crisproff_vars':'epigenetic silencing','num_excision_vars':'excision','num_base_edit_snps':'splice site disruption'}
        with (ROOT/'data_sources/comprehensive/run_multiallelic_summary.csv').open() as f:
            for row in csv.DictReader(f):
                for source,strategy in fields.items():self.assertEqual(by_source[row['hgnc_symbol']]['strategies'][strategy],int(float(row[source] or 0)))

    def test_sptlc_specific_track_regression(self):
        t=self.genes['SPTLC1']['genome_targets']
        self.assertEqual(t['variants'],434)
        self.assertEqual(t['pam_variants'],312)
        self.assertEqual({r['cell_line']:r['variants'] for r in t['cells']},{'WTB':191,'WTC':1,'WTD':121,'KOLF2':126})
        self.assertEqual(self.genes['SPTLC1']['other_genome_targets']['variants'],47)

    def test_current_clinvar_browse_consistency(self):
        idx=json.loads((OUT/'data/index.json').read_text());cols=idx['fields']
        for r in idx['rows']:
            g=self.genes.get(r[cols.index('symbol')],{})
            if 'clinvar_date' in g:self.assertEqual(r[cols.index('clinvar_total')],(g['clinvar'] or {}).get('Total_alleles'))
        self.assertEqual(self.genes['SPTLC1']['clinvar_date'],'2026-11-09')
if __name__=='__main__':unittest.main()
