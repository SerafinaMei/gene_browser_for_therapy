import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {searchGenes,unpackIndex,flatten,csv,fmt,present,esc,hiLabels,genomeUrl} from '../dist/core.js';
const data=JSON.parse(readFileSync(new URL('../dist/data/index.json',import.meta.url)));
const genes=unpackIndex(data);
test('case-insensitive symbols and exact HGNC, Ensembl, versioned Ensembl, NCBI IDs resolve',()=>{
  for(const q of ['aars1','HGNC:20','hgnc:20','ENSG00000090861','ENSG00000090861.12','16'])assert.equal(searchGenes(genes,q)[0].symbol,'AARS1');
  assert.equal(searchGenes(genes,'PMP22')[0].symbol,'PMP22');
});
test('aliases resolve and approved symbols rank before shared aliases',()=>{
  assert.equal(searchGenes(genes,'CMT1A')[0].symbol,'PMP22');
  const fixtures=[{symbol:'ABC',name:'one',aliases:['XYZ']},{symbol:'XYZ',name:'two',aliases:[]}];
  assert.equal(searchGenes(fixtures,'XYZ')[0].symbol,'XYZ');
  assert.equal(searchGenes(genes,'no_such_gene_123').length,0);
  assert.equal(searchGenes(genes,'').length,0);
});
test('missing data, false and zero are distinct',()=>{
  assert.equal(present(0),true);assert.equal(present(false),true);assert.equal(present(null),false);
  assert.equal(fmt(0),'0');assert.equal(fmt(null),'Not available');
  assert.equal(hiLabels[0],'No');assert.equal(hiLabels[40],'Dosage sensitivity unlikely');
});
test('CSV preserves scientific values and defuses spreadsheet formulas',()=>{
  const out=csv([['=CMD()',0,null,'a,"b"\nc',0.00083437]]);
  assert.ok(out.startsWith('"\'=CMD()","0","","a,""b""\nc","0.00083437"'));
});
test('raw export keeps every disease row, arrays, false, null and raw zero',()=>{
  assert.deepEqual(flatten({diseases:[{MOI:'AD'},{MOI:'AR'}],mouse:{has_hom_record:false,het_rows:0},score:null}),{'diseases[0].MOI':'AD','diseases[1].MOI':'AR','mouse.has_hom_record':false,'mouse.het_rows':0,score:null});
});
test('source strings are escaped before HTML rendering',()=>{
  assert.equal(esc('<img src=x onerror="bad()">'), '&lt;img src=x onerror=&quot;bad()&quot;&gt;');
  assert.equal(esc('No evidence — uncertain – review'), 'No evidence - uncertain - review');
});
test('Genome Tracks overrides the saved location and never uses a temporary session ID',()=>{
  for (const g of [{symbol:'SPTLC1'}, {symbol:'AARS1',coordinates:'16:70252298-70289506'}, {symbol:'PMP22',coordinates:'chr17:15324739-15360262'}]) {
    const url = new URL(genomeUrl(g, true));
    assert.equal(url.searchParams.get('db'),'hg38');
    assert.equal(url.searchParams.get('position'),g.coordinates ? 'chr'+g.coordinates.replace(/^chr/,'') : g.symbol);
    assert.equal(url.searchParams.get('ignoreCookie'),'1');
    assert.equal(url.searchParams.get('hgS_otherUserName'),'gramey02');
    assert.equal(url.searchParams.get('hgS_otherUserSessionName'),'SPTLC_gene_session');
    assert.equal(url.searchParams.has('hgsid'),false);
    assert.equal(url.searchParams.has('hubUrl'),false);
  }
  assert.equal(new URL(genomeUrl({symbol:'SPTLC1'})).searchParams.get('singleSearch'),'knownCanonical');
});
