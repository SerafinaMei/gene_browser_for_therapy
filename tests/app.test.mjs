// Application-flow smoke tests with a minimal DOM adapter, not browser/visual QA.
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {createHash} from 'node:crypto';
import {gzipSync} from 'node:zlib';
const elements=new Map(),documentEvents={},windowEvents={},downloads=[];
class Element {
  constructor(){this.innerHTML='';this.textContent='';this.value='';this.open=false;this.attrs={};this.listeners={};this.classList={add(){},remove(){},toggle(){}};}
  setAttribute(k,v){this.attrs[k]=v;}
  getAttribute(k){return this.attrs[k];}
  removeAttribute(k){delete this.attrs[k];}
  addEventListener(k,v){this.listeners[k]=v;}
  showModal(){this.open=true;}
  close(){this.open=false;}
  focus(){}
}
globalThis.document={body:{dataset:{}},querySelector(s){if(!elements.has(s))elements.set(s,new Element());return elements.get(s);},querySelectorAll(){return [];},addEventListener(k,v){documentEvents[k]=v;},activeElement:{tagName:'BODY'},createElement(){return {click(){downloads.push(this);}};}};
globalThis.window={addEventListener(k,v){windowEvents[k]=v;}};
globalThis.localStorage={getItem(){return '{"malformed":"saved selection"}';},setItem(){}};
globalThis.location={hash:'',reload(){}};
globalThis.fetch=async path=>({ok:true,json:async()=>JSON.parse(await readFile(new URL(`../dist/${path}`,import.meta.url),'utf8')),blob:async()=>new Blob([await readFile(new URL(`../dist/${path}`,import.meta.url))])});
await import('../dist/app.js');
const main=()=>elements.get('#main').innerHTML;
async function route(hash){location.hash=hash;await windowEvents.hashchange();return main();}
test('initial page opens search without a preselected gene',()=>{
  assert.equal(document.body.dataset.view,'home');
  assert.match(main(),/28,678 genes available/);
  assert.doesNotMatch(main(),/AARS1|Dosage and constraint/);
  assert.doesNotMatch(main(),/could not load/);
});
test('submitting a gene ID opens its report and Search returns to the home page',async()=>{
  elements.get('#gene-search').value='HGNC:20';
  elements.get('#search-form').listeners.submit({preventDefault(){}});
  const hash=location.hash.startsWith('#')?location.hash:'#'+location.hash;
  assert.match(await route(hash),/AARS1/);
  const evidence=JSON.parse(await readFile(new URL('../dist/data/evidence-A.json',import.meta.url),'utf8'));assert.ok(main().includes(evidence.AARS1.clinvar.Total_alleles.toLocaleString('en-US')));assert.match(main(),/0\.02246/);
  assert.equal(document.body.dataset.view,'report');
  elements.get('#search-nav').onclick();
  await route('#view=search');
  assert.equal(document.body.dataset.view,'home');
  assert.equal(elements.get('#gene-search').value,'');
});
test('all report sections render for a complete record and a sparse record',async()=>{
  for(const symbol of ['AARS1','SPTLC1','SEC22B','TCAF2C'])for(const tab of ['overview','diseases','variants','mouse','raw']){
    const html=await route(`#gene=${symbol}&tab=${tab}`);
    assert.match(html,new RegExp(symbol));assert.doesNotMatch(html,/could not load|undefined|NaN/);
  }
});
test('PMP22 report preserves region distinction and its curated HI evidence',async()=>{
  const html=await route('#gene=PMP22&tab=diseases');
  assert.match(html,/Sufficient evidence/);assert.match(html,/recurrent 17p12/);assert.match(html,/PMP22-associated neuropathies/);
});
test('stable ID URL, invalid symbol, source page and browsing render without exceptions',async()=>{
  assert.match(await route('#gene=HGNC:20'),/AARS1/);
  assert.match(await route('#gene=NO_SUCH_GENE_123'),/Gene not found/);
  assert.match(await route('#view=sources'),/GeneBayes/);
  assert.match(await route('#view=browse'),/28,678/);
  assert.match(elements.get('#browse-results').innerHTML,/ClinVar alleles/);
});
test('plain evidence label, D&D membership and gene-specific Genome Tracks links',async()=>{
  let html=await route('#gene=AARS1');
  assert.match(html,/<td>D&D<\/td><td>Yes<\/td>/);
  assert.match(html,/<td>ClinGen haploinsufficiency evidence<\/td><td>No<\/td>/);
  assert.match(html,/PMC13042096/);
  assert.doesNotMatch(html,/D&D gene set|—/);
  html=await route('#gene=PMP22');
  assert.match(html,/<td>D&D<\/td><td>No<\/td>/);
  html=await route('#gene=SPTLC1');
  assert.match(html,/position=chr9%3A92000087-92115491/);
  assert.match(html,/>Genome Tracks/);
  assert.doesNotMatch(html,/hgsid|hubUrl/);
});
test('cell-line view replaces internal scenarios and offers record downloads',async()=>{
  const html=await route('#gene=AARS1&tab=variants');
  assert.match(html,/Cell-line targeting/);
  assert.match(html,/WTC/);assert.match(html,/WTD/);assert.match(html,/p47/);
  assert.match(html,/p47 confirms positions only/);
  assert.match(html,/AARS1-guides.csv.gz/);
  assert.match(html,/AARS1-variants.csv.gz/);
  assert.match(html,/downloads\/clinvar-packed\/A-\d+\.json.gz/);
  assert.doesNotMatch(html,/November|Recovered candidate|Scenario|—/);
  documentEvents.change({target:{id:'target-window',value:'1-30'}});
  assert.match(elements.get('#cell-line-results').innerHTML,/WTD/);
});
test('download action decompresses actual CSV records with guide sequences',async()=>{
  const blobs=[];
  const originalCreate=URL.createObjectURL;
  URL.createObjectURL=blob=>{blobs.push(blob);return 'blob:test';};
  try {
    const button={dataset:{recordFile:'downloads/targeting/AARS1-guides.csv.gz',filename:'AARS1-grna-sequences.csv'},textContent:'gRNA sequences (CSV)',hasAttribute(){return false;}};
    const target={closest:selector=>selector==='button,[data-result]'?button:null};
    await documentEvents.click({target});
    assert.equal(button.disabled,false);
    assert.equal(downloads.at(-1).download,'AARS1-grna-sequences.csv');
    const text=await blobs.at(-1).text();
    assert.match(text,/gRNA,PAM/);
    assert.match(text,/AARS1/);
    assert.match(text,/[ACGT]{20},[ACGT]GG/);
  } finally {URL.createObjectURL=originalCreate;}
});
test('expanded strategy coverage shows HI results and distinguishes missing assessments',async()=>{
  let html=await route('#gene=PMP22&tab=variants');
  assert.match(html,/Exonic indel, before NMD assessment<\/td><td>2<\/td>/);
  assert.match(html,/Exonic indel, predicted NMD<\/td><td>0<\/td>/);
  assert.match(html,/Not run in this dataset/);
  assert.match(html,/downloads\/expanded\/PMP22-guides.csv.gz/);
  assert.match(html,/Expanded indel\/splice screen/);
  html=await route('#gene=SPTLC1&tab=variants');
  assert.match(html,/Genome-wide targeting candidates/);assert.match(html,/KOLF2/);assert.doesNotMatch(html,/SPTLC1 was not included/);
  assert.doesNotMatch(html,/Exonic indel, before NMD assessment/);
  html=await route('#gene=MECP2&tab=variants');
  assert.match(html,/Exonic indel, before NMD assessment<\/td><td>0<\/td>/);
  html=await route('#gene=SRY&tab=variants');
  assert.match(html,/No chromosome Y allele-frequency input/);
  await route('#gene=AARS1&tab=variants');
  documentEvents.change({target:{id:'strategy-dataset',value:'expanded'}});
  assert.match(elements.get('#strategy-results').innerHTML,/Not run in this dataset/);
  documentEvents.change({target:{id:'strategy-dataset',value:'reference'}});
  assert.match(elements.get('#strategy-results').innerHTML,/D&D reference dataset/);
});
test('packed ClinVar downloads reconstruct the original CSV bytes',async()=>{
  const originalCreate=URL.createObjectURL;
  const blobs=[];
  URL.createObjectURL=blob=>{blobs.push(blob);return 'blob:clinvar-test';};
  try {
    for(const symbol of ['AARS1','SPTLC1','PMP22','TTN']) {
      const detail=JSON.parse(await readFile(new URL(`../dist/data/record-details-${symbol[0]}.json`,import.meta.url),'utf8'))[symbol].clinvar_export;
      const button={dataset:{recordFile:detail.path,recordFormat:detail.format,recordKey:detail.record_key,filename:`${symbol}.csv`},textContent:'ClinVar variants (CSV)',hasAttribute(){return false;}};
      await documentEvents.click({target:{closest:selector=>selector==='button,[data-result]'?button:null}});
      const bytes=Buffer.from(await blobs.at(-1).arrayBuffer());
      assert.equal(createHash('sha256').update(bytes).digest('hex'),detail.sha256,symbol);
      assert.equal(downloads.at(-1).download,`${symbol}.csv`);
      assert.equal(button.disabled,false);
    }
  } finally {URL.createObjectURL=originalCreate;}
});
test('hosted reports load compressed annotation shards',async()=>{
  const originalFetch=globalThis.fetch;
  document.documentElement={dataset:{jsonEncoding:'gzip',buildVersion:'test-revision'}};
  const requested=[];
  globalThis.fetch=async path=>{
    requested.push(path);
    assert.match(path,/\.json\.gz\?v=test-revision$/);
    const raw=await readFile(new URL(`../dist/${path.split('?')[0].replace(/\.gz$/,'')}`,import.meta.url));
    return {ok:true,blob:async()=>new Blob([gzipSync(raw)])};
  };
  try {
    assert.match(await route('#gene=BRCA1&tab=variants'),/BRCA1/);
    assert.doesNotMatch(main(),/could not load/);
    assert.ok(requested.includes('data/genes-B.json.gz?v=test-revision'));
    assert.ok(requested.includes('data/record-details-B.json.gz?v=test-revision'));
  } finally {globalThis.fetch=originalFetch;delete document.documentElement;}
});

test('current IMPC includes manual heterozygous embryo phenotypes and disease details',async()=>{
  const html=await route('#gene=SPTLC1&tab=mouse');
  for(const text of ['abnormal embryo size','embryonic growth retardation','preweaning lethality, complete penetrance','heterozygote','homozygote','E9.5','OMIM:620285','ORPHA:300605','Manually asserted; no statistical p-value'])assert.ok(html.includes(text),text);
  assert.doesNotMatch(html,/No summary is available/);
  documentEvents.change({target:{id:'mouse-zygosity',value:'heterozygote'}});
  const filtered=elements.get('#mouse-phenotype-results').innerHTML;
  assert.match(filtered,/embryonic growth retardation/);
  assert.doesNotMatch(filtered,/preweaning lethality/);
});
test('requested copy cleanup and current ClinVar provenance',async()=>{
  let html=await route('#gene=SPTLC1');
  assert.doesNotMatch(html,/Review the evidence together|No combined haplosufficiency probability|Total: distinct ClinVar alleles/);
  assert.equal(elements.get('#site-footer').textContent,'Test version');
  html=await route('#gene=SPTLC1&tab=variants');
  assert.match(html,/2026-11-09/);
  assert.match(html,/>Pre-PAM<\/th>/);
  assert.match(html,/>Post-PAM<\/th>/);
  assert.doesNotMatch(html,/HIGH quality|LOW only/);
  assert.match(html,/<td>WTD<\/td><td>119<\/td><td>93<\/td>/);
  assert.match(html,/SPTLC1-genome-targets.csv/);
});

test('full multiallelic run is prominent and separate from the smaller genome track',async()=>{
  for(const [gene,count] of [['RARB','2,427'],['RYR2','2,201']]){
    const html=await route(`#gene=${gene}&tab=variants`);
    assert.match(html,/Pipeline targeting candidates/);
    assert.ok(html.includes(`excision</td><td>${count}</td>`));
    assert.ok(html.indexOf('Pipeline targeting candidates')<html.indexOf('Genome-wide targeting candidates'));
    assert.doesNotMatch(html,/RUN_MULTIALLELIC/);
    assert.match(html,new RegExp(`downloads/run-multiallelic/${gene}.csv.gz`));
  }
});
