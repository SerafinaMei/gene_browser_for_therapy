import {present,fmt,esc,unpackIndex,searchGenes,flatten,csv,unpackRecordCsv,hiLabels,moiLabels,genomeUrl} from './core.js?v=6';
import {validateSettings,applySettings,pageTitle} from './settings.js?v=8';
const $=s=>document.querySelector(s);
let siteConfig;
let genes=[],manifest,current=null,activeTab='overview',requestId=0,results=[],highlight=-1;
let compare=[];
try { const saved=JSON.parse(localStorage.getItem('dnd-compare')||'[]');compare=Array.isArray(saved)?saved.filter(s=>typeof s==='string').slice(0,4):[]; } catch {}
const shards=new Map(), recordShards=new Map(), evidenceShards=new Map();
const tabs=[['overview','Summary'],['diseases','Clinical annotations'],['variants','Variants and targeting'],['mouse','Mouse phenotypes'],['raw','Raw data']];
const ext=(url,label)=>`<a href="${esc(url)}" target="_blank" rel="noopener">${label} ↗</a>`;
const tag=(text,type='')=>`<span class="tag ${type}">${esc(text)}</span>`;
const missing=()=>'<span class="quiet">Not available</span>';
const value=v=>present(v)?esc(fmt(v)):missing();
const panel=(title,subtitle,body,foot='',action='')=>`<section class="panel"><div class="panel-head"><div><h2>${title}</h2>${subtitle?`<p>${subtitle}</p>`:''}</div>${action}</div>${body}${foot?`<div class="panel-foot">${foot}</div>`:''}</section>`;
const stat=(label,v)=>`<div class="stat-row"><span>${label}</span><strong>${value(v)}</strong></div>`;
const more=(tab,label='View details →')=>`<button class="text-button" data-tab="${tab}">${label}</button>`;
const sourceLink=(id,label='Source and definitions')=>`<button class="text-button" data-source="${id}">${label}</button>`;
function notify(t){$('#toast').textContent=t;$('#toast').classList.add('show');clearTimeout(notify.timer);notify.timer=setTimeout(()=>$('#toast').classList.remove('show'),3200);}
function selectNav(id) {
  document.querySelectorAll('.site-nav button').forEach(b => {
    b.classList.toggle('selected', b.id === id);
    if (b.id === id) b.setAttribute('aria-current', 'page');
    else b.removeAttribute('aria-current');
  });
}
function setView(view) {
  document.body.dataset.view = view;
  selectNav(view === 'home' || view === 'report' ? 'search-nav' : view === 'browse' ? 'browse-nav' : 'sources-nav');
}
function showHome() {
  setView('home');
  current = null;
  document.title = pageTitle(siteConfig);
  $('#gene-search').value = '';
  $('#main').innerHTML = `<p class="home-status">${fmt(genes.length)} genes available. Data coverage varies by gene.</p>`;
}
function resourceUrl(path){const version=document.documentElement?.dataset.buildVersion;return version?`${path}?v=${encodeURIComponent(version)}`:path;}
async function fetchJSON(path){
  const compressed=document.documentElement?.dataset.jsonEncoding==='gzip';
  const r=await fetch(resourceUrl(path+(compressed?'.gz':'')));
  if(!r.ok)throw new Error(`${path}: HTTP ${r.status}`);
  if(!compressed)return r.json();
  return new Response((await r.blob()).stream().pipeThrough(new DecompressionStream('gzip'))).json();
}
async function loadGene(symbol) {
  const key = /^[A-Za-z]/.test(symbol) ? symbol[0] : '_';
  if (!shards.has(key)) shards.set(key, fetchJSON(`data/genes-${key}.json`).catch(e => { shards.delete(key); throw e; }));
  if (!recordShards.has(key)) recordShards.set(key, fetchJSON(`data/record-details-${key}.json`).catch(e => { recordShards.delete(key); throw e; }));
  if (!evidenceShards.has(key)) evidenceShards.set(key, fetchJSON(`data/evidence-${key}.json`).catch(e => { evidenceShards.delete(key); throw e; }));
  const [data, records, evidence] = await Promise.all([shards.get(key), recordShards.get(key), evidenceShards.get(key)]);
  if (!data[symbol]) throw new Error('Gene record unavailable');
  return {...data[symbol], ...records[symbol], ...evidence[symbol]};
}
function navigate(symbol,tab='overview'){location.hash=`gene=${encodeURIComponent(symbol)}${tab==='overview'?'':`&tab=${tab}`}`;}
function download(name,text,type){const blob=new Blob([text],{type});const url=URL.createObjectURL(blob);const a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
function updateCompare(){ $('#compare-count').textContent=compare.length;try{localStorage.setItem('dnd-compare',JSON.stringify(compare));}catch{} }
function toggleCompare(symbol){if(compare.includes(symbol))compare=compare.filter(s=>s!==symbol);else if(compare.length<4)compare.push(symbol);else return notify('Compare up to four genes. Remove one in Compare genes first.');updateCompare();notify(compare.includes(symbol)?`${symbol} added to comparison`:`${symbol} removed from comparison`);if(current?.symbol===symbol){const b=$('#add-compare');if(b)b.textContent=compare.includes(symbol)?'Remove from comparison':'Add to comparison';}}
function diseaseTable(rows, limit = Infinity) {
  if (!rows.length) return '<div class="empty">No gene-disease record in the ClinGen data used here.</div>';
  return `<div class="table-scroll"><table class="disease-table"><thead><tr><th>Disease</th><th>Inheritance</th><th>ClinGen classification</th><th>Last evaluated</th></tr></thead><tbody>${rows.slice(0,limit).map(r => `<tr>
    <td>${ext(r['ONLINE REPORT'],esc(r['DISEASE LABEL']))}<span class="secondary">${esc(r['DISEASE ID (MONDO)'])}</span></td>
    <td>${esc(moiLabels[r.MOI] || r.MOI)}${moiLabels[r.MOI] ? ` (${esc(r.MOI)})` : ''}</td>
    <td>${esc(r.CLASSIFICATION)}<details><summary>Review details</summary><p>${esc(r.GCEP)}</p></details></td>
    <td>${esc((r['CLASSIFICATION DATE'] || '').slice(0,10))}</td>
  </tr>`).join('')}</tbody></table></div>`;
}
function score(title, n, note) {
  return `<div class="score-block"><div class="score-top"><span>${title}</span><strong>${value(n)}</strong></div><p>${note}</p></div>`;
}
function pmp(g){if(g.symbol!=='PMP22')return '';return `<div class="note amber"><strong>PMP22-associated neuropathies</strong><br>Keep gene-level loss-of-function evidence separate from the recurrent 17p12 duplication associated with CMT1A. The ClinGen gene record and the regional CNV record have different scopes. ${ext('https://search.clinicalgenome.org/kb/gene-dosage/HGNC:9118','ClinGen explanation')}</div>`;}
function overallNote(g){if(g.hi_score===3)return '<strong>Curated evidence supports haploinsufficiency.</strong> Read this alongside the associated disease and mechanism. Computational predictions may differ.';if(g.hi_score===40)return '<strong>ClinGen reports dosage sensitivity unlikely.</strong> This is curated evidence about dosage sensitivity within the source’s scope.';if(g.hi_score===0)return '<strong>ClinGen reports no evidence for haploinsufficiency.</strong> This does not establish that losing one gene copy is tolerated.';return '';}
function overview(g) {
  const d = g.dnd || {}, m = g.mouse, c = g.collins || {}, sh = g.shet;
  const hi = present(g.hi_score) ? hiLabels[g.hi_score] : null;
  const measures = [
    ['D&D', g.dnd ? 'Yes' : 'No', `Membership in the supplied dominant and dispensable gene set: disease alleles proposed as candidates for inactivation while retaining a functional copy. ${ext('https://pmc.ncbi.nlm.nih.gov/articles/PMC13042096/', 'Ramey et al.')}`],
    ['ClinGen haploinsufficiency evidence', hi, 'Curated evidence that loss of one gene copy can cause disease.'],
    ['s_het (GeneBayes)', sh?.post_mean, 'Selection against heterozygous loss of function. Higher values indicate stronger constraint.'],
    ['s_het 95% credible interval', sh ? `${fmt(sh.post_lower_95,4)} - ${fmt(sh.post_upper_95,4)}` : null, 'Uncertainty in the GeneBayes estimate.'],
    ['pHaplo (Collins)', c.pHaplo, 'Probability of haploinsufficiency, from 0 to 1. Higher values indicate greater predicted HI.'],
    ['ClinGen %HI (DECIPHER)', g.hi_percentile, 'Percentile rank, not a probability. Lower values indicate greater predicted HI.']
  ];
  const table = `<div class="table-scroll"><table class="evidence-table"><thead><tr><th>Measure</th><th>Result</th><th>Meaning</th></tr></thead><tbody>${measures.map(([label,n,definition]) => `<tr><td>${label}</td><td>${present(n) ? esc(fmt(n,4)) : missing()}</td><td>${definition}</td></tr>`).join('')}</tbody></table></div>`;
  return `${panel('Dosage and constraint', '', table, overallNote(g), more('diseases', 'Full annotations'))}
    ${panel('ClinGen disease associations', 'ClinGen data dated 2025-11-13.', diseaseTable(g.diseases,3), g.diseases.length > 3 ? `Showing 3 of ${g.diseases.length} records. Open Clinical annotations to see all.` : '', more('diseases', 'All disease records'))}
    ${pmp(g)}
    <div class="two-col">
      ${panel('Variant counts', '', `<div class="panel-body">${stat('Total ClinVar variants', g.clinvar?.Total_alleles)}${stat('Dominant pathogenic / likely pathogenic variants', d.dominant_plp_mutations)}</div>`, '', more('variants', 'Count definitions'))}
      ${strategySummary(g)}
    </div>
    ${g.mouse_current ? currentMouseSummary(g) : panel('Mouse phenotype results', 'Results for mapped mouse homologs from IMPC / MGI.', m ? `<div class="panel-body">${stat('Heterozygous test results (successful analyses)', m.het_rows)}${stat('Statistically significant results', m.het_significant_rows)}${stat('Distinct parameters tested', m.het_parameters_tested)}</div>` : '<div class="empty">No mouse summary is available for this gene in the included datasets.</div>', 'Counts refer to test results, not individual animals.', more('mouse', 'Mouse details'))}`;
}
function dosage(g){return `${pmp(g)}${panel('ClinGen gene-disease validity','Each disease has its own inheritance pattern and evidence classification.',diseaseTable(g.diseases),'AD: autosomal dominant · AR: autosomal recessive · SD: semidominant · XL: X-linked')}${panel('Dosage evidence','Evidence for sensitivity to loss or gain of a gene copy.',g.dosage.length?`<div class="table-scroll"><table><thead><tr><th>Evidence type</th><th>Evidence</th><th>Associated disease</th><th>Evaluated</th></tr></thead><tbody>${g.dosage.map(r=>`<tr><td>Haploinsufficiency</td><td>${esc(hiLabels[g.hi_score] || r['HI Score'])}</td><td>${value(r['Haplo Disease'])}<span class="secondary">${esc(r['Haplo Disease ID']||'')}</span></td><td>${esc(r['LastEvaluated Date'])}</td></tr><tr><td>Triplosensitivity</td><td>${esc(({NoEvidence:'No evidence',LittleEvidence:'Little evidence',SomeEvidence:'Some evidence',SufficientEvidence:'Sufficient evidence',SensitivityUnlikely:'Dosage sensitivity unlikely'})[r['TS Score']] || r['TS Score'])}</td><td>${value(r['Triplo Disease'])}<span class="secondary">${esc(r['Triplo Disease ID']||'')}</span></td><td>${esc(r['LastEvaluated Date'])}</td></tr>`).join('')}</tbody></table></div>`:'<div class="empty">No gene-level dosage record in this snapshot.</div>',sourceLink('dosage'))}<div class="two-col">${panel('Evolutionary constraint · s_het','GeneBayes posterior estimate',`<div class="panel-body">${stat('Posterior mean',g.shet?.post_mean)}${stat('95% credible interval',g.shet?`${fmt(g.shet.post_lower_95)} - ${fmt(g.shet.post_upper_95)}`:null)}${stat('Prior mean',g.shet?.prior_mean)}${stat('Observed LoF variants',g.shet?.obs_lof)}${stat('Expected LoF variants',g.shet?.exp_lof)}<p class="hint">Selection against heterozygous loss of function. Higher values indicate stronger evolutionary constraint. This is not a disease probability.</p></div>`,sourceLink('shet'))}${panel('Computational dosage predictions','Prediction and curation are separate evidence types.',`<div class="panel-body">${score('Collins pHaplo',g.collins?.pHaplo,'Predicted probability of haploinsufficiency.')}${score('Collins pTriplo',g.collins?.pTriplo,'Predicted probability of triplosensitivity.')}${g.collins_records?.length>1?'<div class="note amber">Multiple historical gene records map here. Headline probabilities are withheld; inspect all source values in Raw data.</div>':''}${score('ClinGen %HI',g.hi_percentile,'Percentile rank (not a probability); lower indicates greater predicted HI.')}</div>`,sourceLink('collins'))}</div>${g.pmp22_regions?.length?panel('PMP22 recurrent-region records','Regional evidence is not merged with the single-gene dosage score.',`<div class="panel-body">${g.pmp22_regions.map(r=>`<div class="source-card"><h3>${esc(r['Gene Symbol /Region Name'])}</h3><p>HI: ${esc(r['HI Score'])} · TS: ${esc(({NoEvidence:'No evidence',LittleEvidence:'Little evidence',SomeEvidence:'Some evidence',SufficientEvidence:'Sufficient evidence',SensitivityUnlikely:'Dosage sensitivity unlikely'})[r['TS Score']] || r['TS Score'])}</p><p>${esc(r['Haplo Disease'])} / ${esc(r['Triplo Disease'])}</p><code>GRCh38 ${esc(r.GRCh38)}</code></div>`).join('')}</div>`):''}${panel('OMIM & disease context','Gene catalog links and associations in the ClinGen snapshot',`<div class="panel-body"><div class="link-row">${g.omim_ids.length?g.omim_ids.map(id=>ext(`https://omim.org/entry/${id}`,`OMIM ${esc(id)}`)).join(''):ext(`https://omim.org/search?search=${encodeURIComponent(g.symbol)}`,'Search OMIM')}</div><p class="hint">OMIM numbers are gene-entry identifiers from HGNC. Disease descriptions and inheritance above are from ClinGen.</p>${stat('CMT-labelled ClinGen assertions',g.cmt_assertions.length)}${stat('PMP22 disease group',g.symbol==='PMP22'?'PMP22-associated neuropathies':g.cmt_assertions.length?'CMT association at a different gene':'No PMP22-specific assignment')}</div>`)}`;}
function downloadButton(record, label, name) {
  return record ? `<button class="download-button" data-record-file="${esc(record.path)}" data-record-format="${esc(record.format || 'csv-gzip')}" data-record-key="${esc(record.record_key || '')}" data-filename="${esc(name)}">${label} (CSV)</button>` : '';
}

async function downloadRecords(button) {
  const label = button.textContent;
  button.disabled = true;
  button.textContent = 'Preparing CSV...';
  try {
    const response = await fetch(resourceUrl(button.dataset.recordFile));
    if (!response.ok) throw new Error('Download unavailable');
    const compressed = await response.blob();
    if (typeof DecompressionStream === 'undefined') {
      const packed = button.dataset.recordFormat === 'columnar-gzip';
      download(packed ? button.dataset.filename.replace(/\.csv$/, '.json.gz') : button.dataset.filename + '.gz', compressed, 'application/gzip');
      if (packed) notify('Downloaded compressed source data. Use a current browser for CSV conversion.');
    } else {
      const stream = compressed.stream().pipeThrough(new DecompressionStream('gzip'));
      const data = new Response(stream);
      const file = button.dataset.recordFormat === 'columnar-gzip'
        ? new Blob([unpackRecordCsv(await data.json(), button.dataset.recordKey)], {type:'text/csv;charset=utf-8'})
        : await data.blob();
      download(button.dataset.filename, file, 'text/csv;charset=utf-8');
    }
  } catch (error) {
    notify('Download failed. Check the local server and try again.');
  } finally {
    button.disabled = false;
    button.textContent = label;
  }
}

function cellLineTable(scope) {
  return `<div class="table-scroll"><table class="cell-line-table"><thead><tr><th>Dataset / cell line</th><th>SNPs</th><th>With NGG</th><th>Selected SNPs</th><th>Selected gRNAs</th><th>Editors</th></tr></thead><tbody>${scope.rows.map(r => `<tr><td>${esc(r.cell_line)}</td><td>${fmt(r.variants)}</td><td>${fmt(r.ngg_variants)}</td><td>${fmt(r.selected_variants)}</td><td>${fmt(r.guides)}</td><td>${esc(r.editors.join(', ') || '-')}</td></tr>`).join('')}</tbody></table></div>`;
}

function cellLinePanel(g) {
  if(g.pipeline_targets && g.genome_targets?.source_id==='sptlc_tracks')return publishedCellLinePanel(g)+genomeCellPanel({...g,genome_targets:g.pipeline_targets});
  if(g.pipeline_targets)return genomeCellPanel({...g,genome_targets:g.pipeline_targets})+`<details class="additional-counts"><summary>Published-track overlaps and earlier selected-guide datasets</summary>${publishedCellLinePanel(g)}</details>`;
  return publishedCellLinePanel(g);
}
function publishedCellLinePanel(g) {
  if (g.genome_targets) return genomeCellPanel(g) + (g.targeting || g.expanded_targeting ? `<details class="additional-counts"><summary>Earlier selected-guide dataset</summary>${legacyCellLinePanel(g)}</details>` : '');
  return legacyCellLinePanel(g);
}
function legacyCellLinePanel(g) {
  const target = g.targeting || g.expanded_targeting;
  const expanded = target?.dataset === 'expanded';
  if (!target) return panel('Cell-line targeting', '', '<div class="empty">No variant or gRNA records are available for this gene in the cell-line dataset.</div>');
  const controls = `<div class="panel-controls"><label for="target-window">Splice distance</label><select id="target-window">${target.scopes.map(s => `<option value="${esc(s.window)}">${esc(s.window)} bp + exonic SNPs</option>`).join('')}</select></div>`;
  const downloads = `<div class="download-row">${downloadButton(target.exports.variants, 'Targeting SNPs', `${g.symbol}-targeting-variants.csv`)}${downloadButton(target.exports.guides, 'gRNA sequences', `${g.symbol}-grna-sequences.csv`)}</div>`;
  return panel('Cell-line targeting', expanded ? 'Expanded indel/splice screen.' : 'Cell-line candidate dataset.', `${controls}<div id="cell-line-results">${cellLineTable(target.scopes[0])}</div>
    <div class="panel-body"><p class="hint">Cell-line rows count observed heterozygous SNPs. p47 confirms positions only; its alleles are unknown. Counts overlap across cell lines.</p>
    <details class="method-details"><summary>Definitions and selection criteria</summary>
      <p>SNP: single-nucleotide polymorphism. NGG: the SpCas9 PAM sequence. ABE: adenine base editor. CBE: cytosine base editor. Indel: exon disruption by insertion or deletion.</p>
      <p>Splice distance is the SNP-to-splice-boundary distance, not the editor activity window. Both settings also include exonic SNPs. Counts use distinct variant IDs and distinct genomic guides.</p>
      <p>${expanded ? 'Selected guides also pass the 10 bp SNP-to-PAM rule and the source coding/NMD or splice-coverage criteria. ' : ''}Selected candidates meet the source expression and genotype criteria. WTC and WTD have allele-matched VCF calls; WTD calls retain an unspecified FILTER value. p47 is locus-only evidence. Guides use GRCh38 plus the focal allele; neighboring variants are unphased.</p>
      <p>Coverage is limited to the available candidate export. These counts describe guide candidates, not experimentally validated activity.</p>
    </details></div>`, `${downloads}<p class="hint">Downloads include both distance settings and selection flags. Rows preserve variant/editor/window associations, so a SNP or gRNA may appear more than once.</p>${sourceLink(expanded ? 'expanded_variants' : 'targeting_variants')}`);
}

const strategyFields = [
  ['indel_pre', 'Exonic indel, before NMD assessment', 'num_indel_vars_pre_NMD_assessment', 'indel_targetable_pre_NMD_assessment', 'indel_pam_targetable'],
  ['indel_nmd', 'Exonic indel, predicted NMD', 'num_indel_amenable_vars_inducing_NMD', 'indel_targetable_post_NMD_assessment', 'indel_pam_targetable'],
  ['indel_escape', 'Exonic indel, predicted NMD escape', 'num_indel_amenable_vars_escaping_NMD', null, null],
  ['splice', 'Splice base editing', 'num_base_editable_vars', 'base_editable', 'base_edit_pam_targetable'],
  ['excision', 'Exon / gene excision', 'num_excision_vars', 'excision_targetable', 'excision_pam_targetable'],
  ['crisproff', 'CRISPRoff', 'num_crisproff_vars', 'crisproff_targetable', 'crisproff_pam_targetable'],
  ['excision_pairs', 'Exon / gene excision pairs', 'num_excision_pairs', null, null]
];

function strategyRows(g, dataset) {
  if (dataset === 'expanded') return g.strategy_screen.rows;
  return strategyFields.map(([key,label,count,flag,pam]) => ({key,label,count:g.dnd[count],
    targetable:flag ? g.dnd[flag] : null,pam:pam ? g.dnd[pam] : null,status:'reference'}));
}

function strategyContent(g, dataset) {
  const expanded = dataset === 'expanded';
  if (expanded && !g.strategy_screen.assessed) return `<div class="empty">Not assessed. ${esc(g.strategy_screen.exclusion)}</div>`;
  const bool = v => present(v) ? v ? 'Yes' : 'No' : '-';
  const table = `<div class="table-scroll"><table><thead><tr><th>Strategy</th><th>Candidate variants</th><th>Targetable</th><th>Viable PAM</th></tr></thead><tbody>${strategyRows(g,dataset).map(r => `<tr><td>${esc(r.label)}</td>${['not_run','not_assessed'].includes(r.status) ? `<td colspan="3">${r.status==='not_run' ? 'Not run in this dataset' : 'Not assessed: coding model unavailable'}</td>` : `<td>${value(r.count)}</td><td>${bool(r.targetable)}</td><td>${bool(r.pam)}</td>`}</tr>`).join('')}</tbody></table></div>`;
  const provenance = expanded ? 'Common SNPs (AF 0.1-0.9), before cell-line selection. Splice distance: 1-23 bp. Targetable: candidate present. Viable PAM: NGG within 10 bp of a candidate SNP.' : 'D&D reference dataset. Counts retain the source definitions for each strategy. Targetable and viable-PAM flags apply to the gene.';
  const downloadRow = expanded ? `<div class="download-row">${downloadButton(g.strategy_screen.exports.variants,'Screen SNPs',`${g.symbol}-expanded-screen-variants.csv`)}${downloadButton(g.strategy_screen.exports.guides,'Screen gRNAs',`${g.symbol}-expanded-screen-grnas.csv`)}</div>` : '';
  return table + `<div class="panel-foot"><p>${provenance}</p><p class="hint">0 means assessed with no candidates. NMD: nonsense-mediated decay. PAM: protospacer-adjacent motif. Counts may overlap; NMD induction and escape can leave mixed predictions unclassified.</p>${downloadRow}${sourceLink(expanded ? 'expanded_screen' : 'dnd')}</div>`;
}

function strategyPanel(g) {
  if(g.pipeline_targets && g.genome_targets?.source_id==='sptlc_tracks')return publishedStrategyPanel(g)+genomeTargetPanel({...g,genome_targets:g.pipeline_targets});
  if(g.pipeline_targets)return genomeTargetPanel({...g,genome_targets:g.pipeline_targets})+`<details class="additional-counts"><summary>Published tracks and earlier pipeline results</summary>${publishedStrategyPanel(g)}</details>`;
  return publishedStrategyPanel(g);
}
function publishedStrategyPanel(g) {
  if (g.genome_targets) return genomeTargetPanel(g) + (g.other_genome_targets ? `<details class="additional-counts"><summary>Separate all-protein-coding track</summary>${genomeTargetPanel({...g,genome_targets:g.other_genome_targets})}${genomeCellPanel({...g,genome_targets:g.other_genome_targets})}</details>` : '') + (g.dnd || g.strategy_screen ? `<details class="additional-counts"><summary>Earlier pipeline strategy results</summary>${legacyStrategyPanel(g)}</details>` : '');
  return legacyStrategyPanel(g);
}
function legacyStrategyPanel(g) {
  if (!g.dnd && !g.strategy_screen) return panel('Targeting strategies', '', `<div class="empty">${esc(g.symbol)} was not included in the available targeting runs. Targetability has not been assessed here.</div>`);
  const dataset = g.dnd ? 'reference' : 'expanded';
  const controls = g.dnd && g.strategy_screen ? '<div class="panel-controls"><label for="strategy-dataset">Dataset</label><select id="strategy-dataset"><option value="reference">D&D reference</option><option value="expanded">Expanded indel/splice screen</option></select></div>' : '';
  return panel('Targeting strategies', '', `${controls}<div id="strategy-results">${strategyContent(g,dataset)}</div>`);
}

function strategySummary(g) {
  if(g.pipeline_targets && g.genome_targets?.source_id!=='sptlc_tracks')g={...g,genome_targets:g.pipeline_targets};
  if (g.genome_targets) return panel('Targeting evidence',esc(g.genome_targets.dataset_title),`<div class="panel-body">${stat('Distinct candidate variants',g.genome_targets.variants)}${Object.entries(g.genome_targets.strategies).map(([s,n])=>stat(s,n)).join('')}</div>`,'',more('variants','Candidates and cell lines'));
  if (!g.dnd && (!g.strategy_screen || !g.strategy_screen.assessed)) return panel('Targeting evidence', '', `<div class="empty">${g.strategy_screen ? esc(g.strategy_screen.exclusion) : 'This gene was not included in the available targeting runs.'}</div>`);
  const dataset = g.dnd ? 'reference' : 'expanded';
  const rows = strategyRows(g,dataset).filter(r=>['indel_nmd','splice','excision','crisproff'].includes(r.key) && r.status!=='not_run');
  return panel('Targeting evidence', dataset==='reference' ? 'D&D reference dataset.' : 'Expanded indel/splice screen.', `<div class="panel-body">${rows.map(r=>stat(esc(r.label),r.count)).join('')}</div>`, '', more('variants','All strategies'));
}

function variants(g) {
  const d = g.dnd || {}, c = g.clinvar || {};
  const counts = [
    ['Total ClinVar variants', c.Total_alleles, 'Distinct alleles in or overlapping the gene', g.clinvar_date || '2026-10-24'],
    ['ClinVar pathogenic / likely pathogenic variants', c.Alleles_reported_Pathogenic_Likely_pathogenic, 'Excludes overlapping structural variants', g.clinvar_date || '2026-10-24'],
    ['Dominant pathogenic / likely pathogenic variants', d.dominant_plp_mutations, 'Dominant subset in the supplied D&D dataset', '2025-12-23']
  ];
  const extraCounts = [
    ['ClinVar submissions', c.Total_submissions, 'Multiple submissions can describe one variant', g.clinvar_date || '2026-10-24'],
    ['Dominant missense P/LP variants', d.dominant_mis_plp_mutations, 'D&D dataset', '2025-12-23'],
    ['Dominant nonsense P/LP variants', d.dominant_ns_plp_mutations, 'D&D dataset', '2025-12-23'],
    ['Total P/LP variants in the D&D dataset', d.total_plp, 'All inheritance modes', '2025-12-23'],
    ['ClinVar-OMIM dominant count', g.legacy_dominant_count, 'Separate source; inclusion rules unspecified', 'Unspecified']
  ];
  const countTable = rows => `<div class="table-scroll"><table class="variant-count-table"><thead><tr><th>Measure</th><th>Count</th><th>Definition</th><th>Source date</th></tr></thead><tbody>${rows.map(([label,n,definition,date]) => `<tr><td>${label}</td><td>${value(n)}</td><td>${definition}</td><td>${date}</td></tr>`).join('')}</tbody></table></div>`;
  const exportInfo = g.clinvar_export;
  const exportNote = exportInfo ? `<p class="hint">Download: ${fmt(exportInfo.alleles)} GRCh38 alleles from a separate local ClinVar file with an unspecified snapshot date. This may differ from the totals above. Dominant-only variant assignments are unavailable.</p>` : '<p class="hint">No individual ClinVar records are available in the local export for this gene.</p>';
  return panel('ClinVar variants', '', countTable(counts) + `<details class="additional-counts"><summary>Additional counts</summary>${countTable(extraCounts)}</details>`,
    `<p>P/LP: pathogenic or likely pathogenic. Counts use different source dates; no dominant fraction is calculated.</p>${exportNote}<div class="download-row">${downloadButton(exportInfo, 'ClinVar variants', `${g.symbol}-clinvar-GRCh38.csv`)}${sourceLink(g.clinvar_date?'clinvar_current':'clinvar', 'Count source')}${exportInfo ? sourceLink('clinvar_variants', 'Download source') : ''}</div>`)
    + strategyPanel(g)
    + cellLinePanel(g);
}

function legacyMouseView(g){const m=g.mouse;if(!m)return panel('Mouse phenotype evidence','IMPC / MGI · mapped human-mouse homologs','<div class="empty">No summary is available for this gene in the supplied cohorts. This is not a zero-result or no-phenotype finding.</div>');const total=m.het_rows||0;return `<div class="note"><strong>Mouse model evidence</strong> describes mapped mouse homologs. Non-significant test results do not prove biological tolerance, and absence of a homozygous assertion does not establish that homozygotes were untested.</div>${panel('Heterozygous results','Unique test results from successful IMPC analyses, combined across the included gene sets.',`<div class="panel-body">${stat('Human gene',g.symbol)}${stat('Mapped mouse symbols with het results',m.mouse_symbols.join(', ')||null)}${stat('All heterozygous test results (successful analyses)',m.het_rows)}${stat('Statistically significant heterozygous results',m.het_significant_rows)}${stat('Non-significant heterozygous results',m.het_nonsignificant_rows)}${stat('Results with significance not recorded',m.het_significance_unknown_rows)}${total?`<div class="stack" role="img" aria-label="${m.het_significant_rows} significant and ${m.het_nonsignificant_rows} non-significant test rows"><span class="sig" style="width:${m.het_significant_rows/total*100}%"></span><span class="nonsig" style="width:${m.het_nonsignificant_rows/total*100}%"></span></div><div class="legend"><span><i style="background:#4c5b6c"></i>Significant</span><span><i style="background:#c0c8d1"></i>Non-significant</span></div>`:''}</div>`,'The source significant flag is retained. No new p-value threshold is applied.')}<div class="two-col">${panel('Testing coverage','The number of different measurements and experimental procedures.',`<div class="panel-body">${stat('Distinct parameters tested',m.het_parameters_tested)}${stat('Distinct significant parameters',m.het_parameters_significant)}${stat('Distinct procedures tested',m.het_procedures_tested)}${stat('Homozygote assertion in source',m.has_hom_record?'Present':'Not found in supplied assertions')}<p class="hint">A homozygous record alone does not establish normal heterozygous phenotypes or breeding ability.</p></div>`)}${panel('Mouse result summary','Summary of the included heterozygous test results.',`<div class="panel-body"><p>${esc(m.het_class)}</p><p class="hint">This summary includes all measured phenotypes. The project’s mortality and neurobehavioral categories use a narrower set of phenotypes.</p><code>het_class</code></div>`)}</div>${panel('Mouse data fields','Field names preserved for analysis and export.',rawTable(Object.entries(m)),m.sources.map(id=>sourceLink(id)).join(' · '))}`;}
function genomeTargetPanel(g) {
  const t=g.genome_targets;
  if(!t)return '';
  return panel(t.source_id==='run_multiallelic'?'Pipeline targeting candidates':t.source_id==='sptlc_tracks'?'SPTLC targeting candidates':'Genome-wide targeting candidates',`${esc(t.dataset_title)} · GRCh38`,`<div class="panel-body">${stat('Distinct candidate variants',t.variants)}${stat('Variants marked PAM-targetable in source',t.pam_variants)}${present(t.excision_pairs)?stat('Excision pairs (separate from variants)',t.excision_pairs):''}</div><div class="table-scroll"><table><thead><tr><th>Editing strategy</th><th>Distinct variants</th></tr></thead><tbody>${Object.entries(t.strategies).map(([s,n])=>`<tr><td>${esc(s)}</td><td>${fmt(n)}</td></tr>`).join('')}</tbody></table></div>`, `<p>${t.source_id==='run_multiallelic'?'All per-gene variant records reconcile to the RUN_MULTIALLELIC summary. Exon-disruption candidates use the source NMD assessment. ':''}Strategies overlap. These are common-variant editing candidates, not ClinVar pathogenic variants or selected gRNAs. Excision partners and source PAM flags are included in the download.</p><div class="download-row">${downloadButton(t.export,'All candidate variants and cell-line matches',`${g.symbol}-genome-targets.csv`)}${sourceLink(t.source_id)}</div>`);
}
function genomeCellPanel(g) {
  const t=g.genome_targets;
  if(!t)return '';
  return panel('Cell-line heterozygous candidates',`Allele-matched overlap · ${esc(t.dataset_title)}`,`<div class="table-scroll"><table><thead><tr><th>Cell line</th><th>All matched variants</th><th>HIGH quality</th><th>LOW only</th><th>HIGH with source PAM flag</th></tr></thead><tbody>${t.cells.map(r=>`<tr><td>${esc(r.cell_line)}</td><td>${fmt(r.variants)}</td><td>${fmt(r.high_quality)}</td><td>${fmt(r.low_only)}</td><td>${fmt(r.pam_high_quality)}</td></tr>`).join('')}</tbody></table></div>`,`<p>${t.source_id==='run_multiallelic'?'RUN_MULTIALLELIC candidates intersected with the current published cell-line tracks. ':''}Matches require the same GRCh38 interval and allele pair. HIGH/LOW are the original track quality labels. A zero means no matching candidate in these tracks. Expression, phasing and selected-guide criteria have not been applied to these counts.</p><details><summary>Cell-line sources</summary>${t.cells.map(r=>sourceLink('track_'+r.cell_line,r.cell_line)).join(' · ')}</details>`);
}
function mousePhenotypeTable(m,zygosity='all') {
  const rows=m.phenotypes.filter(r=>zygosity==='all'||r.zygosity===zygosity);
  if(!rows.length)return '<div class="empty">No significant phenotype assertions in this snapshot for the selected zygosity.</div>';
  const list=xs=>[...new Set(xs.flat().filter(Boolean))].join(', ');
  return `<div class="table-scroll"><table class="mouse-table"><thead><tr><th>Phenotype</th><th>System</th><th>Allele / zygosity</th><th>Life stage</th><th>Sex</th><th>Evidence</th></tr></thead><tbody>${rows.map(r=>`<tr><td>${esc(r.term)}<span class="secondary">${esc(r.term_id)}</span></td><td>${esc(list(r.support.map(x=>x.top_level_mp_term_name)))}</td><td>${esc(r.allele)}<span class="secondary">${esc(r.zygosity)}</span></td><td>${esc(r.life_stage||'Unspecified')}</td><td>${esc(list(r.support.map(x=>x.sex==='not_considered'?'Not assessed by sex':x.sex))||'Unspecified')}</td><td><details><summary>${r.support.length} supporting dataset${r.support.length===1?'':'s'}</summary>${r.support.map(x=>`<p><strong>${esc(x.procedure_name)}</strong><br>${esc(x.parameter_name)}<br>${esc(x.phenotyping_center)} · ${x.assertion_type==='manual'||x.statistical_method==='Supplied as data'?'Manually asserted; no statistical p-value':`p = ${esc(x.p_value??'not reported')}`}</p>`).join('')}</details></td></tr>`).join('')}</tbody></table></div>`;
}
function mouseView(g) {
  const m=g.mouse_current;
  if(!m)return legacyMouseView(g);
  const links=m.models.map(r=>ext(`https://www.mousephenotype.org/data/genes/${encodeURIComponent(r.mgi_accession_id)}`,esc(r.marker_symbol)+' · '+esc(r.mgi_accession_id))).join(' · ');
  const phenotypes=panel('Significant mouse phenotypes',`IMPC release 24.0 · ${m.unique_terms} distinct phenotype terms · ${m.supporting_datasets} supporting datasets`,`<div class="panel-body">${links||'Mapped IMPC / MGI evidence'}<p>${m.models.map(r=>`${esc(r.marker_symbol)}: ${esc(r.phenotype_status||'Testing status unavailable')}`).join(' · ')}</p></div><div class="panel-controls"><label for="mouse-zygosity">Zygosity</label><select id="mouse-zygosity"><option value="all">All</option>${[...new Set(m.phenotypes.map(r=>r.zygosity))].sort().map(z=>`<option value="${esc(z)}">${esc(z)}</option>`).join('')}</select><button id="export-mouse">Download phenotype evidence CSV</button></div><div id="mouse-phenotype-results">${mousePhenotypeTable(m)}</div>`,`${sourceLink('impc_genotype-phenotype','IMPC phenotype source')} · Grouped by mouse gene, phenotype, allele, zygosity and life stage. Manual assertions are included; statistical significance is not inferred from a placeholder p-value.`);
  const diseases=panel('Associated human diseases','IMPC / PhenoDigm gene–disease associations; phenotype similarity is separate from causation.',m.diseases.length?`<div class="table-scroll"><table><thead><tr><th>Disease</th><th>Identifier</th><th>Matching human phenotypes</th><th>Mouse evidence sources</th></tr></thead><tbody>${m.diseases.map(r=>`<tr><td>${esc(r.disease_term)}</td><td>${r.disease_id.startsWith('OMIM:')?ext('https://omim.org/entry/'+r.disease_id.split(':')[1],esc(r.disease_id)):r.disease_id.startsWith('ORPHA:')?ext('https://www.orpha.net/en/disease/detail/'+r.disease_id.split(':')[1],esc(r.disease_id)):esc(r.disease_id)}</td><td>${r.matched_phenotypes.length?r.matched_phenotypes.map(esc).join('<br>'):'Not supplied in this API snapshot'}</td><td>${esc(r.model_sources.join(', ')||'No model comparison recorded')}</td></tr>`).join('')}</tbody></table></div>`:'<div class="empty">No associated disease record in this IMPC snapshot.</div>',`${sourceLink('impc_diseases','Disease associations')} · ${sourceLink('impc_disease-models','Phenotype matches')}`);
  const systems=panel('Physiological systems','Coverage reported by IMPC for each mapped mouse gene.',`<div class="table-scroll"><table><thead><tr><th>Mouse gene</th><th>Significant phenotypes observed</th><th>Tested without a significant phenotype</th></tr></thead><tbody>${m.models.map(r=>`<tr><td>${esc(r.marker_symbol)}</td><td>${(r.significant_top_level_mp_terms||[]).map(esc).join('<br>')||'None recorded'}</td><td>${(r.not_significant_top_level_mp_terms||[]).map(esc).join('<br>')||'None recorded'}</td></tr>`).join('')}</tbody></table></div>`,sourceLink('impc_gene','IMPC gene coverage'));
  return phenotypes+diseases+systems+(g.mouse?`<details class="additional-counts"><summary>Earlier cohort analysis</summary>${legacyMouseView(g)}</details>`:'');
}
function currentMouseSummary(g) {
  const m=g.mouse_current;
  if(!m)return '';
  return panel('Mouse phenotype results','IMPC release 24.0',`<div class="panel-body">${stat('Distinct significant phenotype terms',m.unique_terms)}${stat('Associated disease records',m.diseases.length)}<p>${[...new Set(m.phenotypes.map(r=>r.term))].map(esc).join('; ')||'No significant assertion in this snapshot.'}</p></div>`,'',more('mouse','Phenotypes and diseases'));
}

function rawTable(entries){return `<div class="table-scroll"><table class="raw-table"><thead><tr><th>Field</th><th>Value</th></tr></thead><tbody>${entries.map(([k,v])=>`<tr><td>${esc(k)}</td><td>${v===null?'<span class="quiet">null · not available</span>':typeof v==='object'?esc(JSON.stringify(v)):esc(String(v))}</td></tr>`).join('')}</tbody></table></div>`;}
function rawView(g){return `<div class="panel-head"><div><h2>Raw gene data</h2><p>Original source fields, ClinGen records and calculated mouse counts.</p></div><button data-export="csv">Download CSV ↓</button></div><label class="sr-only" for="raw-search">Filter raw fields</label><input id="raw-search" class="raw-search" placeholder="Filter fields, e.g. het_rows, pHaplo, dominant…"><div class="panel" id="raw-table">${rawTable(Object.entries(flatten(g)))}</div>`;}
function header(g) {
  return `<div class="gene-heading"><div>
    <h1>${esc(g.symbol)}</h1><p class="gene-name">${esc(g.name)}</p>
    <div class="gene-ids"><span>${esc(g.hgnc_id)}</span>${g.ensembl_id ? `<span>${esc(g.ensembl_id)}</span>` : ''}<span>${esc(g.cytoband || 'Location unavailable')}</span></div>
    ${g.cmt_assertions.length || g.symbol === 'PMP22' ? `<div class="gene-flags">${g.cmt_assertions.length ? tag('CMT disease association') : ''}${g.symbol === 'PMP22' ? tag('PMP22-associated neuropathies') : ''}</div>` : ''}
    </div><div class="actions"><button id="add-compare" data-compare="${esc(g.symbol)}">${compare.includes(g.symbol) ? 'Remove from comparison' : 'Add to comparison'}</button><button data-export="json">Export report</button></div></div>
    <div class="external-links">${ext(`https://search.clinicalgenome.org/kb/genes/${g.hgnc_id}`, 'ClinGen')}${ext(`https://www.ncbi.nlm.nih.gov/clinvar/?term=${encodeURIComponent(g.entrez_id ? g.entrez_id+'[geneid]' : g.symbol+'[gene]')}`, 'ClinVar')}${ext(g.omim_ids.length ? `https://omim.org/entry/${g.omim_ids[0]}` : `https://omim.org/search?search=${encodeURIComponent(g.symbol)}`, 'OMIM')}${ext(genomeUrl(g), 'UCSC (GRCh38)')}${ext(genomeUrl(g, true), 'Genome Tracks')}</div>
    <div class="tabs" role="tablist" aria-label="Gene report sections">${tabs.map(([id,label]) => `<button id="tab-${id}" role="tab" aria-controls="tab-panel" aria-selected="${id===activeTab}" tabindex="${id===activeTab ? 0 : -1}" data-tab="${id}">${label}</button>`).join('')}</div>`;
}
function renderGene() {
  const g = current;
  setView('report');
  document.title = pageTitle(siteConfig,g.symbol);
  $('#gene-search').value = g.symbol;
  $('#main').innerHTML = header(g) + `<div id="tab-panel" role="tabpanel" aria-labelledby="tab-${activeTab}">${({overview,diseases:dosage,variants,mouse:mouseView,raw:rawView}[activeTab] || overview)(g)}</div>`;
}
function sourceCard(s){return `<div class="source-card" id="source-${esc(s.id)}"><h3>${s.url?ext(s.url,esc(s.title)):esc(s.title)}</h3><p>${esc(s.note)}</p><p><strong>Snapshot:</strong> ${esc(s.date||'Date not encoded; per-record evaluation dates retained where available.')}</p><details><summary>Source path and integrity hash</summary><code>${esc(s.path)}</code><code>SHA-256 ${esc(s.sha256)}</code></details></div>`;}
function showSources(){setView('sources');document.title=pageTitle(siteConfig,'Data sources');$('#main').innerHTML=`<p class="eyebrow">DATA PROVENANCE</p><h1 class="section-title">Data sources and methods</h1><p class="subheading">Source files, dates and definitions for each dataset are listed below.</p><div class="note">${manifest.notes.map(n=>esc(n)).join('<br>')}</div><div class="panel">${manifest.sources.map(sourceCard).join('')}</div>${panel('Join audit','Stable HGNC IDs are preferred; historical aliases are used only when unambiguous.',`<div class="panel-body">${Object.entries(manifest.audit.unmatched).map(([key,v])=>stat(`${esc(key)} · unmatched input records`,v.length)).join('')}${stat('Historical prediction collisions',manifest.audit.conflicts.length)}<p class="hint">Unmatched rows are excluded from gene reports and recorded in the manifest. Ambiguous Collins predictions are retained in Raw data with the headline value withheld.</p><button data-download-manifest>Download source manifest ↓</button></div>`)}`;}
let browse={q:'',group:'all',sort:'symbol',page:0};
function filteredGenes(){const q=browse.q.trim();let rows=q?searchGenes(genes,q,genes.length):[...genes];rows=rows.filter(g=>browse.group==='all'||browse.group==='dnd'&&g.dnd_member||browse.group==='cmt'&&g.cmt||browse.group==='hi'&&g.hi_score===3);const sort=browse.sort;rows.sort((a,b)=>sort==='symbol'?a.symbol.localeCompare(b.symbol):(present(b[sort])?b[sort]:-1)-(present(a[sort])?a[sort]:-1)||a.symbol.localeCompare(b.symbol));return rows;}
function showBrowse(){setView('browse');document.title=pageTitle(siteConfig,'Browse genes');$('#main').innerHTML=`<p class="eyebrow">EXPLORE THE DATASET</p><h1 class="section-title">Browse genes</h1><p class="subheading">Search ${fmt(genes.length)} genes or focus on a source cohort. Open a gene for its diseases, score definitions and raw data.</p><div class="filters"><input id="browse-search" aria-label="Filter genes" placeholder="Filter by gene name or ID" value="${esc(browse.q)}"><label>Gene set <select id="browse-group"><option value="all">All source genes</option><option value="dnd">D&D source set</option><option value="cmt">CMT-labelled assertions</option><option value="hi">ClinGen HI = 3</option></select></label><label>Sort by <select id="browse-sort"><option value="symbol">Symbol · A-Z</option><option value="phaplo">pHaplo · high to low</option><option value="shet">s_het · high to low</option><option value="clinvar_total">ClinVar total · high to low</option><option value="dominant_plp">Dominant P/LP · high to low</option></select></label><button id="export-browse">Export filtered CSV ↓</button></div><div class="panel" id="browse-results"></div><p class="hint">CMT filters match ClinGen disease labels, including limited or disputed assertions. The filter is not a pathogenicity conclusion. Blank metrics mean unavailable.</p>`;$('#browse-group').value=browse.group;$('#browse-sort').value=browse.sort;renderBrowseRows();}
function renderBrowseRows(){const rows=filteredGenes(),pages=Math.max(1,Math.ceil(rows.length/25));browse.page=Math.min(browse.page,pages-1);$('#browse-results').innerHTML=`<div class="table-scroll"><table class="browse-table"><thead><tr><th>Gene</th><th>ClinGen HI</th><th>s_het</th><th>pHaplo</th><th>ClinVar alleles</th><th>Dominant P/LP</th><th>Compare</th></tr></thead><tbody>${rows.slice(browse.page*25,(browse.page+1)*25).map(g=>`<tr><td><button class="gene-button" data-gene="${esc(g.symbol)}">${esc(g.symbol)} ↗</button><span class="secondary">${esc(g.name)}</span></td><td>${present(g.hi_score)?esc(hiLabels[g.hi_score]):'-'}</td><td>${present(g.shet)?fmt(g.shet):'-'}</td><td>${present(g.phaplo)?fmt(g.phaplo):'-'}</td><td>${present(g.clinvar_total)?fmt(g.clinvar_total):'-'}</td><td>${present(g.dominant_plp)?fmt(g.dominant_plp):'-'}</td><td><button data-compare="${esc(g.symbol)}" aria-label="Compare ${esc(g.symbol)}">${compare.includes(g.symbol)?'✓':'+'}</button></td></tr>`).join('')||'<tr><td colspan="7" class="empty">No genes match these filters.</td></tr>'}</tbody></table></div><div class="pagination"><span>${fmt(rows.length)} genes · page ${browse.page+1} of ${pages}</span><div><button data-page="-1" ${browse.page===0?'disabled':''}>← Previous</button> <button data-page="1" ${browse.page>=pages-1?'disabled':''}>Next →</button></div></div>`;}
function openDialog(html){$('#dialog-content').innerHTML=html;if(!$('#dialog').open)$('#dialog').showModal();}
async function showCompare(){openDialog('<h2>Compare genes</h2><p class="hint">Loading comparison…</p>');try{const records=await Promise.all(compare.map(loadGene));const rows=[['Gene name',g=>g.name],['ClinGen HI evidence',g=>present(g.hi_score)?hiLabels[g.hi_score]:null],['GeneBayes s_het',g=>g.shet?.post_mean],['95% credible interval',g=>g.shet?`${fmt(g.shet.post_lower_95)} - ${fmt(g.shet.post_upper_95)}`:null],['Collins pHaplo · P(HI)',g=>g.collins?.pHaplo],['ClinGen %HI · percentile',g=>g.hi_percentile],['ClinVar total alleles · 2026',g=>g.clinvar?.Total_alleles],['Dominant P/LP · older D&D',g=>g.dnd?.dominant_plp_mutations],['Het rows · mouse',g=>g.mouse?.het_rows],['Significant het rows · mouse',g=>g.mouse?.het_significant_rows],['CMT-labelled assertions',g=>g.cmt_assertions.length]];openDialog(`<p class="eyebrow">SIDE BY SIDE</p><h2 class="section-title">Compare genes</h2><p class="subheading">Up to four genes. Score directions and snapshot definitions remain separate.</p>${records.length?`<div class="table-scroll"><table class="comparison-table"><thead><tr><th>Measure</th>${records.map(g=>`<th>${esc(g.symbol)} <button class="text-button" data-remove-compare="${esc(g.symbol)}" aria-label="Remove ${esc(g.symbol)}">×</button></th>`).join('')}</tr></thead><tbody>${rows.map(([label,get])=>`<tr><td>${label}</td>${records.map(g=>`<td>${value(get(g))}</td>`).join('')}</tr>`).join('')}</tbody></table></div><button class="spacer-top" id="download-comparison">Download comparison CSV ↓</button>`:'<div class="empty">Use “Add to comparison” on a gene report, or the + button in Browse genes.</div>'}`);if($('#download-comparison'))$('#download-comparison').onclick=()=>download('gene-comparison.csv',csv([['Metric',...records.map(g=>g.symbol)],...rows.map(([label,get])=>[label,...records.map(get)])]),'text/csv');}catch(e){openDialog('<h2>Comparison could not load</h2><p class="hint">Check the connection and try again.</p>');}}
function help(){openDialog(`<p class="eyebrow">A QUICK GUIDE</p><h2 class="section-title">Using the gene browser</h2><div class="help-grid"><div><h3>Start with a symbol or stable ID</h3><p>Search a gene name, HGNC: ID, Ensembl ENSG ID, NCBI Gene ID, or recorded alias. Use Browse genes to filter the dataset and Compare to review up to four genes.</p></div><div><h3>Read the direction of each score</h3><p>Collins pHaplo is P(haploinsufficiency), so higher means greater predicted HI. ClinGen %HI is a percentile rank with the opposite direction. GeneBayes s_het is evolutionary selection, with a 95% credible interval. They are not interchangeable probabilities.</p></div><div><h3>Keep disease and inheritance together</h3><p>A gene can have several diseases and modes of inheritance. Each ClinGen assertion stays intact. “No evidence” for HI does not establish haplosufficiency.</p></div><div><h3>Check what each count measures</h3><p>ClinVar alleles, ClinVar submissions, dominant P/LP variants, targeting candidates and mouse result rows are different units. Each section identifies its snapshot and definition.</p></div><div><h3>Export data</h3><p>Download ClinVar records, targeting SNPs and gRNA sequences from Variants and targeting when available. Reports, raw fields and comparisons can also be exported. Source paths and SHA-256 hashes are available under Data sources.</p></div></div>`);}
async function route() {
  const id = ++requestId;
  hideResults();
  const p = new URLSearchParams(location.hash.slice(1));
  if (p.get('view') === 'browse') { showBrowse(); return; }
  if (p.get('view') === 'sources') { showSources(); return; }
  const q = (p.get('gene') || '').trim();
  if (!q) { showHome(); return; }
  setView('report');
  const matches = searchGenes(genes, q, 20);
  const normalized = q.toUpperCase().replace(/^(ENSG\d+)\.\d+$/, '$1');
  const exact = matches.find(g => g.symbol.toUpperCase() === normalized || [g.hgnc_id,g.ensembl_id,g.entrez_id].includes(normalized));
  const match = exact || (matches.length === 1 ? matches[0] : null);
  if (!match) {
    current = null;
    $('#gene-search').value = q;
    $('#main').innerHTML = `<div class="not-found"><h1>${matches.length ? 'Select a gene' : 'Gene not found'}</h1><p>${esc(q)} ${matches.length ? 'matches more than one gene.' : 'has no matching record in this dataset.'}</p>${matches.map(g => `<button data-gene="${esc(g.symbol)}">${esc(g.symbol)} (${esc(g.hgnc_id)})</button>`).join('')}<p class="hint">Enter an approved symbol, a gene name, or a stable gene ID above.</p><button id="browse-empty">Browse all genes</button></div>`;
    return;
  }
  activeTab = tabs.some(t => t[0] === p.get('tab')) ? p.get('tab') : 'overview';
  try {
    const g = await loadGene(match.symbol);
    if (id !== requestId) return;
    current = g;
    renderGene();
  } catch(e) {
    if (id !== requestId) return;
    current = null;
    $('#main').innerHTML = '<div class="not-found"><h1>Gene data could not load</h1><p>The data file could not be read. Check the local server and try again.</p><button id="retry">Try again</button></div>';
  }
}
function hideResults(){$('#search-results').hidden=true;$('#gene-search').setAttribute('aria-expanded','false');$('#gene-search').removeAttribute('aria-activedescendant');highlight=-1;}
function renderSearch(){const q=$('#gene-search').value;results=searchGenes(genes,q);highlight=-1;const box=$('#search-results');if(!q.trim())return hideResults();box.innerHTML=results.map((g,i)=>`<div class="search-option" id="result-${i}" role="option" aria-selected="false" data-result="${i}"><strong>${esc(g.symbol)} <span class="quiet">${esc(g.hgnc_id)}</span></strong><small>${esc(g.name)}</small></div>`).join('')||'<div class="empty">No matching gene in this snapshot.</div>';box.hidden=false;$('#gene-search').setAttribute('aria-expanded','true');}
document.addEventListener('click',async e=>{const b=e.target.closest('button,[data-result]');if(b){if(b.dataset.recordFile)await downloadRecords(b);if(b.dataset.gene){$('#dialog').close();navigate(b.dataset.gene);hideResults();}if(b.dataset.tab&&current){navigate(current.symbol,b.dataset.tab);}if(b.dataset.compare){toggleCompare(b.dataset.compare);if($('#browse-results'))renderBrowseRows();}if(b.dataset.removeCompare){compare=compare.filter(s=>s!==b.dataset.removeCompare);updateCompare();showCompare();if(current&&$('#add-compare'))$('#add-compare').textContent=compare.includes(current.symbol)?'Remove from comparison':'Add to comparison';}if(b.dataset.result!==undefined){navigate(results[Number(b.dataset.result)].symbol);hideResults();}if(b.dataset.source){const s=manifest.sources.find(s=>s.id===b.dataset.source);if(s)openDialog(sourceCard(s));}if(b.dataset.export&&current){const pack={gene:current,exported_at:new Date().toISOString(),source_manifest:manifest};download(`${current.symbol}-evidence.${b.dataset.export}`,b.dataset.export==='json'?JSON.stringify(pack,null,2):csv([['field','value'],...Object.entries(flatten(pack))]),b.dataset.export==='json'?'application/json':'text/csv');}if(b.hasAttribute('data-download-manifest'))download('source-manifest.json',JSON.stringify(manifest,null,2),'application/json');if(b.dataset.page){browse.page+=Number(b.dataset.page);renderBrowseRows();}if(b.id==='export-mouse'&&current?.mouse_current){const rows=current.mouse_current.phenotypes.flatMap(p=>p.support);const fields=[...new Set(rows.flatMap(Object.keys))];download(`${current.symbol}-mouse-phenotypes.csv`,csv([fields,...rows.map(r=>fields.map(k=>Array.isArray(r[k])?r[k].join('; '):r[k]))]),'text/csv');}if(b.id==='browse-empty')location.hash='view=browse';if(b.id==='retry')route();if(b.id==='export-browse'){const f=['symbol','name','hgnc_id','ensembl_id','hi_score','shet','phaplo','clinvar_total','dominant_plp'];download('filtered-genes.csv',csv([f,...filteredGenes().map(g=>f.map(k=>g[k]))]),'text/csv');}}if(!e.target.closest('#search-form'))hideResults();});
document.addEventListener('input',e=>{if(e.target.id==='raw-search'){const q=e.target.value.toLowerCase();$('#raw-table').innerHTML=rawTable(Object.entries(flatten(current)).filter(([k,v])=>`${k} ${v??''}`.toLowerCase().includes(q)));}if(e.target.id==='browse-search'){browse.q=e.target.value;browse.page=0;renderBrowseRows();}});
document.addEventListener('change',e=>{if(e.target.id==='mouse-zygosity'&&current?.mouse_current){$('#mouse-phenotype-results').innerHTML=mousePhenotypeTable(current.mouse_current,e.target.value);}if(e.target.id==='strategy-dataset'&&current){$('#strategy-results').innerHTML=strategyContent(current,e.target.value);}if(e.target.id==='target-window'&&(current?.targeting||current?.expanded_targeting)){const scope=(current.targeting||current.expanded_targeting).scopes.find(s=>s.window===e.target.value);if(scope)$('#cell-line-results').innerHTML=cellLineTable(scope);}if(e.target.id==='browse-group'){browse.group=e.target.value;browse.page=0;renderBrowseRows();}if(e.target.id==='browse-sort'){browse.sort=e.target.value;browse.page=0;renderBrowseRows();}});
$('#gene-search').addEventListener('input',renderSearch);
$('#gene-search').addEventListener('keydown',e=>{if(e.key==='Escape')hideResults();if(['ArrowDown','ArrowUp'].includes(e.key)&&results.length){e.preventDefault();$('#search-results').hidden=false;highlight=highlight<0?(e.key==='ArrowDown'?0:results.length-1):(highlight+(e.key==='ArrowDown'?1:-1)+results.length)%results.length;$('#gene-search').setAttribute('aria-expanded','true');$('#gene-search').setAttribute('aria-activedescendant',`result-${highlight}`);document.querySelectorAll('.search-option').forEach((el,i)=>el.setAttribute('aria-selected',i===highlight));document.getElementById(`result-${highlight}`)?.scrollIntoView({block:'nearest'});}});
$('#search-form').addEventListener('submit',e=>{e.preventDefault();const q=$('#gene-search').value.trim();if(!q)return;if(highlight>=0&&results[highlight])navigate(results[highlight].symbol);else navigate(q);hideResults();});
document.addEventListener('keydown',e=>{if(e.key==='/'&&!['INPUT','TEXTAREA','SELECT'].includes(document.activeElement.tagName)&&!$('#dialog').open){e.preventDefault();$('#gene-search').focus();}if(e.target.getAttribute('role')==='tab'&&['ArrowLeft','ArrowRight','Home','End'].includes(e.key)){e.preventDefault();const i=tabs.findIndex(t=>t[0]===activeTab),n=e.key==='Home'?0:e.key==='End'?tabs.length-1:(i+(e.key==='ArrowRight'?1:-1)+tabs.length)%tabs.length;navigate(current.symbol,tabs[n][0]);setTimeout(()=>$(`#tab-${tabs[n][0]}`)?.focus(),50);}});
$('.skip').onclick=e=>{e.preventDefault();$('#gene-search').focus();};
$('#search-nav').onclick=()=>{location.hash='view=search';hideResults();$('#gene-search').focus();};$('#browse-nav').onclick=()=>{location.hash='view=browse';};$('#sources-nav').onclick=()=>{location.hash='view=sources';};$('#compare-nav').onclick=showCompare;$('#help-button').onclick=help;$('.dialog-close').onclick=()=>$('#dialog').close();$('#dialog').addEventListener('click',e=>{if(e.target===$('#dialog')){const r=$('#dialog').getBoundingClientRect();if(e.clientX<r.left||e.clientX>r.right||e.clientY<r.top||e.clientY>r.bottom)$('#dialog').close();}});window.addEventListener('hashchange',route);
try{const [idx,m,r,evidenceSources]=await Promise.all([fetchJSON('data/index.json'),fetchJSON('data/manifest.json'),fetchJSON('data/record-sources.json'),fetchJSON('data/evidence-sources.json'),fetch(resourceUrl('site-config.json'),{cache:'no-cache'}).then(r=>{if(!r.ok)throw new Error('Settings unavailable');return r.json();}).then(config=>{siteConfig=validateSettings(config);applySettings(siteConfig,document);})]);genes=unpackIndex(idx);manifest={...m,sources:[...m.sources.map(s=>s.id==='library'?{...s,title:'Cell-line targeting summary',note:'Gene-level counts for splice-distance and exonic candidate sets. Counts overlap across distance settings; source coverage is incomplete.'}:s),...r.sources,...evidenceSources.sources]};compare=compare.filter(s=>genes.some(g=>g.symbol===s));updateCompare();$('#gene-search').disabled=false;$('#search-submit').disabled=false;await route();}catch(e){$('#main').innerHTML='<div class="not-found"><h1>The evidence index could not load</h1><p>Serve the dist folder over HTTP using the project’s start command, then refresh.</p><button id="reload-app">Reload</button></div>';$('#reload-app').onclick=()=>location.reload();}
