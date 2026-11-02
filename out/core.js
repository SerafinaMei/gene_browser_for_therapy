export const present = v => v !== null && v !== undefined && v !== '';
export const fmt = (v, digits = 3) => !present(v) ? 'Not available' : typeof v === 'number' ? (Number.isInteger(v) ? v.toLocaleString('en-US') : v.toLocaleString('en-US', {maximumSignificantDigits:digits})) : String(v);
export const esc = v => String(v ?? '').replace(/[\u2013\u2014]/g, '-').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
export function unpackIndex(data) { return data.rows.map(row => Object.fromEntries(data.fields.map((f,i) => [f,row[i]]))); }
export function searchGenes(genes, input, limit=12) {
  const q = input.trim().toUpperCase().replace(/^(ENSG\d+)\.\d+$/, '$1');
  if (!q) return [];
  return genes.map(g => {
    const s=g.symbol.toUpperCase(), ids=[g.hgnc_id,g.ensembl_id,g.entrez_id].filter(Boolean).map(String);
    const aliases=(g.aliases||[]).map(a=>a.toUpperCase());
    const rank = s===q ? 100 : ids.includes(q) ? 95 : aliases.includes(q) ? 90 : s.startsWith(q) ? 75 : ids.some(i=>i.startsWith(q)) ? 65 : aliases.some(a=>a.startsWith(q)) ? 50 : s.includes(q) ? 40 : g.name.toUpperCase().includes(q) ? 30 : 0;
    return {g,rank};
  }).filter(x=>x.rank).sort((a,b)=>b.rank-a.rank||a.g.symbol.localeCompare(b.g.symbol)).slice(0,limit).map(x=>x.g);
}
export function flatten(value, prefix='', out={}) {
  if (Array.isArray(value)) {
    if (!value.length) out[prefix]='[]';
    value.forEach((v,i)=>flatten(v,`${prefix}[${i}]`,out));
  } else if (value && typeof value === 'object') {
    Object.entries(value).forEach(([k,v])=>flatten(v,prefix?`${prefix}.${k}`:k,out));
  } else out[prefix]=value;
  return out;
}
export function csv(rows) {
  return rows.map(row=>row.map(v=>{
    let s=present(v)?String(v):'';
    if(typeof v==='string' && /^[=+\-@\t\r]/.test(s))s="'"+s;
    return '"'+s.replaceAll('"','""')+'"';
  }).join(',')).join('\r\n');
}
export function unpackRecordCsv(packet, symbol) {
  const headers = packet.headers, columns = packet.genes?.[symbol];
  if (!Array.isArray(headers) || !Array.isArray(columns) || columns.length !== headers.length ||
      !columns.length || !columns.every(c => Array.isArray(c) && c.length === columns[0].length)) {
    throw new Error('Invalid download record');
  }
  const cell = value => {
    if (typeof value !== 'string') throw new Error('Invalid CSV value');
    if (/^[=+\-@\t\r]/.test(value)) value = "'" + value;
    return /[",\r\n]/.test(value) ? '"' + value.replaceAll('"', '""') + '"' : value;
  };
  const lines = [headers.map(cell).join(',')];
  for (let i = 0; i < columns[0].length; i++) lines.push(columns.map(c => cell(c[i])).join(','));
  return '\ufeff' + lines.join('\r\n') + '\r\n';
}
export const hiLabels={0:'No evidence',1:'Little evidence',2:'Some evidence',3:'Sufficient evidence',30:'Autosomal recessive',40:'Dosage sensitivity unlikely'};
export const moiLabels={AD:'Autosomal dominant',AR:'Autosomal recessive',SD:'Semidominant',XL:'X-linked',XLD:'X-linked dominant',XLR:'X-linked recessive',MT:'Mitochondrial'};

export function genomeUrl(g, tracks = false) {
  const url = new URL('https://genome.ucsc.edu/cgi-bin/hgTracks');
  url.searchParams.set('db', 'hg38');
  const coordinates = (g.coordinates || '').replace(/^chr/, '').replaceAll(',', '');
  const valid = /^([0-9]+|X|Y|M|MT):\d+-\d+$/.test(coordinates);
  url.searchParams.set('position', valid ? `chr${coordinates}` : g.symbol);
  if (!valid) url.searchParams.set('singleSearch', 'knownCanonical');
  url.searchParams.set('virtModeType', 'default');
  url.searchParams.set('ignoreCookie', '1');
  if (tracks) {
    // Load the supplied track configuration, then override its saved location.
    url.searchParams.set('hgS_doOtherUser', 'submit');
    url.searchParams.set('hgS_otherUserName', 'gramey02');
    url.searchParams.set('hgS_otherUserSessionName', 'SPTLC_gene_session');
  }
  return url.href;
}
