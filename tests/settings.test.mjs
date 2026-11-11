import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {validateSettings,applySettings,pageTitle} from '../dist/settings.js';
import {unpackIndex,searchGenes} from '../dist/core.js';

const original=JSON.parse(await readFile(new URL('../dist/site-config.json',import.meta.url),'utf8'));
test('published settings are valid and example genes exist',async()=>{
  validateSettings(original);
  const index=JSON.parse(await readFile(new URL('../dist/data/index.json',import.meta.url),'utf8'));
  const genes=unpackIndex(index);
  for(const query of original.example_genes)assert.equal(searchGenes(genes,query)[0]?.symbol,'SPTLC1',query);
});
test('custom text, hidden status and colors apply without interpreting markup',()=>{
  const config=structuredClone(original);
  config.title='My <b>Gene Browser</b>';
  config.status_label='';
  config.example_genes=[];
  config.theme.link='#123456';
  const elements=new Map(),style=new Map();
  const document={querySelector(selector){if(!elements.has(selector))elements.set(selector,{setAttribute(k,v){this[k]=v;}});return elements.get(selector);},documentElement:{style:{setProperty(k,v){style.set(k,v);}}}};
  applySettings(config,document);
  assert.equal(elements.get('#site-name').textContent,config.title);
  assert.equal(elements.get('#site-name').innerHTML,undefined);
  assert.equal(elements.get('#site-status').hidden,true);
  assert.equal(elements.get('#search-examples').hidden,true);
  assert.equal(style.get('--link'),'#123456');
  assert.equal(pageTitle(config,'PMP22'),'PMP22 · My <b>Gene Browser</b>');
});
test('settings reject invalid CSS, sizes and gene-button markup',()=>{
  for(const mutate of [c=>c.theme.link='red;display:none',c=>c.body_font_size=100,c=>c.content_width='wide',c=>c.example_genes=['<img src=x>'],c=>c.font_family='url(https://example.com)']) {
    const config=structuredClone(original);mutate(config);assert.throws(()=>validateSettings(config));
  }
});
