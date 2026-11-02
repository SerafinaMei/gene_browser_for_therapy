import {esc} from './core.js?v=6';

const fonts = {
  system: '-apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif',
  arial: 'Arial,Helvetica,sans-serif',
  verdana: 'Verdana,Geneva,sans-serif'
};
const textKeys = ['title','organization','status_label','footer_note','search_heading','description','search_label','search_placeholder','search_help'];
const colors = {text:'--ink',muted:'--muted',border:'--line',surface:'--soft',link:'--link',background:'--background',primary:'--primary'};

export function validateSettings(config) {
  if (!config || typeof config !== 'object' || Array.isArray(config)) throw new Error('Settings must be a JSON object.');
  for (const key of textKeys) {
    if (typeof config[key] !== 'string' || config[key].length > 300) throw new Error(`Invalid setting: ${key}`);
  }
  if (!config.title.trim() || !config.search_heading.trim()) throw new Error('The title and search heading cannot be empty.');
  if (!Array.isArray(config.example_genes) || config.example_genes.length > 8 ||
      !config.example_genes.every(s => typeof s === 'string' && /^[A-Za-z0-9_.-]{1,32}$/.test(s))) throw new Error('Use up to eight gene symbols for examples.');
  if (!Object.hasOwn(fonts, config.font_family)) throw new Error('font_family must be system, arial or verdana.');
  if (!Number.isInteger(config.body_font_size) || config.body_font_size < 12 || config.body_font_size > 20) throw new Error('body_font_size must be an integer from 12 to 20.');
  if (!Number.isInteger(config.content_width) || config.content_width < 900 || config.content_width > 1600) throw new Error('content_width must be an integer from 900 to 1600.');
  for (const key of Object.keys(colors)) {
    if (!/^#[0-9a-f]{6}$/i.test(config.theme?.[key])) throw new Error(`theme.${key} must be a six-digit hex color.`);
  }
  return config;
}

export function pageTitle(config, page='') {
  return [page, config.title, config.status_label].filter(Boolean).join(' · ');
}

export function applySettings(config, document) {
  validateSettings(config);
  const setText = (selector, value) => { const element=document.querySelector(selector); if(element) element.textContent=value; };
  setText('#site-name', config.title);
  setText('#site-status', config.status_label);
  const status=document.querySelector('#site-status');
  if (status) status.hidden=!config.status_label;
  setText('#search-heading', config.search_heading);
  setText('#search-description', config.description);
  setText('#search-label', config.search_label);
  setText('#search-help', config.search_help);
  setText('#site-footer', [config.organization,config.status_label,config.footer_note].filter(Boolean).join(' · '));
  document.querySelector('#gene-search')?.setAttribute('placeholder',config.search_placeholder);
  document.querySelector('meta[name="description"]')?.setAttribute('content',config.description);
  const examples=document.querySelector('#search-examples');
  if (examples) {
    examples.hidden=config.example_genes.length===0;
    examples.innerHTML='<span>Examples:</span>'+config.example_genes.map(s=>`<button data-gene="${esc(s)}">${esc(s)}</button>`).join('');
  }
  const style=document.documentElement?.style;
  if (style) {
    for (const [key,variable] of Object.entries(colors)) style.setProperty(variable,config.theme[key]);
    style.setProperty('--font',fonts[config.font_family]);
    style.setProperty('--body-size',`${config.body_font_size}px`);
    style.setProperty('--content-width',`${config.content_width}px`);
  }
  document.title=pageTitle(config);
}
