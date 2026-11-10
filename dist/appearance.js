// Appearance is browser-local; changing it does not reload the selected gene.
export function initAppearance(document, environment) {
  const key = 'dnd-appearance';
  const select = document.querySelector('#appearance-select');
  const sidebar = document.querySelector('#shiny-sidebar');
  const search = document.querySelector('.search-section');
  const nav = document.querySelector('.site-nav');
  let saved;
  try { saved = environment.localStorage.getItem(key); } catch {}

  function apply(mode) {
    const appearance = mode === 'classic' ? 'classic' : 'shiny';
    document.body.dataset.appearance = appearance;
    if (select) select.value = appearance;
    if (sidebar) sidebar.hidden = appearance !== 'shiny';
    const searchSlot = document.querySelector(appearance === 'shiny' ? '#sidebar-search-slot' : '#search-slot');
    const navSlot = document.querySelector(appearance === 'shiny' ? '#sidebar-nav-slot' : '#header-nav-slot');
    if (search && searchSlot?.append) searchSlot.append(search);
    if (nav && navSlot?.append) navSlot.append(nav);
    try { environment.localStorage.setItem(key, appearance); } catch {}
  }

  apply(saved);
  select?.addEventListener('change', () => apply(select.value));
}
