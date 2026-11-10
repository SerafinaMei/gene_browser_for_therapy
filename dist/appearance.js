// Use one public layout; the Classic version remains preserved in Git.
export function initAppearance(document) {
  const sidebar = document.querySelector('#shiny-sidebar');
  const search = document.querySelector('.search-section');
  const nav = document.querySelector('.site-nav');
  document.body.dataset.appearance = 'shiny';
  if (sidebar) sidebar.hidden = false;
  const navSlot = document.querySelector('#sidebar-nav-slot');
  if (nav && navSlot?.append) navSlot.append(nav);

  function placeSearch() {
    const inSidebar = document.body.dataset.view !== 'home';
    const slot = document.querySelector(inSidebar ? '#sidebar-search-slot' : '#search-slot');
    if (search && slot?.append && search.parentElement !== slot) slot.append(search);
  }

  placeSearch();
  return placeSearch;
}
