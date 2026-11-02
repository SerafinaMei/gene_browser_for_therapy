"""Build a compact static deployment without changing the local data files."""
import gzip
import hashlib
import json
import re
import shutil
from html import escape
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / 'dist'
OUTPUT = PROJECT / 'out'


def main():
    if OUTPUT.exists():
        shutil.rmtree(OUTPUT)
    shutil.copytree(SOURCE, OUTPUT, ignore=shutil.ignore_patterns('.DS_Store', '.git', '.openai'))
    fingerprint = hashlib.sha256()
    for path in sorted(list(SOURCE.glob('*.js')) + list(SOURCE.glob('*.css')) + list(SOURCE.glob('*.json')) + list((SOURCE / 'data').glob('*.json'))):
        fingerprint.update(str(path.relative_to(SOURCE)).encode())
        fingerprint.update(path.read_bytes())
    revision = fingerprint.hexdigest()[:16]
    compressed = 0
    for path in sorted((OUTPUT / 'data').glob('*.json')):
        original = path.read_bytes()
        packed = gzip.compress(original, compresslevel=9, mtime=0)
        assert gzip.decompress(packed) == original
        path.with_suffix('.json.gz').write_bytes(packed)
        path.unlink()
        compressed += 1
    index = OUTPUT / 'index.html'
    html = index.read_text()
    assert '<html lang="en">' in html
    html = html.replace('<html lang="en">', f'<html lang="en" data-json-encoding="gzip" data-build-version="{revision}">', 1)
    config = json.loads((SOURCE / 'site-config.json').read_text())
    title = ' · '.join(filter(None, [config['title'], config['status_label']]))
    html = re.sub(r'<title>.*?</title>', lambda _: '<title>' + escape(title) + '</title>', html, count=1)
    for key, value in {'site-name': config['title'], 'site-status': config['status_label'], 'search-heading': config['search_heading'], 'search-description': config['description'], 'search-label': config['search_label'], 'search-help': config['search_help'], 'site-footer': ' · '.join(filter(None, [config['organization'], config['status_label'], config['footer_note']]))}.items():
        html = re.sub(r'(<[^>]+id="' + re.escape(key) + r'"[^>]*>)[^<]*(</[^>]+>)', lambda m: m[1] + escape(value) + m[2], html, count=1)
    html = re.sub(r'(<meta name="description" content=")[^"]*(">)', lambda m: m[1] + escape(config['description'], quote=True) + m[2], html, count=1)
    index.write_text(re.sub(r'\?v=[A-Za-z0-9]+', '?v=' + revision, html))
    for path in OUTPUT.glob('*.js'):
        path.write_text(re.sub(r'\?v=[A-Za-z0-9]+', '?v=' + revision, path.read_text()))
    (OUTPUT / '.nojekyll').write_text('')
    size = sum(p.stat().st_size for p in OUTPUT.rglob('*') if p.is_file())
    assert size < 900 * 1024**2, f'Hosting output too large: {size / 1024**2:.1f} MiB'
    print(f'Hosting output: {size / 1024**2:.1f} MiB; {compressed} annotation files compressed and verified.', flush=True)


if __name__ == '__main__':
    main()
