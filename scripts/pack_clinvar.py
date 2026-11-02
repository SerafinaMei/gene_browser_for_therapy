"""Losslessly pack generated ClinVar CSVs for static hosting and CSV downloads."""
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / 'dist'


def csv_bytes(headers, columns):
    stream = io.StringIO(newline='')
    writer = csv.writer(stream)
    writer.writerow(headers)
    writer.writerows(zip(*columns))
    return stream.getvalue().encode('utf-8-sig')


def pack_clinvar():
    shards = {p: json.loads(p.read_text()) for p in sorted((OUT / 'data').glob('record-details-*.json'))}
    replacements, originals, used = {}, [], set()
    old_size = new_size = count = 0
    for shard, details in shards.items():
        pending, pending_size, batch = {}, 0, 0
        headers = None

        def flush():
            nonlocal pending, pending_size, batch, new_size
            if not pending:
                return
            batch += 1
            key = shard.stem.removeprefix('record-details-')
            relative = f'downloads/clinvar-packed/{key}-{batch:04d}.json.gz'
            path = OUT / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            payload = json.dumps({'headers': headers, 'genes': pending}, separators=(',', ':')).encode()
            packed = gzip.compress(payload, compresslevel=9, mtime=0)
            decoded = json.loads(gzip.decompress(packed))
            for symbol, columns in decoded['genes'].items():
                digest = hashlib.sha256(csv_bytes(decoded['headers'], columns)).hexdigest()
                assert digest == replacements[symbol]['sha256'], symbol
                replacements[symbol].update(path=relative, format='columnar-gzip', record_key=symbol)
            path.write_bytes(packed)
            used.add(path)
            new_size += len(packed)
            pending, pending_size = {}, 0

        for symbol, detail in sorted(details.items()):
            record = detail.get('clinvar_export')
            if not record:
                continue
            if record.get('format') == 'columnar-gzip':
                used.add(OUT / record['path'])
                continue
            path = OUT / record['path']
            assert path.is_relative_to(OUT / 'downloads/clinvar')
            raw = gzip.decompress(path.read_bytes())
            rows = list(csv.reader(io.StringIO(raw.decode('utf-8-sig'), newline='')))
            assert len(rows) - 1 == record['rows']
            if headers is None:
                headers = rows[0]
            assert headers == rows[0]
            assert all(len(row) == len(headers) for row in rows)
            columns = list(map(list, zip(*rows[1:])))
            assert csv_bytes(headers, columns) == raw, symbol
            size = len(json.dumps(columns, separators=(',', ':')))
            if pending and pending_size + size > 1024 * 1024:
                flush()
            pending[symbol] = columns
            pending_size += size
            replacements[symbol] = dict(record, sha256=hashlib.sha256(raw).hexdigest())
            originals.append(path)
            old_size += path.stat().st_size
            count += 1
        flush()
    if not count:
        print('ClinVar downloads are already packed.', flush=True)
        return
    for path, details in shards.items():
        for symbol, detail in details.items():
            if symbol in replacements:
                detail['clinvar_export'] = replacements[symbol]
        temporary = path.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(details, separators=(',', ':')))
        temporary.replace(path)
    # Only generated CSVs verified byte-for-byte above are removed. Original
    # ClinVar inputs outside this website are never modified.
    for path in originals:
        path.unlink()
    for path in (OUT / 'downloads/clinvar-packed').glob('*.json.gz'):
        if path not in used:
            path.unlink()
    print(f'Packed {count:,} gene exports in {len(used):,} files: {old_size/1024**2:.1f} -> {new_size/1024**2:.1f} MiB. Every CSV verified byte-for-byte.', flush=True)


if __name__ == '__main__':
    pack_clinvar()
