"""Serve the prebuilt explorer using only Python's standard library."""
import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--port', type=int, default=8765)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1] / 'dist'
handler = partial(SimpleHTTPRequestHandler, directory=str(root))
with ThreadingHTTPServer(('127.0.0.1', args.port), handler) as server:
    print(f'D&D Gene Explorer: http://127.0.0.1:{args.port}/', flush=True)
    print('Press Ctrl+C to stop. Serving only the dist directory.', flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
