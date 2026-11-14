from pathlib import Path
import hashlib,json,ast
r=Path(__file__).resolve().parent
m=json.loads((r/'SOURCE_MANIFEST.json').read_text())
for name,h in m['files'].items():
 p=r/name
 assert p.is_file(),name
 assert hashlib.sha256(p.read_bytes()).hexdigest()==h,name
for name in m['new_python_sources']:ast.parse((r/name).read_text(),filename=name)
print('All source hashes verified; new authoring Python parsed. No analysis executed.')
