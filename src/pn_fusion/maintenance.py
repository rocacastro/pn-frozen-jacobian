"""File integrity manifest for a local release candidate."""
from pathlib import Path
import hashlib
from .io import ROOT

EXCLUDED={'.git','.venv','venv','__pycache__','.pytest_cache','build','dist'}

def release_files(root: Path=ROOT):
    for p in sorted(root.rglob('*')):
        rel=p.relative_to(root)
        if p.is_file() and not any(part in EXCLUDED or part.endswith('.egg-info') for part in rel.parts) and p.name!='SHA256SUMS.txt' and p.suffix not in {'.pyc','.aux','.log','.out'}:
            yield p

def write_manifest(root: Path=ROOT):
    lines=[hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(root).as_posix() for p in release_files(root)]
    (root/'metadata/SHA256SUMS.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return {'files':len(lines)}

def verify_manifest(root: Path=ROOT):
    expected={}
    for line in (root/'metadata/SHA256SUMS.txt').read_text().splitlines():
        digest,name=line.split('  ',1);expected[name]=digest
    found={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in release_files(root)}
    missing=sorted(set(expected)-set(found));new=sorted(set(found)-set(expected))
    changed=sorted(name for name in set(expected)&set(found) if expected[name]!=found[name])
    if missing or new or changed:
        raise RuntimeError(f'Release integrity mismatch: missing={missing}, new={new}, changed={changed}')
    return {'verified_files':len(expected),'passed':True}
