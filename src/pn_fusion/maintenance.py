"""File-integrity manifest for repository release artifacts."""
from pathlib import Path
import hashlib
import json
import re
import tomllib
from .io import ROOT

EXCLUDED={'.git','.venv','venv','__pycache__','.pytest_cache','build','dist'}

def _matched_version(path: Path, pattern: str):
    match=re.search(pattern,path.read_text(encoding='utf-8'),re.MULTILINE)
    if not match:
        raise RuntimeError(f'Cannot determine release version from {path.name}')
    return match.group(1)

def verify_version_metadata(root: Path=ROOT):
    with (root/'pyproject.toml').open('rb') as stream:
        project_version=str(tomllib.load(stream)['project']['version'])
    versions={
        'pyproject.toml':project_version,
        'CITATION.cff':_matched_version(root/'CITATION.cff',r'^version:\s*["\']?([^"\'\s#]+)["\']?\s*(?:#.*)?$'),
        'CHANGELOG.md':_matched_version(root/'CHANGELOG.md',r'^##\s+(\d+\.\d+\.\d+(?:[-+][0-9A-Za-z.-]+)?)\b'),
        'metadata/zenodo_metadata_template.json':str(json.loads((root/'metadata/zenodo_metadata_template.json').read_text(encoding='utf-8'))['version']),
        'src/pn_fusion/__init__.py':_matched_version(root/'src/pn_fusion/__init__.py',r'^__version__\s*=\s*["\']([^"\']+)["\']\s*$'),
    }
    if len(set(versions.values())) != 1:
        details=', '.join(f'{name}={version}' for name,version in versions.items())
        raise RuntimeError(f'Release version mismatch: {details}')
    return project_version

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
    version=verify_version_metadata(root)
    expected={}
    for line in (root/'metadata/SHA256SUMS.txt').read_text().splitlines():
        digest,name=line.split('  ',1);expected[name]=digest
    found={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in release_files(root)}
    missing=sorted(set(expected)-set(found));new=sorted(set(found)-set(expected))
    changed=sorted(name for name in set(expected)&set(found) if expected[name]!=found[name])
    if missing or new or changed:
        raise RuntimeError(f'Release integrity mismatch: missing={missing}, new={new}, changed={changed}')
    return {'verified_files':len(expected),'version':version,'passed':True}
