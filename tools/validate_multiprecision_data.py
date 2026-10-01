#!/usr/bin/env python3
from pathlib import Path
import csv,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1]
manifest=ROOT/'metadata/multiprecision/manifest.sha256'
errors=[]
for line in manifest.read_text(encoding='utf-8').splitlines():
    if not line.strip(): continue
    expected,rel=line.split('  ',1); p=ROOT/rel
    if not p.is_file(): errors.append(f'missing: {rel}'); continue
    if hashlib.sha256(p.read_bytes()).hexdigest()!=expected: errors.append(f'hash mismatch: {rel}')
def count(rel):
    with (ROOT/rel).open(newline='',encoding='utf-8') as f: return sum(1 for _ in csv.reader(f))-1
checks={
'data/multiprecision/families/audit/observations.csv':496,
'data/multiprecision/families/audit/precision_controls.csv':16,
'data/multiprecision/external/audit/raw_observations.csv':736,
'data/multiprecision/external/audit/precision_checks.csv':48}
for rel,n in checks.items():
    got=count(rel)
    if got!=n: errors.append(f'row count {rel}: expected {n}, got {got}')
if errors:
    print('VALIDATION FAILED'); [print('-',e) for e in errors]; sys.exit(1)
print('OK: multiprecision data hashes and published counts validated.')
print('Family observations: 496; controls: 16.')
print('External observations: 736; controls: 48.')
