"""UTF-8 CSV and JSON I/O; reference data are always read-only."""
from __future__ import annotations
from pathlib import Path
import csv
import json

ROOT=Path(__file__).resolve().parents[2]

def read_csv(name: str, root: Path|None=None) -> list[dict[str,str]]:
    root=ROOT if root is None else root
    path=root/'data/reference'/name
    with path.open(encoding='utf-8-sig',newline='') as stream:
        return list(csv.DictReader(stream))

def write_csv(path: Path, rows: list[dict], fields: list[str]|None=None) -> None:
    path=Path(path)
    if (ROOT/'data/reference').resolve() in path.resolve().parents:
        raise ValueError('Refusing to overwrite archived reference data')
    path.parent.mkdir(parents=True,exist_ok=True)
    if fields is None:
        if not rows: raise ValueError('fields must be supplied for an empty CSV')
        fields=list(rows[0])
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)

def write_json(path: Path, value) -> None:
    path = Path(path)
    if (ROOT/'data/reference').resolve() in path.resolve().parents:
        raise ValueError('Refusing to overwrite archived reference data')
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
