"""Verify the distributed payload independently of private repository history."""
from pathlib import Path
import csv,hashlib
ROOT=Path(__file__).resolve().parents[1]
def main():
    errors=[]
    with (ROOT/'FILE_MANIFEST_SHA256.csv').open(newline='',encoding='utf-8') as f:rows=list(csv.DictReader(f))
    for row in rows:
        p=ROOT/row['path']
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:errors.append(row['path'])
    if errors:raise SystemExit('Missing/modified release payload: '+', '.join(errors))
    print(f'Verified {len(rows)} release files. Locally regenerated ignored data are not part of this manifest.')
if __name__=='__main__':main()
