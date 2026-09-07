"""Explicitly acquire/check required public inputs and parse local GEO metadata."""
from pathlib import Path
import argparse,csv,hashlib,json,urllib.request
from query_geo_metadata import parse_family_soft,flatten
ROOT=Path(__file__).resolve().parents[1]
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--download-geo',action='store_true');parser.add_argument('--check',action='store_true');parser.add_argument('--parse-geo',action='store_true');args=parser.parse_args()
    sources=json.loads((ROOT/'provenance/required_public_inputs.json').read_text())
    problems=[]
    for row in sources:
        p=ROOT/row['local_file'];p.parent.mkdir(parents=True,exist_ok=True)
        if args.download_geo and row.get('download_url') and not p.exists():
            tmp=p.with_suffix(p.suffix+'.partial')
            request=urllib.request.Request(row['download_url'],headers={'User-Agent':'IBD-mucosal-response-reproduction/1.0'})
            print('Downloading',row['accession'],flush=True)
            with urllib.request.urlopen(request,timeout=60) as response,tmp.open('wb') as target:
                for b in iter(lambda:response.read(1024*1024),b''):target.write(b)
            if digest(tmp)!=row['sha256']:raise RuntimeError(f'Changed source hash; partial file retained for inspection: {tmp.name}')
            tmp.replace(p)
        if not p.is_file():problems.append(f'MISSING {row["local_file"]} — obtain from {row["record_url"]}');continue
        if digest(p)!=row['sha256']:problems.append(f'HASH MISMATCH {row["local_file"]}');continue
        print('Verified',row['local_file'],flush=True)
        if args.parse_geo and row['accession'].startswith('GSE'):
            parsed=parse_family_soft(p);samples=parsed['samples']
            out=p.parent/f'{row["accession"]}_samples.csv'
            fields=['sample_id']+sorted({k for r in samples for k in r if k!='sample_id'})
            with out.open('w',newline='',encoding='utf-8-sig') as f:
                writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows({k:flatten(r.get(k)) for k in fields} for r in samples)
    if problems:raise SystemExit('\n'.join(problems))
    print('All six required source files verified.')
if __name__=='__main__':main()
