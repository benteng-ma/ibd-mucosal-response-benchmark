"""Portable staged wrapper; original scientific scripts remain unchanged."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
def run(script,*args):subprocess.run([sys.executable,str(ROOT/'scripts'/script),*args],cwd=ROOT,check=True)
def verify_frozen_inputs():
    expected=json.loads((ROOT/'provenance/frozen_input_sha256.json').read_text())
    for name,digest in expected.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Frozen input changed: '+name)
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['phase1','phase3','figures','full'],required=True);args=parser.parse_args()
    verify_frozen_inputs()
    if args.stage in ['phase1','full']:
        run('prepare_public_inputs.py','--check','--parse-geo')
        for s in ['extract_phase1_expression.py','build_phase1_analysis_manifests.py','run_phase1_benchmark.py','verify_phase1_inputs.py','build_phase1_sensitivity_audit.py']:run(s)
        expected=json.loads((ROOT/'provenance/frozen_result_sha256.json').read_text())
        for name,digest in expected.items():
            if name.startswith('results/performance/') and hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Recreated aggregate differs from frozen reference: '+name)
        print('Phase 1 aggregate files match frozen reference exactly.',flush=True)
    if args.stage in ['phase3','full']:run('run_phase3_robustness.py')
    if args.stage in ['figures','full']:
        run('plot_phase3_results.py');run('build_phase3_supplementary_tables.py')
    print('Completed',args.stage,'. Individual-level regenerated outputs remain Git-ignored.',flush=True)
if __name__=='__main__':main()
