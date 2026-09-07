"""Independent checks that run without distributing individual-level data."""
from pathlib import Path
import csv,hashlib,json,math,random,sys,unittest
from sklearn.metrics import average_precision_score,roc_auc_score
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from run_phase1_benchmark import auroc,average_precision,permutation_p,holm_adjust

def rows(path):
    with (ROOT/path).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

class PublicSnapshotTests(unittest.TestCase):
    def test_synthetic_metrics_against_independent_library(self):
        labels=[0,1,0,1,1,0,0,1]
        for scores in ([0,1,2,3,4,5,6,7],[1]*8,[0,1,1,2,2,2,3,3],[7,6,5,4,3,2,1,0]):
            with self.subTest(scores=scores):
                self.assertAlmostEqual(auroc(labels,scores),roc_auc_score(labels,scores),places=12)
                self.assertAlmostEqual(average_precision(labels,scores),average_precision_score(labels,scores),places=12)
        self.assertTrue(math.isnan(auroc([1,1],[0,1])))

    def test_permutation_is_two_sided_with_add_one(self):
        labels=[0,0,0,1,1,1];scores=[0,1,2,3,4,5];rng=random.Random(20260831)
        shuffled=labels[:];exceed=0;observed=abs(roc_auc_score(labels,scores)-0.5)
        for _ in range(199):
            rng.shuffle(shuffled)
            exceed+=abs(roc_auc_score(shuffled,scores)-0.5)>=observed-1e-15
        self.assertEqual(permutation_p(labels,scores,199,20260831),(exceed+1)/200)
        self.assertEqual(permutation_p(labels,scores,199,20260831),permutation_p(labels,[-s for s in scores],199,20260831))

    def test_holm_known_example(self):
        values=[{'permutation_p_value':str(p)} for p in [.04,.01,.03]]
        holm_adjust(values)
        self.assertEqual([float(v['holm_adjusted_p_value']) for v in values],[.06,.03,.06])

    def test_frozen_input_and_result_checksums(self):
        for record in ['frozen_input_sha256.json','frozen_result_sha256.json']:
            expected=json.loads((ROOT/'provenance'/record).read_text())
            for name,digest in expected.items():
                with self.subTest(path=name):
                    self.assertEqual(hashlib.sha256((ROOT/name).read_bytes()).hexdigest(),digest)

    def test_aggregate_partition_and_confirmatory_family(self):
        folder='results/performance/phase1_v1/'
        summary=json.loads((ROOT/folder/'benchmark_summary.json').read_text())
        self.assertEqual(summary['n_executed_cells'],46)
        self.assertEqual(summary['transport_classifications'],{'SUPPORTED_DIRECTION':19,'UNCERTAIN_OR_NONTRANSPORTABLE':27})
        strict=rows(folder+'strict_external_evaluations.csv')
        self.assertEqual(len(strict),3)
        self.assertTrue(all(float(r['holm_adjusted_p_value'])>=.05 for r in strict))
        self.assertEqual(len({r['dataset_id'] for r in rows(folder+'study_specific_evaluations.csv')}),5)

    def test_expected_inputs_are_hash_identified(self):
        inputs=json.loads((ROOT/'provenance/required_public_inputs.json').read_text())
        self.assertEqual(len(inputs),6)
        self.assertEqual(len({r['local_file'] for r in inputs}),6)
        for row in inputs:
            self.assertEqual(len(row['sha256']),64)
            self.assertTrue(row['record_url'].startswith('https://'))

    def test_no_claim_of_public_preregistration(self):
        provenance=json.loads((ROOT/'provenance/source_commits.json').read_text())
        self.assertFalse(provenance['history_exported'])
        self.assertFalse(provenance['public_preregistration_claim'])

if __name__=='__main__':unittest.main(verbosity=2)
