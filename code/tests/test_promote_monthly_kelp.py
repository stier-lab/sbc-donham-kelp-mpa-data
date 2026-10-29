"""Fail-closed promotion tests, including preservation of untouched CSV bytes."""
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location("promote", Path(__file__).resolve().parents[1] / "promote_monthly_kelp.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class PromotionTests(unittest.TestCase):
    baseline = (b'CA_MPA_Name_Short,year,taxon_name,source,resp,status,value,area\n'
                b'"A",2007,"Macrocystis pyrifera","KFM","Bio","mpa",80,10\n'
                b'"A",2007,"Other taxon","KFM","Bio","mpa",12.0000000000,10\n')

    def candidate(self, old, new):
        return self.baseline.replace(old, new)

    def test_scoped_change_preserves_other_line(self):
        result, report, changes = module.prepare(self.baseline, self.candidate(b',80,10', b',100,10'), True)
        self.assertEqual(result.splitlines()[-1], self.baseline.splitlines()[-1])
        self.assertEqual(report['changed_rows'], 1)
        self.assertEqual(changes[0]['after']['value'], '100')

    def test_nonkelp_change_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Non-target'):
            module.prepare(self.baseline, self.candidate(b'12.0000000000', b'13'), True)

    def test_metadata_change_rejected(self):
        with self.assertRaisesRegex(ValueError, 'metadata'):
            module.prepare(self.baseline, self.candidate(b',80,10', b',80,11'), True)

    def test_identity_change_rejected(self):
        with self.assertRaisesRegex(ValueError, 'identity'):
            module.prepare(self.baseline, self.candidate(b'2007', b'2008'), True)

    def test_nonfinite_change_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Nonfinite'):
            module.prepare(self.baseline, self.candidate(b',80,10', b',nan,10'), True)

    def test_sheephead_scope_preserves_density_and_other_programmes(self):
        baseline = (b'CA_MPA_Name_Short,year,taxon_name,source,resp,status,value\n'
                    b'"A",2007,"Semicossyphus pulcher","PISCO","Bio","mpa",100\n'
                    b'"A",2007,"Semicossyphus pulcher","PISCO","Den","mpa",2\n'
                    b'"A",2007,"Semicossyphus pulcher","LTER","Bio","mpa",200\n')
        options = dict(taxon_name="Semicossyphus pulcher", sources=("PISCO",), responses=("Bio",))
        candidate = baseline.replace(b',100\n', b',110\n')
        result, report, _ = module.prepare(baseline, candidate, True, **options)
        self.assertEqual(result.splitlines()[2:], baseline.splitlines()[2:])
        self.assertEqual(report["non_target_rows_preserved"], 2)
        for old, new in [(b',2\n', b',3\n'), (b',200\n', b',210\n')]:
            with self.assertRaisesRegex(ValueError, "Non-target"):
                module.prepare(baseline, candidate.replace(old, new), True, **options)


if __name__ == '__main__':
    unittest.main()
