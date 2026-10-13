"""Small checks for population matching, localization, and dosage LD."""
import csv
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = Path(os.environ.get("PAIR_LD_SCRIPT", ROOT / "scripts/pair_ld.py"))
sys.path.insert(0, str(SCRIPT.parent))


class PairLDTest(unittest.TestCase):
    def test_cloud_uri_is_a_localization_error(self):
        from pair_ld import require_local
        with self.assertRaisesRegex(ValueError, "localization"):
            require_local("gs://bucket/cohort.pgen")

    def test_population_table_must_cover_the_cohort(self):
        from pair_ld import read_groups
        with tempfile.TemporaryDirectory() as directory:
            psam = Path(directory) / "test.psam"
            populations = Path(directory) / "populations.tsv"
            psam.write_text("#IID\tSEX\nS1\t1\nS2\t2\n")
            populations.write_text("FID\tIID\tpopulation\n0\tS1\tA\n")
            with self.assertRaisesRegex(ValueError, "sample.*population"):
                read_groups(psam, populations)

    def test_dosage_ld_for_full_cohort_and_populations(self):
        plink = os.environ.get("PLINK2", shutil.which("plink2"))
        self.assertTrue(plink, "Set PLINK2 or put plink2 on PATH for this tiny check")
        with tempfile.TemporaryDirectory(prefix="pair ld '") as directory:
            work = Path(directory)
            samples = [f"S{i}" for i in range(120)]
            dose_a = [(.25, 1., 1.75)[i % 3] for i in range(120)]
            dose_b = [x if i < 60 else 2 - x for i, x in enumerate(dose_a)]
            vcf = work / "tiny.vcf"
            with vcf.open("w") as output:
                output.write('##fileformat=VCFv4.2\n')
                output.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                output.write('##FORMAT=<ID=DS,Number=1,Type=Float,Description="ALT dosage">\n')
                output.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
                for chrom, pos, variant, doses in [(1, 100, "v1", dose_a), (1, 2000100, "v2", dose_b), (1, 3000100, "v3", [0.] * 120), (2, 100, "v4", dose_a)]:
                    output.write(f"{chrom}\t{pos}\t{variant}\tA\tG\t.\tPASS\t.\tGT:DS\t" + "\t".join(f"0/0:{dose}" for dose in doses) + "\n")
            subprocess.run([plink, "--vcf", str(vcf), "dosage=DS", "--make-pgen", "--out", str(work / "tiny")], check=True, capture_output=True)
            # Keep known parents in the source PSAM. Every sample must still be used.
            (work / "tiny.psam").write_text("#IID\tPAT\tMAT\n" + "".join(f"{sample}\t{'0' if i < 2 else 'S0'}\t{'0' if i < 2 else 'S1'}\n" for i, sample in enumerate(samples)))
            populations = work / "populations.tsv"
            populations.write_text("FID\tIID\tpopulation\n" + "".join(f"0\t{sample}\t{'A' if i < 60 else 'B'}\n" for i, sample in enumerate(samples)))
            pairs = work / "pairs.tsv"
            pairs.write_text("variant_a\tvariant_b\nv1\tv2\nv1\tv3\nv1\tabsent\nv2\tv1\nv1\tv4\n")
            result = subprocess.run([sys.executable, str(SCRIPT), "--pgen", str(work / "tiny.pgen"), "--pvar", str(work / "tiny.pvar"), "--psam", str(work / "tiny.psam"), "--pairs", str(pairs), "--populations", str(populations), "--plink", plink, "--out-dir", str(work / "result")], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with (work / "result/ld.tsv").open() as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(len(rows), 15)
            for group, count, expected in [("ALL", 120, 0.), ("A", 60, 1.), ("B", 60, 1.)]:
                chosen = [r for r in rows if r["population"] == group]
                self.assertEqual([r["n_samples"] for r in chosen], [str(count)] * 5)
                self.assertAlmostEqual(float(chosen[0]["r2"]), expected, places=6)
                self.assertEqual(chosen[1]["status"], "undefined_ld")
                self.assertEqual(chosen[2]["status"], "missing_variant")
                self.assertAlmostEqual(float(chosen[3]["r2"]), expected, places=6)
                self.assertAlmostEqual(float(chosen[4]["r2"]), 1., places=6)


if __name__ == "__main__":
    unittest.main()
