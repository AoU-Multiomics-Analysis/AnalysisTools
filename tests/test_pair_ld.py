"""Small checks for population matching, localization, and dosage LD."""
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = Path(os.environ.get("PAIR_LD_SCRIPT", ROOT / "scripts/pair_ld.py"))
sys.path.insert(0, str(SCRIPT.parent))


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def reference_per_anchor_ld(plink, prefix, work, pairs, groups, supported, allow_small):
    """Run real PLINK with the old per-anchor calculation, without production helpers."""
    work.mkdir()
    founders = work / "founders.psam"
    founders.write_text("#FID\tIID\n" + "".join(f"0\t{sample}\n" for sample in groups["ALL"]))
    eligible = {tuple(sorted((a, b))) for a, b in pairs if a in supported and b in supported}
    anchors = {}
    for a, b in sorted(eligible):
        anchors.setdefault(a, []).append(b)
    results = {}
    for group_index, (label, samples) in enumerate(groups.items()):
        if len(samples) < (3 if allow_small else 50):
            continue
        keep = work / f"group_{group_index}.keep"
        keep.write_text("#FID\tIID\n" + "".join(f"0\t{sample}\n" for sample in samples))
        values = {}
        for anchor_index, (anchor, partners) in enumerate(anchors.items()):
            extract = work / f"anchor_{anchor_index}.txt"
            extract.write_text("\n".join([anchor, *partners]) + "\n")
            output = work / f"reference_{group_index}_{anchor_index}"
            command = [plink, "--pgen", str(prefix.with_suffix(".pgen")),
                       "--pvar", str(prefix.with_suffix(".pvar")), "--psam", str(founders),
                       "--keep", str(keep), "--extract", str(extract),
                       "--r2-unphased", "inter-chr", "cols=chrom,pos,id", "--ld-snp", anchor,
                       "--ld-window-r2", "-1", "--threads", "1", "--memory", "640", "--out", str(output)]
            if allow_small and len(samples) < 50:
                command.append("--bad-ld")
            result = subprocess.run(command, capture_output=True, text=True)
            if result.returncode:
                raise AssertionError("Reference PLINK failed:\n" + result.stdout + result.stderr)
            for row in read_rows(output.with_suffix(".vcor")):
                values[tuple(sorted((row["ID_A"], row["ID_B"])))] = float(row["UNPHASED_R2"])
        if set(values) != eligible:
            raise AssertionError(f"Reference report omitted pairs: {eligible - set(values)}")
        results[label] = values
    return results


class PairLDTest(unittest.TestCase):
    def make_regression_fixture(self, plink, work):
        samples = [f"S{i}" for i in range(126)]
        groups = {"ALL": samples, "A": samples[:60], "B": samples[60:120],
                  "TINY": samples[120:124], "TOO_SMALL": samples[124:]}
        dose = [(.25, 1., 1.75)[i % 3] for i in range(126)]
        orthogonal = [(.25, 1., 1.75)[(i // 3) % 3] for i in range(126)]
        second_anchor = [value if i < 60 else 2 - value for i, value in enumerate(dose)]
        missing = [None if i % 5 == 0 else value for i, value in enumerate(dose)]
        variants = [
            (1, 100, "b_early", "A", "G", dose),
            (1, 200, "c_low", "A", "G", orthogonal),
            (1, 500, "d_second_anchor", "A", "G", second_anchor),
            (1, 700, "e_monomorphic", "A", "G", [0.] * 126),
            (1, 800, "f_missing", "A", "G", missing),
            (1, 1000, "h_indel", "A", "AT", dose),
            (1, 2100000, "a_late", "A", "G", dose),
            (2, 100, "g_other_chr", "A", "G", orthogonal),
            (0, 150, "x_unsupported", "A", "G", dose),
        ]
        vcf = work / "regression.vcf"
        with vcf.open("w") as stream:
            stream.write('##fileformat=VCFv4.2\n')
            stream.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            stream.write('##FORMAT=<ID=DS,Number=A,Type=Float,Description="ALT dosage">\n')
            stream.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
            for chrom, pos, variant, ref, alt, values in variants:
                fields = ["./.:." if value is None else f"0/0:{value}" for value in values]
                stream.write(f"{chrom}\t{pos}\t{variant}\t{ref}\t{alt}\t.\tPASS\t.\tGT:DS\t"
                             + "\t".join(fields) + "\n")
        prefix = work / "regression"
        imported = subprocess.run([plink, "--vcf", str(vcf), "dosage=DS", "--make-pgen",
                                   "--out", str(prefix)], capture_output=True, text=True)
        self.assertEqual(imported.returncode, 0, imported.stdout + imported.stderr)
        prefix.with_suffix(".psam").write_text("#IID\tPAT\tMAT\n" + "".join(
            f"{sample}\t{'0' if i < 2 else 'S0'}\t{'0' if i < 2 else 'S1'}\n" for i, sample in enumerate(samples)))
        populations = work / "populations.tsv"
        populations.write_text("FID\tIID\tpopulation\n" + "".join(
            f"0\t{sample}\t{label}\n" for label, chosen in groups.items() if label != "ALL" for sample in chosen))
        supported = {"a_late", "b_early", "c_low", "d_second_anchor", "e_monomorphic", "f_missing", "g_other_chr", "h_indel"}
        return prefix, populations, groups, supported

    def test_batched_ld_matches_per_anchor_reference_with_one_real_launch_per_group(self):
        # This catches per-anchor process launches, loss of low/undefined LD, or
        # wrong anchor direction/window bounds in the real batched calculation.
        plink = os.environ.get("PLINK2", shutil.which("plink2"))
        self.assertTrue(plink, "Set PLINK2 or put plink2 on PATH for this regression")
        with tempfile.TemporaryDirectory(prefix="pair ld batch '") as directory:
            work = Path(directory)
            prefix, populations, groups, supported = self.make_regression_fixture(plink, work)
            base_pairs = [("a_late", "b_early"), ("a_late", "c_low"),
                          ("a_late", "e_monomorphic"), ("a_late", "f_missing"),
                          ("d_second_anchor", "b_early"), ("a_late", "d_second_anchor"),
                          ("a_late", "h_indel"), ("b_early", "a_late"),
                          ("a_late", "b_early"), ("a_late", "absent"),
                          ("a_late", "x_unsupported"),
                          ("c_low", "f_missing")]
            for cross_chromosome, allow_small in [(False, False), (True, True)]:
                with self.subTest(cross_chromosome=cross_chromosome, allow_small=allow_small):
                    label = f"cross_{cross_chromosome}_small_{allow_small}"
                    pairs = base_pairs + ([("a_late", "g_other_chr")] if cross_chromosome else [])
                    pairs_path = work / f"{label}.tsv"
                    pairs_path.write_text("variant_a\tvariant_b\n" + "".join(f"{a}\t{b}\n" for a, b in pairs))
                    reference = reference_per_anchor_ld(plink, prefix, work / f"reference_{label}", pairs,
                                                        groups, supported, allow_small)
                    output = work / f"result_{label}"
                    command = [sys.executable, str(SCRIPT), "--pgen", str(prefix.with_suffix(".pgen")),
                               "--pvar", str(prefix.with_suffix(".pvar")), "--psam", str(prefix.with_suffix(".psam")),
                               "--pairs", str(pairs_path), "--populations", str(populations), "--plink", plink,
                               "--threads", "1", "--memory-mb", "640", "--out-dir", str(output)]
                    if allow_small:
                        command.append("--allow-small-groups")
                    result = subprocess.run(command, capture_output=True, text=True)
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    rows = read_rows(output / "ld.tsv")
                    self.assertEqual(len(rows), len(pairs) * len(groups))
                    for group, samples in groups.items():
                        chosen = [row for row in rows if row["population"] == group]
                        self.assertEqual([(row["variant_a"], row["variant_b"]) for row in chosen], pairs)
                        self.assertEqual([row["pair_index"] for row in chosen], [str(i) for i in range(1, len(pairs) + 1)])
                        self.assertEqual([row["n_samples"] for row in chosen], [str(len(samples))] * len(pairs))
                        self.assertEqual([row["small_group"] for row in chosen], [str(len(samples) < 50).lower()] * len(pairs))
                        for row, (a, b) in zip(chosen, pairs):
                            with self.subTest(group=group, pair=(a, b)):
                                if b == "absent":
                                    expected_status = "missing_variant"
                                elif b == "x_unsupported":
                                    expected_status = "unsupported_chromosome"
                                elif group not in reference:
                                    expected_status = "too_few_samples"
                                else:
                                    value = reference[group][tuple(sorted((a, b)))]
                                    expected_status = "ok" if math.isfinite(value) else "undefined_ld"
                                    if math.isfinite(value):
                                        self.assertAlmostEqual(float(row["r2"]), value, delta=1e-10)
                                self.assertEqual(row["status"], expected_status)
                                if expected_status != "ok":
                                    self.assertEqual(row["r2"], "NA")
                        if group in reference:
                            self.assertAlmostEqual(float(chosen[0]["r2"]), 1., places=6)
                            if len(samples) >= 50:
                                self.assertAlmostEqual(float(chosen[1]["r2"]), 0., places=6)
                            self.assertEqual(chosen[2]["status"], "undefined_ld")
                            self.assertAlmostEqual(float(chosen[3]["r2"]), 1., places=6)
                    enabled = [group for group in groups if group in reference]
                    with tarfile.open(output / "plink_logs.tar.gz", "r:gz") as archive:
                        ld_logs = [(member.name, archive.extractfile(member).read().decode())
                                   for member in archive.getmembers() if member.isfile()]
                    ld_logs = [(name, text) for name, text in ld_logs if "--r2-unphased" in text]
                    self.assertEqual(len(ld_logs), len(enabled),
                                     f"Expected one real PLINK LD launch per enabled group {enabled}; logs: {[name for name, _ in ld_logs]}")

    def test_cloud_uri_is_a_localization_error(self):
        from pair_ld import require_local
        with self.assertRaisesRegex(ValueError, "localization"):
            require_local("gs://bucket/cohort.pgen")

    def test_multiallelic_pairs_keep_the_unsupported_status(self):
        from pair_ld import pair_status
        variants = {"biallelic": {"CHROM": "chr1", "ALT": "G"},
                    "multiallelic": {"CHROM": "chr1", "ALT": "C,G"}}
        self.assertEqual(pair_status("biallelic", "multiallelic", variants), "unsupported_alleles")
        self.assertEqual(pair_status("multiallelic", "biallelic", variants), "unsupported_alleles")

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
