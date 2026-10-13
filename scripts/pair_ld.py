#!/usr/bin/env python3
"""Calculate unphased dosage r2 for requested pairs and sample groups."""
import argparse
import csv
import gzip
import io
import logging
import math
import os
from pathlib import Path
import subprocess
import tarfile

LOG = logging.getLogger("pair_ld")


def require_local(value):
    if "://" in str(value):
        raise ValueError(f"File localization error: unresolved URI {value}")
    path = Path(value)
    if not path.is_file() or not os.access(path, os.R_OK) or path.stat().st_size == 0:
        raise ValueError(f"File localization error: input is not readable or is empty: {value}")
    return path


def read_tsv(path, required):
    with path.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or not set(required).issubset(reader.fieldnames):
            raise ValueError(f"{path}: required TSV columns: {', '.join(required)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path}: table is empty")
    for row in rows:
        if None in row or any(not row.get(key) for key in required):
            raise ValueError(f"{path}: row has missing values or extra fields")
    return rows


def token(value):
    if not value or value == "." or any(c.isspace() or ord(c) < 32 for c in value):
        raise ValueError(f"Invalid sample or variant ID: {value!r}")
    return value


def read_groups(psam, populations):
    samples = []
    with psam.open() as stream:
        header = None
        for line in stream:
            if line.startswith("##") or not line.strip():
                continue
            fields = line.split()
            if header is None:
                header = fields
                header[0] = header[0].lstrip("#")
                if "IID" not in header:
                    raise ValueError("PSAM must have a header with an IID column")
                continue
            row = dict(zip(header, fields))
            if len(fields) != len(header):
                raise ValueError("PSAM row does not match the header")
            samples.append((token(row.get("FID", "0")), token(row["IID"])))
    if not samples or len(samples) != len(set(samples)):
        raise ValueError("PSAM must contain unique FID/IID sample IDs; duplicate IDs with SID are unsupported")
    labels = {}
    for row in read_tsv(populations, ["FID", "IID", "population"]):
        sample = (token(row["FID"]), token(row["IID"]))
        if sample in labels:
            raise ValueError(f"Duplicate sample in population table: {sample}")
        label = row["population"]
        if label == "ALL" or any(ord(c) < 32 for c in label):
            raise ValueError("Population label must be printable and must not be ALL")
        labels[sample] = label
    if set(samples) != set(labels):
        raise ValueError("Each PSAM sample must have one population; the population table must contain no extra samples")
    groups = {"ALL": samples}
    for label in sorted(set(labels.values())):
        groups[label] = [s for s in samples if labels[s] == label]
    return groups


def read_variants(pvar, requested):
    if pvar.suffix == ".zst":
        import zstandard
        stream = io.TextIOWrapper(zstandard.ZstdDecompressor().stream_reader(pvar.open("rb")))
    elif pvar.suffix == ".gz":
        stream = gzip.open(pvar, "rt")
    else:
        stream = pvar.open()
    variants = {}
    with stream:
        header = None
        for line in stream:
            if line.startswith("##") or not line.strip():
                continue
            fields = line.split()
            if header is None:
                header = fields
                header[0] = header[0].lstrip("#")
                if not {"CHROM", "POS", "ID", "REF", "ALT"}.issubset(header):
                    raise ValueError("PVAR must have a header with CHROM, POS, ID, REF, and ALT")
                continue
            if len(fields) != len(header):
                raise ValueError("PVAR row does not match the header")
            row = dict(zip(header, fields))
            if row["ID"] in requested:
                if row["ID"] in variants:
                    raise ValueError(f"Requested variant ID is not unique: {row['ID']}")
                variants[row["ID"]] = row
    if header is None:
        raise ValueError("PVAR has no header")
    return variants


def pair_status(a, b, variants):
    if a not in variants or b not in variants:
        return "missing_variant"
    for variant in (variants[a], variants[b]):
        if variant["CHROM"].removeprefix("chr") not in {str(i) for i in range(1, 23)}:
            return "unsupported_chromosome"
        if "," in variant["ALT"] or variant["ALT"] == ".":
            return "unsupported_alleles"
    return "ok"


def run_plink(executable, arguments):
    LOG.info("Run PLINK: %s", " ".join(map(str, arguments)))
    subprocess.run([executable, *map(str, arguments)], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("pgen", "pvar", "psam", "pairs", "populations"):
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument("--pgi")
    parser.add_argument("--plink", default="plink2")
    parser.add_argument("--out-dir", default="results")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--memory-mb", type=int, default=6000)
    parser.add_argument("--allow-small-groups", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    if args.threads < 1 or args.memory_mb < 640:
        parser.error("threads must be positive and memory-mb must be at least 640")
    paths = {name: require_local(getattr(args, name)) for name in ("pgen", "pvar", "psam", "pairs", "populations")}
    pgi = require_local(args.pgi) if args.pgi else None
    pairs = [(token(r["variant_a"]), token(r["variant_b"])) for r in read_tsv(paths["pairs"], ["variant_a", "variant_b"])]
    if any(a == b for a, b in pairs):
        raise ValueError("A variant pair must contain two different variant IDs")
    groups = read_groups(paths["psam"], paths["populations"])
    variants = read_variants(paths["pvar"], {v for pair in pairs for v in pair})
    eligible = {(min(a, b), max(a, b)) for a, b in pairs if pair_status(a, b, variants) == "ok"}
    output = Path(args.out_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    work = output / "plink"
    work.mkdir()
    # Omitting parent columns makes every input sample eligible for LD.
    # Keep the original sample order, including samples with known parents.
    founders = work / "all_samples.psam"
    founders.write_text("#FID\tIID\n" + "".join(f"{fid}\t{iid}\n" for fid, iid in groups["ALL"]))
    version = subprocess.check_output([args.plink, "--version"], text=True)
    (output / "plink_version.txt").write_text(version)
    resources = ["--threads", args.threads, "--memory", args.memory_mb]
    if eligible:
        extract = work / "requested.txt"
        extract.write_text("".join(v + "\n" for v in sorted({v for pair in eligible for v in pair})))
        inputs = ["--pgen", paths["pgen"], "--pvar", paths["pvar"], "--psam", founders]
        if pgi:
            inputs += ["--pgi", pgi]
        run_plink(args.plink, inputs + ["--extract", extract, "--make-pgen", "--out", work / "subset"] + resources)
    anchors = {}
    for a, b in sorted(eligible):
        anchors.setdefault(a, []).append(b)
    fields = ["pair_index", "variant_a", "variant_b", "population", "n_samples", "r2", "status", "small_group"]
    with (output / "ld.tsv").open("w") as stream, (output / "sample_counts.tsv").open("w") as counts:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        count_writer = csv.writer(counts, delimiter="\t", lineterminator="\n")
        count_writer.writerow(["population", "n_samples", "ld_enabled"])
        for group_index, (label, samples) in enumerate(groups.items()):
            n = len(samples)
            enabled = n >= (3 if args.allow_small_groups else 50)
            LOG.info("Population %s: %s samples; LD enabled: %s", label, n, enabled)
            count_writer.writerow([label, n, str(enabled).lower()])
            values = {}
            if enabled and eligible:
                keep = work / f"group_{group_index}.keep"
                keep.write_text("#FID\tIID\n" + "".join(f"{fid}\t{iid}\n" for fid, iid in samples))
                for anchor_index, (anchor, partners) in enumerate(anchors.items()):
                    extract = work / f"anchor_{anchor_index}.txt"
                    extract.write_text("\n".join([anchor, *partners]) + "\n")
                    prefix = work / f"group_{group_index}_anchor_{anchor_index}"
                    arguments = ["--pfile", work / "subset", "--keep", keep, "--extract", extract, "--r2-unphased", "inter-chr", "cols=chrom,pos,id", "--ld-snp", anchor, "--ld-window-r2", "-1", "--out", prefix] + resources
                    if args.allow_small_groups and n < 50:
                        arguments.append("--bad-ld")
                    run_plink(args.plink, arguments)
                    with prefix.with_suffix(".vcor").open() as report:
                        for row in csv.DictReader(report, delimiter="\t"):
                            values[tuple(sorted((row["ID_A"], row["ID_B"])))] = float(row["UNPHASED_R2"])
            for index, (a, b) in enumerate(pairs, 1):
                status = pair_status(a, b, variants)
                value = "NA"
                if status == "ok":
                    if not enabled:
                        status = "too_few_samples"
                    else:
                        key = tuple(sorted((a, b)))
                        if key not in values:
                            raise ValueError(f"PLINK did not report requested pair {a}, {b} in {label}")
                        if math.isfinite(values[key]):
                            value = format(values[key], ".10g")
                        else:
                            status = "undefined_ld"
                writer.writerow(dict(zip(fields, [index, a, b, label, n, value, status, str(n < 50).lower()])))
    # Export logs only. Temporary genotype and sample files stay in task storage.
    with tarfile.open(output / "plink_logs.tar.gz", "w:gz") as archive:
        for path in sorted(work.glob("*.log")):
            archive.add(path, arcname=path.name)
    LOG.info("Complete: %s requested pairs across %s groups", len(pairs), len(groups))


if __name__ == "__main__":
    main()
