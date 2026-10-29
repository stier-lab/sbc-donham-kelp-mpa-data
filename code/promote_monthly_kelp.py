#!/usr/bin/env python3
"""Promote only verified kelp changes from a complete staged harmonization."""

import argparse
import csv
import hashlib
import io
import json
import math
import tempfile
from pathlib import Path


TAXON = "Macrocystis pyrifera"
FILES = ("harmonized_raw_responses.csv", "harmonized_response_ratios.csv")


def md5(data):
    return hashlib.md5(data).hexdigest()


def prepare(baseline, candidate, raw, *, taxon_name=TAXON,
            sources=("KFM", "LTER", "PISCO"), responses=None):
    if b"\r" in baseline or not baseline.endswith(b"\n"):
        raise ValueError("Expected one LF-terminated line per canonical record")
    lines = baseline.decode().splitlines(keepends=True)
    before = list(csv.DictReader(io.StringIO(baseline.decode())))
    after = list(csv.DictReader(io.StringIO(candidate.decode())))
    if not before or len(lines) != len(before) + 1 or len(before) != len(after):
        raise ValueError("Record count or physical-line count changed")
    columns = list(before[0])
    if list(after[0]) != columns:
        raise ValueError("Schema changed")
    taxon = "taxon_name" if raw else "y"
    taxa = (taxon_name,) if isinstance(taxon_name, str) else tuple(taxon_name)
    keys = ("CA_MPA_Name_Short", "year", taxon, "source", "resp") + (("status",) if raw else ())
    fields = ("value",) if raw else ("mpa", "reference", "Diff", "lnDiff")
    key = lambda r: tuple(r[k] for k in keys)
    old = {key(r): r for r in before}
    new = {key(r): r for r in after}
    if len(old) != len(before) or len(new) != len(after) or old.keys() != new.keys():
        raise ValueError("Duplicate, added or missing record identity")
    changes = {}
    records = []
    for i, row in enumerate(before):
        proposed = new[key(row)]
        target = (row[taxon] in taxa and row["source"] in sources
                  and (responses is None or row["resp"] in responses))
        for field in columns:
            if field in fields:
                a, b = float(row[field]), float(proposed[field])
                if not math.isfinite(a) or not math.isfinite(b):
                    raise ValueError("Nonfinite response value")
                if not target and not math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-10):
                    raise ValueError("Non-target value changed: " + repr(key(row)))
            elif row[field] != proposed[field]:
                raise ValueError("Record metadata changed: " + repr((key(row), field)))
        changed = target and any(not math.isclose(float(row[f]), float(proposed[f]),
                                                 rel_tol=1e-12, abs_tol=1e-12) for f in fields)
        if not changed:
            continue
        label = row["source"] + "/" + row["resp"]
        changes[label] = changes.get(label, 0) + 1
        updated = dict(row)
        for field in fields:
            updated[field] = float(proposed[field])
        stream = io.StringIO(newline="")
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n", quoting=csv.QUOTE_NONNUMERIC)
        writer.writerow(updated)
        lines[i + 1] = stream.getvalue()
        records.append({"key": key(row), "before": {f: row[f] for f in fields},
                        "after": {f: proposed[f] for f in fields}})
    result = "".join(lines).encode()
    return result, {"rows": len(before), "changed_rows": len(records), "changes": changes,
                    "non_target_rows_preserved": sum(not (r[taxon] in taxa and r["source"] in sources
                        and (responses is None or r["resp"] in responses)) for r in before),
                    "before_md5": md5(baseline), "after_md5": md5(result)}, records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--target", type=Path, action="append", required=True)
    parser.add_argument("--calibration", type=Path, required=True)
    parser.add_argument("--record", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    calibration = args.calibration.read_bytes()
    if md5(calibration) != "50eda5d3664bf70780b3bead9674fac1":
        raise ValueError("Wrong calibration source")
    staged, results, details = {}, {}, {}
    for name in FILES:
        original = (args.baseline / name).read_bytes()
        staged[name], results[name], details[name] = prepare(
            original, (args.candidate / name).read_bytes(), "raw" in name)
        for target in args.target:
            if (target / name).read_bytes() != original:
                raise ValueError("Destination differs from preserved baseline: " + str(target / name))
    report = {"package": "knb-lter-sbc.112.15", "model": 1, "calibration_md5": md5(calibration),
              "mass_basis": "dry foliar", "conversion": "recorded month before annual aggregation",
              "files": results, "applied": args.apply}
    print(json.dumps(report, indent=2))
    if not args.apply:
        return
    temporary = []
    try:
        for target in args.target:
            for name, data in staged.items():
                with tempfile.NamedTemporaryFile(dir=target, delete=False) as f:
                    f.write(data)
                    temporary.append((Path(f.name), target / name))
        for source, target in temporary:
            source.replace(target)
        for target in args.target:
            (target / "kelp_monthly_calibration.json").write_text(json.dumps(report, indent=2) + "\n")
            csvs = sorted(target.glob("*.csv"))
            (target / "checksums.md5").write_text("".join(
                md5(p.read_bytes()) + " " + p.name + "\n" for p in csvs))
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(json.dumps({**report, "changed_records": details}, indent=2) + "\n")
    finally:
        for source, _target in temporary:
            source.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
