#!/usr/bin/env python3
"""Promote only PISCO sheephead biomass from a complete staged harmonization."""

import argparse
import json
import tempfile
from pathlib import Path

from promote_monthly_kelp import FILES, md5, prepare


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--target", type=Path, action="append", required=True)
    parser.add_argument("--record", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    source = json.loads((args.candidate / "pisco_sheephead_allometry_source.json").read_text())
    expected = {"reference_doi": "10.1098/rspb.2014.1817", "a": 0.0144, "b": 3.04,
                "length_type": "TL", "length_unit": "cm", "weight_unit": "g",
                "count_weighted": True, "original_lookup_md5": "6ccfacd18fff07208f510384e914a3b2"}
    if any(source.get(k) != v for k, v in expected.items()):
        raise ValueError("Unverified sheephead source or units")
    staged, results, details = {}, {}, {}
    for name in FILES:
        original = (args.baseline / name).read_bytes()
        staged[name], results[name], details[name] = prepare(
            original, (args.candidate / name).read_bytes(), "raw" in name,
            taxon_name="Semicossyphus pulcher", sources=("PISCO",), responses=("Bio",))
        for target in args.target:
            if target.resolve() == args.baseline.resolve():
                raise ValueError("Baseline cannot be a promotion destination")
            if (target / name).read_bytes() != original:
                raise ValueError("Destination differs from preserved baseline: " + str(target / name))
    for name in ("harmonized_landsat_rr.csv", "harmonized_site_metadata.csv"):
        if (args.baseline / name).read_bytes() != (args.candidate / name).read_bytes():
            raise ValueError("Unrelated export changed: " + name)
    report = {"release": "hamilton_caselle_2015_sheephead", "source": source,
              "scope": {"programme": "PISCO", "taxon": "Semicossyphus pulcher", "response": "Bio"},
              "files": results, "applied": args.apply}
    print(json.dumps(report, indent=2))
    if not args.apply:
        return
    temporary = []
    try:
        for target in args.target:
            for name, data in staged.items():
                with tempfile.NamedTemporaryFile(dir=target, delete=False) as handle:
                    handle.write(data)
                    temporary.append((Path(handle.name), target / name))
        for old, new in temporary:
            old.replace(new)
        for target in args.target:
            (target / "sheephead_allometry.json").write_text(json.dumps(report, indent=2) + "\n")
            (target / "checksums.md5").write_text("".join(
                md5(p.read_bytes()) + " " + p.name + "\n" for p in sorted(target.glob("*.csv"))))
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(json.dumps({**report, "changed_records": details}, indent=2) + "\n")
    finally:
        for path, _ in temporary:
            path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
