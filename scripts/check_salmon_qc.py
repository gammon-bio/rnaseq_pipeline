#!/usr/bin/env python3
"""
Check Salmon QC metrics and flag samples below thresholds.
Run after salmon_pipeline.sh, before run_deseq2.R

Usage: python scripts/check_salmon_qc.py --salmon_dir out/salmon --min_mapping 0.6
"""

import argparse
import json
from pathlib import Path
import sys


def check_salmon_qc(salmon_dir: Path, min_mapping: float = 0.6) -> dict:
    """Check mapping rates for all samples in salmon output directory."""
    results = {"passed": [], "failed": [], "warnings": []}

    sample_dirs = [d for d in salmon_dir.iterdir() if d.is_dir()]

    if not sample_dirs:
        print(f"ERROR: No sample directories found in {salmon_dir}", file=sys.stderr)
        sys.exit(1)

    for sample_dir in sorted(sample_dirs):
        meta_file = sample_dir / "aux_info" / "meta_info.json"

        if not meta_file.exists():
            results["warnings"].append(f"{sample_dir.name}: meta_info.json not found")
            continue

        with open(meta_file) as f:
            meta = json.load(f)

        # Salmon stores mapping rate as a fraction (0-1) or percentage
        mapping_rate = meta.get("percent_mapped", meta.get("mapping_rate", None))

        # Some Salmon versions use different keys
        if mapping_rate is None:
            num_mapped = meta.get("num_mapped", 0)
            num_processed = meta.get("num_processed", 1)
            mapping_rate = num_mapped / num_processed if num_processed > 0 else 0

        # Convert to fraction if stored as percentage
        if mapping_rate > 1:
            mapping_rate = mapping_rate / 100

        sample_result = {
            "sample": sample_dir.name,
            "mapping_rate": mapping_rate,
            "num_processed": meta.get("num_processed", "N/A")
        }

        if mapping_rate < min_mapping:
            results["failed"].append(sample_result)
        else:
            results["passed"].append(sample_result)

    return results


def main():
    parser = argparse.ArgumentParser(description="Check Salmon QC metrics")
    parser.add_argument("--salmon_dir", type=Path, default=Path("out/salmon"),
                        help="Path to Salmon output directory (default: out/salmon)")
    parser.add_argument("--min_mapping", type=float, default=0.6,
                        help="Minimum mapping rate threshold (default: 0.6 = 60%%)")
    parser.add_argument("--out", type=Path, default=None,
                        help="Optional: write summary to file")
    args = parser.parse_args()

    if not args.salmon_dir.exists():
        print(f"ERROR: Salmon directory not found: {args.salmon_dir}", file=sys.stderr)
        sys.exit(1)

    results = check_salmon_qc(args.salmon_dir, args.min_mapping)

    # Build summary
    lines = []
    lines.append("=" * 60)
    lines.append("SALMON QC SUMMARY")
    lines.append("=" * 60)
    lines.append(f"Threshold: {args.min_mapping:.0%} mapping rate")
    lines.append(f"Samples checked: {len(results['passed']) + len(results['failed'])}")
    lines.append("")

    if results["passed"]:
        lines.append(f"PASSED ({len(results['passed'])} samples):")
        for s in results["passed"]:
            lines.append(f"  {s['sample']}: {s['mapping_rate']:.1%}")

    if results["failed"]:
        lines.append("")
        lines.append(f"FAILED ({len(results['failed'])} samples):")
        for s in results["failed"]:
            lines.append(f"  {s['sample']}: {s['mapping_rate']:.1%} [BELOW THRESHOLD]")

    if results["warnings"]:
        lines.append("")
        lines.append(f"WARNINGS ({len(results['warnings'])}):")
        for w in results["warnings"]:
            lines.append(f"  {w}")

    lines.append("")
    lines.append("=" * 60)

    summary = "\n".join(lines)
    print(summary)

    # Write to file if requested
    if args.out:
        args.out.write_text(summary)
        print(f"\nSummary written to: {args.out}")

    # Exit with error code if any samples failed
    if results["failed"]:
        print("\nACTION REQUIRED: Remove failed samples from sample_table.csv before running DESeq2")
        sys.exit(1)

    sys.exit(0)


if __name__ == "__main__":
    main()
