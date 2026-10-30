#!/usr/bin/env python3
"""Sequence-complexity audit of designed binders.

Interface-confidence scores (ipTM) from single-sequence co-folding can be inflated by
low-complexity, alanine-rich helical sequences, a known failure mode of ProteinMPNN at low
sampling temperature. This script measures that for a ranked design table and reports how many
confident designs survive a complexity filter.

Per design it adds:
    ala_fraction         fraction of alanine residues
    shannon_entropy      Shannon entropy of the amino-acid composition (bits; 4.32 is the maximum)
    longest_homopolymer  longest run of one repeated residue
    low_complexity       True if any threshold below is violated

Thresholds are heuristics, not calibrated cut-offs. Typical natural proteins have entropy ~4.0-4.2
bits and alanine ~8%; the defaults are deliberately permissive for short (50-65 aa) designs.

Usage:
    python sequence_complexity.py --designs results/production/production_all_arms_ranked.csv \
        --out-dir results/production [--plot figures/sequence_complexity.png]
"""
import argparse
import csv
import math
import os
import re
from collections import Counter

MAX_ALA_FRACTION = 0.30
MIN_ENTROPY_BITS = 3.2
MAX_HOMOPOLYMER = 5
IPTM_MIN = 0.5
PLDDT_MIN = 0.7


def ala_fraction(seq):
    """Fraction of residues that are alanine."""
    return seq.count("A") / len(seq)


def shannon_entropy(seq):
    """Shannon entropy (bits) of the residue composition."""
    n = len(seq)
    return -sum(c / n * math.log2(c / n) for c in Counter(seq).values())


def longest_homopolymer(seq):
    """Length of the longest run of a single repeated residue."""
    return max(len(m.group(0)) for m in re.finditer(r"(.)\1*", seq))


def is_low_complexity(seq):
    """True if the sequence violates the alanine, entropy or homopolymer threshold."""
    return (
        ala_fraction(seq) > MAX_ALA_FRACTION
        or shannon_entropy(seq) < MIN_ENTROPY_BITS
        or longest_homopolymer(seq) > MAX_HOMOPOLYMER
    )


def annotate(rows):
    """Return one dict per design with complexity metrics and the key confidence scores."""
    out = []
    for row in rows:
        seq = row["binder_seq"]
        out.append(
            {
                "arm": row["arm"],
                "design": row["design"],
                "iptm": float(row["iptm"]),
                "complex_plddt": float(row["complex_plddt"]),
                "passes_boltz_filters": row["pass"] == "True",
                "ala_fraction": round(ala_fraction(seq), 3),
                "shannon_entropy": round(shannon_entropy(seq), 3),
                "longest_homopolymer": longest_homopolymer(seq),
                "low_complexity": is_low_complexity(seq),
            }
        )
    return out


def median(values):
    values = sorted(values)
    mid = len(values) // 2
    return values[mid] if len(values) % 2 else (values[mid - 1] + values[mid]) / 2


def summarise(designs):
    """Funnel counts with and without the complexity filter, plus alanine medians."""
    confident = [d for d in designs if d["iptm"] >= IPTM_MIN]
    folded = [d for d in confident if d["complex_plddt"] >= PLDDT_MIN]
    diverse = [d for d in designs if not d["low_complexity"]]
    return [
        ("designs scored", len(designs)),
        ("ipTM >= 0.5", len(confident)),
        ("ipTM >= 0.5 and complex pLDDT >= 0.7", len(folded)),
        ("passing all Boltz-2 filters", sum(d["passes_boltz_filters"] for d in designs)),
        ("not low-complexity", len(diverse)),
        ("not low-complexity, ipTM >= 0.5", sum(d["iptm"] >= IPTM_MIN for d in diverse)),
        (
            "not low-complexity, ipTM >= 0.5 and pLDDT >= 0.7",
            sum(d["iptm"] >= IPTM_MIN and d["complex_plddt"] >= PLDDT_MIN for d in diverse),
        ),
        ("median alanine fraction, all designs", round(median([d["ala_fraction"] for d in designs]), 3)),
        ("median alanine fraction, ipTM >= 0.5", round(median([d["ala_fraction"] for d in confident]), 3)),
    ]


def write_csv(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot(designs, path):
    """Alanine fraction against ipTM, coloured by the complexity flag."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.2))
    for flag, colour, label in ((False, "#2b6cb0", "not low-complexity"), (True, "#c05050", "low-complexity")):
        subset = [d for d in designs if d["low_complexity"] is flag]
        left.scatter(
            [d["ala_fraction"] for d in subset], [d["iptm"] for d in subset], s=6, alpha=0.4, c=colour, label=label
        )
    left.axhline(IPTM_MIN, ls="--", c="grey", lw=1)
    left.axvline(MAX_ALA_FRACTION, ls=":", c="grey", lw=1)
    left.set_xlabel("Alanine fraction of binder")
    left.set_ylabel("Boltz-2 ipTM")
    left.set_title("Confident designs are alanine-rich")
    left.legend(frameon=False, markerscale=2)

    confident = [d["ala_fraction"] for d in designs if d["iptm"] >= IPTM_MIN]
    right.hist([d["ala_fraction"] for d in designs], bins=30, alpha=0.6, label="all designs", color="#888888")
    right.hist(confident, bins=30, alpha=0.7, label="ipTM >= 0.5", color="#c05050")
    right.axvline(0.08, c="k", lw=1, ls=":")
    right.text(0.09, right.get_ylim()[1] * 0.9, "typical protein ~8%", fontsize=8)
    right.set_xlabel("Alanine fraction of binder")
    right.set_ylabel("Designs")
    right.legend(frameon=False)
    for ax in (left, right):
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fig.savefig(path, dpi=160)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--designs", required=True, help="ranked design CSV with binder_seq, iptm, complex_plddt")
    parser.add_argument("--out-dir", default=".", help="where to write the CSV outputs")
    parser.add_argument("--plot", help="optional PNG path (needs matplotlib)")
    args = parser.parse_args()

    with open(args.designs) as handle:
        designs = annotate(csv.DictReader(handle))
    os.makedirs(args.out_dir, exist_ok=True)

    write_csv(os.path.join(args.out_dir, "sequence_complexity.csv"), list(designs[0]), designs)
    summary = summarise(designs)
    write_csv(
        os.path.join(args.out_dir, "complexity_funnel.csv"),
        ["step", "value"],
        [{"step": step, "value": value} for step, value in summary],
    )
    for step, value in summary:
        print(f"{step:55s} {value}")
    if args.plot:
        plot(designs, args.plot)
        print("wrote", args.plot)


if __name__ == "__main__":
    main()
