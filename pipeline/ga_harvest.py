#!/usr/bin/env python
"""
ga_harvest.py — collate the partial-diffusion GA refinement outputs.

For each arm reads runs/refine_ga/<arm>/hall.csv (the global-best designs across
all generations, seeds included) and:
  - writes a per-arm ranked hall,
  - reports whether the GA IMPROVED on the seed best (gen 0),
  - assembles a campaign summary (best composite/iptm/plddt per arm, seed vs final),
  - copies each arm's top-N refined complex PDBs into ga_top_candidates/<arm>/.

Stdlib only. Run in any env with python3.
"""
import argparse, csv, os, glob, shutil

def read_csv(p):
    if not os.path.exists(p): return []
    with open(p) as fh: return list(csv.DictReader(fh))

def fnum(x, d=0.0):
    try: return float(x)
    except (TypeError, ValueError): return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ga-dir", required=True)
    ap.add_argument("--arms", required=True, help="space-separated arm names")
    ap.add_argument("--out", required=True)
    ap.add_argument("--topn", type=int, default=10)
    args = ap.parse_args()
    arms = args.arms.split()

    summary = []
    all_rows = []
    for arm in arms:
        hall = read_csv(os.path.join(args.ga_dir, arm, "hall.csv"))
        if not hall:
            summary.append(dict(arm=arm, status="NO_HALL", n_designs=0)); continue
        hall = sorted(hall, key=lambda r: -fnum(r.get("composite_score")))
        # seed best = best gen-0 row; final best = best overall
        seed_rows = [r for r in hall if str(r.get("gen")) == "0"]
        seed_best = max((fnum(r["composite_score"]) for r in seed_rows), default=0.0)
        final_best_row = hall[0]
        final_best = fnum(final_best_row["composite_score"])
        improved = final_best > seed_best + 1e-6
        gain = round(final_best - seed_best, 4)
        # per-arm ranked hall
        arm_out = os.path.join(args.ga_dir, arm, f"{arm}_ga_ranked.csv")
        cols = list(hall[0].keys())
        with open(arm_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(hall)
        # winners bundle
        topdir = os.path.join(args.ga_dir, "ga_top_candidates", arm)
        os.makedirs(topdir, exist_ok=True)
        for r in hall[:args.topn]:
            p = r.get("pdb", "")
            if p and os.path.exists(p):
                shutil.copyfile(p, os.path.join(topdir, f"{r['design']}.pdb"))
        summary.append(dict(arm=arm, status="ok", n_designs=len(hall),
            seed_best_composite=round(seed_best, 4), final_best_composite=round(final_best, 4),
            composite_gain=gain, improved=int(improved),
            final_best_design=final_best_row["design"],
            final_best_iptm=final_best_row.get("iptm", ""),
            final_best_plddt=final_best_row.get("complex_plddt", ""),
            final_best_gen=final_best_row.get("gen", "")))
        for r in hall[:args.topn]:
            r2 = dict(r); r2["arm"] = arm; all_rows.append(r2)

    # write campaign summary
    scols = ["arm", "status", "n_designs", "seed_best_composite", "final_best_composite",
             "composite_gain", "improved", "final_best_design", "final_best_iptm",
             "final_best_plddt", "final_best_gen"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=scols); w.writeheader()
        for s in summary: w.writerow({k: s.get(k, "") for k in scols})
    # combined top table
    if all_rows:
        comb = os.path.join(os.path.dirname(args.out), "ga_all_top.csv")
        cols = sorted({k for r in all_rows for k in r})
        with open(comb, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(all_rows)

    print("=== GA REFINEMENT SUMMARY ===")
    for s in summary:
        if s["status"] == "ok":
            print(f"  {s['arm']:14s} seed={s['seed_best_composite']:.4f} -> "
                  f"final={s['final_best_composite']:.4f} (gain {s['composite_gain']:+.4f}, "
                  f"improved={'YES' if s['improved'] else 'no'}, gen {s['final_best_gen']})")
        else:
            print(f"  {s['arm']:14s} {s['status']}")
    print(f"[ga_harvest] summary -> {args.out}")

if __name__ == "__main__":
    main()
