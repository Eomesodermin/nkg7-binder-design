#!/usr/bin/env python
"""
harvest.py — collect Boltz-2 validation outputs across all arms of a campaign,
run the filter, and produce a single ranked candidate table + per-arm summary.

Runs after submit_campaign.sh finishes (or on whatever has completed so far).
Reusable: reads the same campaign config; nothing target-specific.

Usage:
  python harvest.py --config configs/nkg7_campaign.yaml --mode smoketest \
      --ws /lustre/.../nkg7_binder_design --out harvest_smoketest
"""
import argparse, os, glob, subprocess, sys, csv, json, yaml


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--mode", default="smoketest")
    ap.add_argument("--ws", required=True)
    ap.add_argument("--out", default="harvest")
    ap.add_argument("--filter-script", default=os.path.join(os.path.dirname(__file__), "filter_designs.py"))
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    f = cfg["defaults"]["filters"]
    os.makedirs(args.out, exist_ok=True)
    arms = list(cfg["arms"].keys())

    combined = []
    per_arm = []
    for arm in arms:
        a = cfg["arms"][arm]
        boltz_dir = f"{args.ws}/runs/{args.mode}/{arm}/03_boltz"
        bb_dir = f"{args.ws}/runs/{args.mode}/{arm}/01_rfdiffusion"
        ranked = os.path.join(args.out, f"{arm}_ranked.csv")
        nconf = len(glob.glob(os.path.join(boltz_dir, "**", "confidence_*_model_0.json"), recursive=True))
        if nconf == 0:
            per_arm.append(dict(arm=arm, n_scored=0, n_pass=0, best_iptm="", best_score="", status="no boltz output"))
            continue
        cmd = [sys.executable, args.filter_script,
               "--boltz-dir", boltz_dir, "--backbone-dir", bb_dir,
               "--binder-chain", a["binder_chain"], "--target-chain", a["target_chain"],
               "--iptm-min", str(f["iptm_min"]), "--plddt-min", str(f["complex_plddt_min"]),
               "--binder-ptm-min", str(f["binder_ptm_min"]), "--rmsd-max", str(f["binder_rmsd_max"]),
               "--arm", arm, "--out", ranked]
        subprocess.run(cmd, check=True)
        rows = list(csv.DictReader(open(ranked)))
        for r in rows:
            r["species"] = a["species"]; r["epitope"] = a["epitope"]
            combined.append(r)
        npass = sum(1 for r in rows if r["pass"] == "True")
        best = rows[0] if rows else {}
        per_arm.append(dict(arm=arm, n_scored=len(rows), n_pass=npass,
                            best_iptm=best.get("iptm", ""), best_score=best.get("composite_score", ""),
                            status="ok"))

    # combined ranked table (passers first, then by composite score)
    combined.sort(key=lambda r: (r.get("pass") != "True", -float(r.get("composite_score") or 0)))
    comb_path = os.path.join(args.out, "all_arms_ranked.csv")
    if combined:
        cols = ["species", "epitope", "arm", "design", "pass", "composite_score",
                "iptm", "complex_plddt", "binder_ptm", "binder_rmsd", "confidence_score", "pdb"]
        with open(comb_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader(); w.writerows(combined)

    sum_path = os.path.join(args.out, "campaign_summary.csv")
    with open(sum_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["arm", "n_scored", "n_pass", "best_iptm", "best_score", "status"])
        w.writeheader(); w.writerows(per_arm)

    print(f"[harvest] {len(combined)} designs across {len(arms)} arms")
    print(f"[harvest] per-arm summary -> {sum_path}")
    print(f"[harvest] combined ranked -> {comb_path}")
    for r in per_arm:
        print(f"  {r['arm']:18s} scored={r['n_scored']:3} pass={r['n_pass']:3} best_iptm={r['best_iptm']} ({r['status']})")


if __name__ == "__main__":
    main()
