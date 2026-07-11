#!/usr/bin/env python
"""
filter_designs.py — parse Boltz-2 outputs for binder:target complexes,
compute interface metrics, and rank designs.

Metrics (standard de novo binder filters, cf. Bennett et al. 2023, adapted to
Boltz-2 confidence fields which are 0-1 scaled):
  - iptm          : interface pTM (complex/interface confidence; >0.5 pass line)
  - complex_plddt : whole-complex pLDDT 0-1 (fold confidence; >0.7 pass)
  - binder_ptm    : per-chain pTM for the binder chain (fold quality of binder)
  - binder_rmsd   : CA-RMSD between the RFdiffusion backbone and the Boltz-refolded
                    binder (self-consistency; requires the design backbone).

A design PASSES if it clears all configured thresholds. Ranking is by a
composite score (see composite_score) among passers.

Reusable: chain identities and thresholds are CLI args.

Boltz output layout (per input <name>):
  <boltz_out>/boltz_results_<name>/predictions/<name>/
      confidence_<name>_model_0.json   (iptm, complex_plddt, ptm, chains_ptm, ...)
      <name>_model_0.pdb  (or .cif)

Usage:
  python filter_designs.py \
      --boltz-dir 03_boltz --backbone-dir 01_rfdiffusion \
      --binder-chain B --target-chain A \
      --iptm-min 0.5 --plddt-min 0.7 --binder-ptm-min 0.5 --rmsd-max 2.0 \
      --out 04_filtered/ranked.csv
"""
import argparse, json, os, glob, re
import numpy as np


def ca_coords(pdb, chain):
    xs = []
    with open(pdb) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] == chain:
                xs.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(xs)


def ca_rmsd(pdb1, chain1, pdb2, chain2):
    """CA-RMSD (Kabsch) between the binder chain in the design backbone and the
    Boltz-refolded complex. Same sequence -> same CA count expected."""
    A, B = ca_coords(pdb1, chain1), ca_coords(pdb2, chain2)
    if len(A) == 0 or len(A) != len(B):
        return np.nan
    A = A - A.mean(0); B = B - B.mean(0)
    U, S, Vt = np.linalg.svd(A.T @ B)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    A2 = (R @ A.T).T
    return float(np.sqrt(((A2 - B) ** 2).sum(1).mean()))


def binder_ptm(conf, binder_chain_idx):
    cp = conf.get("chains_ptm", {})
    return cp.get(str(binder_chain_idx), None)


def composite_score(row):
    """Higher is better. Boltz fields are 0-1 already."""
    iptm = row["iptm"] or 0
    plddt = row["complex_plddt"] or 0
    bptm = row["binder_ptm"] or 0
    rmsd = row["binder_rmsd"]
    rmsd_term = max(0, 1 - (rmsd / 5.0)) if rmsd == rmsd else 0
    return round(0.45 * iptm + 0.30 * plddt + 0.15 * bptm + 0.10 * rmsd_term, 4)


def find_conf_and_pdb(boltz_dir, name):
    """Locate the Boltz confidence json + model_0 structure for one design."""
    base = os.path.join(boltz_dir, f"boltz_results_{name}", "predictions", name)
    conf = glob.glob(os.path.join(base, "confidence_*_model_0.json"))
    if not conf:  # some layouts drop the boltz_results_ prefix
        conf = glob.glob(os.path.join(boltz_dir, "**", f"confidence_{name}_model_0.json"), recursive=True)
    pdb = glob.glob(os.path.join(base, f"{name}_model_0.pdb")) or           glob.glob(os.path.join(base, f"{name}_model_0.cif"))
    return (conf[0] if conf else None), (pdb[0] if pdb else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boltz-dir", required=True)
    ap.add_argument("--backbone-dir", default=None,
                    help="RFdiffusion backbone dir for self-consistency RMSD (optional)")
    ap.add_argument("--binder-chain", default="B")
    ap.add_argument("--target-chain", default="A")
    ap.add_argument("--binder-chain-idx", type=int, default=1,
                    help="0-based chain index of the binder in the YAML order (target=0, binder=1)")
    ap.add_argument("--iptm-min", type=float, default=0.5)
    ap.add_argument("--plddt-min", type=float, default=0.7)
    ap.add_argument("--binder-ptm-min", type=float, default=0.5)
    ap.add_argument("--rmsd-max", type=float, default=2.0)
    ap.add_argument("--arm", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # enumerate designs by their confidence files
    conf_files = glob.glob(os.path.join(args.boltz_dir, "**", "confidence_*_model_0.json"), recursive=True)
    rows = []
    for cf in sorted(conf_files):
        m = re.search(r"confidence_(.+)_model_0\.json", os.path.basename(cf))
        name = m.group(1) if m else os.path.basename(cf)
        conf = json.load(open(cf))
        pdb = None
        pdir = os.path.dirname(cf)
        for ext in ("pdb", "cif"):
            hits = glob.glob(os.path.join(pdir, f"{name}_model_0.{ext}"))
            if hits:
                pdb = hits[0]; break
        rmsd = np.nan
        if args.backbone_dir and pdb and pdb.endswith(".pdb"):
            bb = glob.glob(os.path.join(args.backbone_dir, f"{name}.pdb"))
            if not bb:
                # design name is <arm>_<backbone#>_s<sample>; the backbone file is
                # <arm>_<backbone#>.pdb -> strip the trailing _s<sample> suffix.
                stem = re.sub(r"_s\d+$", "", name)
                bb = glob.glob(os.path.join(args.backbone_dir, f"{stem}.pdb"))
                if not bb:  # fallback: strip last underscore token
                    bb = glob.glob(os.path.join(args.backbone_dir, f"{stem.rsplit('_',1)[0]}.pdb"))
            if bb:
                rmsd = ca_rmsd(bb[0], args.binder_chain, pdb, args.binder_chain)
        iptm = conf.get("iptm", conf.get("protein_iptm"))
        plddt = conf.get("complex_plddt")
        bptm = binder_ptm(conf, args.binder_chain_idx)
        row = dict(
            arm=args.arm, design=name,
            iptm=round(iptm, 4) if iptm is not None else np.nan,
            complex_plddt=round(plddt, 4) if plddt is not None else np.nan,
            binder_ptm=round(bptm, 4) if bptm is not None else np.nan,
            binder_rmsd=round(rmsd, 2) if rmsd == rmsd else np.nan,
            confidence_score=conf.get("confidence_score"),
            pdb=os.path.basename(pdb) if pdb else "",
        )
        row["pass"] = bool(
            (row["iptm"] >= args.iptm_min if row["iptm"] == row["iptm"] else False)
            and (row["complex_plddt"] >= args.plddt_min if row["complex_plddt"] == row["complex_plddt"] else False)
            and (row["binder_ptm"] >= args.binder_ptm_min if row["binder_ptm"] == row["binder_ptm"] else True)
            and (row["binder_rmsd"] <= args.rmsd_max if row["binder_rmsd"] == row["binder_rmsd"] else True)
        )
        row["composite_score"] = composite_score(row)
        rows.append(row)

    rows.sort(key=lambda r: (-int(r["pass"]), -r["composite_score"]))
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    import csv
    cols = ["arm", "design", "pass", "composite_score", "iptm", "complex_plddt",
            "binder_ptm", "binder_rmsd", "confidence_score", "pdb"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    npass = sum(r["pass"] for r in rows)
    print(f"[filter] {len(rows)} designs, {npass} pass -> {args.out}")


if __name__ == "__main__":
    main()
